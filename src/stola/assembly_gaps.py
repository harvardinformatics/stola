"""Query-genome assembly-gap (contig boundary) inference, with no AGP file.

Scaffolded assemblies typically mark the join between two contigs with a run
of N bases. Some assemblies pad every "gap of unknown size" to one fixed
length (a real, but never fully standardized, convention), while others
record estimated gap sizes instead. Without an AGP file, the only evidence
available is the empirical distribution of N-run lengths in the genome file
itself (FASTA or .2bit -- see genome_io.py). Only N-runs of at least
min_gap_length (default 10bp) are considered at all -- anything shorter is
more likely an ambiguous basecall than a gap -- and among those:

- one (or a small handful of) length(s) clearly dominate: a fixed-length
  convention, and only N-runs of those lengths are contig boundaries;
- plenty of N-runs but no dominant length: the assembly records estimated
  gap sizes, so every N-run is a contig boundary;
- otherwise (no N-runs at all, or too few to infer a convention): nothing
  is treated as a contig boundary.
"""

from __future__ import annotations

import re
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from .genome_io import iter_named_sequences

# Sequences are already uppercased by iter_named_sequences.
N_RUN_RE = re.compile(r"N+")

# Below this length, an isolated N is more likely a single low-confidence
# basecall embedded in otherwise-real sequence than a deliberate gap marker;
# such runs are ignored entirely (never a gap, and not counted below).
DEFAULT_MIN_GAP_LENGTH = 10

# A single length must cover at least this fraction of (length-filtered)
# N-runs to call it a confident fixed-length convention outright.
HIGH_CONFIDENCE_SHARE = 0.7
# The top few lengths together must cover at least this fraction to call a
# small candidate set "medium confidence".
MEDIUM_CONFIDENCE_SHARE = 0.8
MEDIUM_CONFIDENCE_TOP_N = 3
# Below this many (length-filtered) N-runs, including none at all, there's
# not enough evidence to infer any convention, and no N-run is a gap.
MIN_RUNS_FOR_INFERENCE = 10


@dataclass(frozen=True)
class NRun:
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class GapConvention:
    """Result of looking for the query assembly's gap-marker convention.

    confidence is one of:
      "high"      one dominant length (a fixed-length gap marker);
      "medium"    a small set of lengths together dominate;
      "estimated" enough qualifying N-runs, but no dominant length -- the
                  assembly records estimated gap sizes, so every qualifying
                  length is trusted as a gap marker;
      "none"      no qualifying N-runs, or too few to infer a convention:
                  nothing is a gap marker.
    Qualifying N-runs are those of at least min_gap_length. boundary_lengths
    is the set of lengths trusted as deliberate gap markers (empty for
    "none"). length_counts is the full, unfiltered length -> count
    histogram, for reporting regardless of the verdict.
    """

    confidence: str
    boundary_lengths: frozenset[int]
    length_counts: Counter


def iter_n_runs(genome_path: str | Path) -> Iterator[NRun]:
    """Yield every maximal run of N in a genome file (FASTA or .2bit).

    Delegates to iter_named_sequences (genome_io.py) for one record's full
    sequence at a time (never holds more than one sequence in memory), so it
    scales to whole-genome multi-GB inputs regardless of format. Only ever
    finds gaps *within* a record (a contig join inside a scaffold/chromosome)
    -- a gap *between* scaffolds is just two separate records and leaves no
    N-run to find.
    """
    for chrom, sequence in iter_named_sequences(genome_path):
        for match in N_RUN_RE.finditer(sequence):
            yield NRun(chrom, match.start(), match.end())


def infer_gap_convention(
    runs: list[NRun], min_gap_length: int = DEFAULT_MIN_GAP_LENGTH
) -> GapConvention:
    # Classify the gap-marker convention from runs of at least
    # min_gap_length; report the full histogram regardless.
    all_length_counts = Counter(run.length for run in runs)
    filtered_counts = Counter(
        {length: count for length, count in all_length_counts.items() if length >= min_gap_length}
    )
    total = sum(filtered_counts.values())

    if total < MIN_RUNS_FOR_INFERENCE:
        return GapConvention("none", frozenset(), all_length_counts)

    ranked = filtered_counts.most_common()
    top_length, top_count = ranked[0]
    if top_count / total >= HIGH_CONFIDENCE_SHARE:
        return GapConvention("high", frozenset({top_length}), all_length_counts)

    top_n = ranked[:MEDIUM_CONFIDENCE_TOP_N]
    if sum(count for _, count in top_n) / total >= MEDIUM_CONFIDENCE_SHARE:
        return GapConvention(
            "medium", frozenset(length for length, _ in top_n), all_length_counts
        )

    # Spread-out lengths: estimated gap sizes, so every qualifying run is a gap.
    return GapConvention("estimated", frozenset(filtered_counts), all_length_counts)


def describe_gap_convention(
    convention: GapConvention, min_gap_length: int = DEFAULT_MIN_GAP_LENGTH
) -> str:
    """One line saying what was inferred and how lost CNEEs get categorized."""
    n = sum(c for length, c in convention.length_counts.items() if length >= min_gap_length)
    runs = f"{n} N-run(s) >= {min_gap_length}bp"
    if convention.confidence in ("high", "medium"):
        lengths = "/".join(f"{length}bp" for length in sorted(convention.boundary_lengths))
        return (f"fixed-length gap marker(s) {lengths} ({convention.confidence} confidence, "
                f"{runs}); lost CNEEs in those N-runs -> assembly_gap_loss, "
                "others -> probable_loss")
    if convention.confidence == "estimated":
        return (f"estimated gap sizes ({runs}, no dominant length); lost CNEEs in any "
                f"N-run >= {min_gap_length}bp -> assembly_gap_loss, others -> probable_loss")
    if n == 0:
        return f"no N-runs >= {min_gap_length}bp (gapless); all lost CNEEs -> probable_loss"
    return (f"only {runs}, too few to infer a gap convention; all lost CNEEs "
            "-> probable_loss")


class NRunIndex:
    """Overlap-only lookup over a set of N-runs, grouped by chromosome."""

    def __init__(self, runs: list[NRun]) -> None:
        by_chrom: dict[str, list[NRun]] = defaultdict(list)
        for run in runs:
            by_chrom[run.chrom].append(run)
        self.by_chrom: dict[str, list[NRun]] = {}
        self.ends_by_chrom: dict[str, list[int]] = {}
        for chrom, chrom_runs in by_chrom.items():
            ordered = sorted(chrom_runs, key=lambda r: (r.start, r.end))
            self.by_chrom[chrom] = ordered
            self.ends_by_chrom[chrom] = [r.end for r in ordered]

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        runs = self.by_chrom.get(chrom)
        if not runs:
            return False
        ends = self.ends_by_chrom[chrom]
        index = bisect_right(ends, start)
        while index < len(runs):
            run = runs[index]
            if run.start >= end:
                break
            if run.end > start:
                return True
            index += 1
        return False


def scan_query_assembly_gaps(
    genome_path: str | Path, min_gap_length: int = DEFAULT_MIN_GAP_LENGTH
) -> tuple[GapConvention, NRunIndex]:
    """Scan a query genome (FASTA or .2bit) and return its gap convention
    plus a boundary index.

    The returned NRunIndex is built only from N-runs whose length is one of
    the convention's trusted boundary_lengths -- empty (matches nothing) when
    confidence is "none", by construction.
    """
    runs = list(iter_n_runs(genome_path))
    convention = infer_gap_convention(runs, min_gap_length)
    boundary_runs = [run for run in runs if run.length in convention.boundary_lengths]
    return convention, NRunIndex(boundary_runs)

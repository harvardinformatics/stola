"""Chain-file QC summaries: chain size and score distributions, and per-base
aligned-block depth across the reference and query genomes.

A chain's size is its aligned bases (the sum of its ungapped block sizes).
Per-base depth at a base is the number of chains with an aligned block over
it -- not the chain's whole [start, end) span, which would count its internal
gaps as covered. Depth covers every base of each genome, so a base no chain
aligns, including all of any sequence no chain touches, has depth 0.

Every chain body is parsed, so this is vectorized throughout: each batch of
chains is located through the chain index's byte offsets and parsed with
pandas' C reader, block coordinates come from per-chain cumulative sums, and
depth is a single sweep over the sorted block starts and ends. Memory grows
with the number of blocks, not with genome size.
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from operator import attrgetter
from pathlib import Path

import numpy as np
import pandas as pd

from .index import ChainIndex
from .models import ChainRecord

# Chains are read and parsed in batches of about this many bytes.
BATCH_BYTES = 64 * 1024 * 1024
# Events handled per step of the depth sweep, to bound its temporaries.
SWEEP_CHUNK = 16 * 1024 * 1024


@dataclass(frozen=True)
class DepthSummary:
    histogram: np.ndarray  # histogram[d] = number of bases at depth d
    total_bases: int
    mean: float
    median: int
    covered_fraction: float  # bases at depth >= 1 / total_bases


@dataclass(frozen=True)
class ChainQc:
    aligned_bases: np.ndarray  # per chain, in chain-file order
    scores: np.ndarray
    ref_depth: DepthSummary
    query_depth: DepthSummary


class _CoordinateSpace:
    """One genome laid end to end, so every base has a single integer coordinate.

    Sequence lengths come from the genome file when given, so sequences no
    chain touches still count (as depth 0). A chain chromosome missing from
    the genome, or with no genome given at all, falls back to the size in its
    chain headers.
    """

    def __init__(
        self,
        chain_sizes: dict[str, int],
        genome_lengths: dict[str, int] | None,
        aliases: dict[str, str] | None,
    ) -> None:
        lengths = dict(genome_lengths or {})
        aliases = aliases or {}
        resolved_by_chain_chrom = {}
        for chrom, size in chain_sizes.items():
            resolved = aliases.get(chrom, chrom)
            if resolved not in lengths:
                lengths[resolved] = size
            resolved_by_chain_chrom[chrom] = resolved
        self.total_bases = sum(lengths.values())

        # Space sequences by at least their chain-header size too, so a block
        # past the end of a shorter genome sequence (a chain built against a
        # different assembly version) can't spill into the next sequence.
        span = dict(lengths)
        for chrom, size in chain_sizes.items():
            resolved = resolved_by_chain_chrom[chrom]
            span[resolved] = max(span[resolved], size)
        offsets = {}
        position = 0
        for name, length in span.items():
            offsets[name] = position
            position += length
        self.offset_by_chain_chrom = {
            chrom: offsets[resolved] for chrom, resolved in resolved_by_chain_chrom.items()
        }


def compute_chain_qc(
    chain_index: ChainIndex,
    ref_lengths: dict[str, int] | None = None,
    query_lengths: dict[str, int] | None = None,
    ref_aliases: dict[str, str] | None = None,
    query_aliases: dict[str, str] | None = None,
) -> ChainQc:
    """Chain sizes and scores, plus per-base depth histograms for both genomes.

    ref_lengths/query_lengths are each genome's sequence lengths (see
    _CoordinateSpace); *_aliases map chain chromosome names to those keys
    (e.g. a genome accessor's version-stripping aliases).
    """
    records = chain_index.records
    ref_space = _CoordinateSpace(_chain_sizes(records, "ref"), ref_lengths, ref_aliases)
    query_space = _CoordinateSpace(_chain_sizes(records, "query"), query_lengths, query_aliases)

    aligned_parts: list[np.ndarray] = []
    ref_keys: list[np.ndarray] = []
    query_keys: list[np.ndarray] = []
    with chain_index.indexed_chain_path.open("rb") as handle:
        for batch in _batches(records):
            sizes, dt, dq, counts = _parse_batch(handle, batch)
            first = np.cumsum(counts) - counts
            aligned_parts.append(np.add.reduceat(sizes, first))
            ref_keys.append(
                _event_keys(sizes, dt, counts, first, *_side_columns(batch, "ref", ref_space))
            )
            query_keys.append(
                _event_keys(sizes, dq, counts, first, *_side_columns(batch, "query", query_space))
            )

    aligned_bases = np.concatenate(aligned_parts) if aligned_parts else np.zeros(0, np.int64)
    scores = np.fromiter(map(attrgetter("score"), records), dtype=np.int64, count=len(records))
    ref_depth = _depth_summary(_concatenate_and_free(ref_keys), ref_space.total_bases)
    query_depth = _depth_summary(_concatenate_and_free(query_keys), query_space.total_bases)
    return ChainQc(aligned_bases, scores, ref_depth, query_depth)


def _chain_sizes(records: list[ChainRecord], side: str) -> dict[str, int]:
    # Each chain chromosome's size, from the chain headers of one side.
    chroms = map(attrgetter(f"{side}_chrom"), records)
    return dict(zip(chroms, map(attrgetter(f"{side}_size"), records)))


def _side_columns(batch: list[ChainRecord], side: str, space: _CoordinateSpace):
    # One side's per-chain span, strand and genome offset, as arrays. map()
    # with attrgetter keeps the per-chain work in C.
    def column(name: str) -> np.ndarray:
        values = map(attrgetter(f"{side}_{name}"), batch)
        return np.fromiter(values, dtype=np.int64, count=len(batch))

    minus = np.array(list(map(attrgetter(f"{side}_strand"), batch))) == "-"
    chroms = map(attrgetter(f"{side}_chrom"), batch)
    offsets = np.fromiter(
        map(space.offset_by_chain_chrom.__getitem__, chroms), dtype=np.int64, count=len(batch)
    )
    return column("start"), column("end"), minus, offsets


def _batches(records: list[ChainRecord]):
    # Consecutive runs of chains spanning about BATCH_BYTES of the file. The
    # index lists chains in file order, so each run is one contiguous read.
    start = 0
    while start < len(records):
        end = start
        batch_bytes = 0
        while end < len(records) and (end == start or batch_bytes < BATCH_BYTES):
            batch_bytes += records[end].byte_length
            end += 1
        yield records[start:end]
        start = end


def _parse_batch(handle, batch: list[ChainRecord]):
    # Every block's (size, dt, dq), plus each chain's block count. A chain's
    # last body line holds only its final block's size; it is padded with
    # zero gaps so every block parses as one 3-field row.
    base = batch[0].start_byte
    handle.seek(base)
    data = handle.read(batch[-1].start_byte + batch[-1].byte_length - base)
    bodies: list[bytes] = []
    counts = np.empty(len(batch), dtype=np.int64)
    for i, record in enumerate(batch):
        start = record.start_byte - base
        end = start + record.byte_length
        body = data[data.index(b"\n", start, end) + 1 : end].strip()
        if len(body[body.rfind(b"\n") + 1 :].split()) == 1:
            body += b" 0 0"
        bodies.append(body)
        counts[i] = body.count(b"\n") + 1
    blocks = pd.read_csv(
        io.BytesIO(b"\n".join(bodies)), sep=r"\s+", header=None, dtype=np.int64, engine="c"
    ).to_numpy()
    if blocks.shape != (int(counts.sum()), 3):
        raise ValueError(
            f"Could not parse the aligned blocks of chains {batch[0].chain_id} to "
            f"{batch[-1].chain_id}: every chain body line should be 'size dt dq', "
            "except each chain's last, which is 'size'"
        )
    return blocks[:, 0], blocks[:, 1], blocks[:, 2], counts


def _event_keys(sizes, gaps, counts, first, start, end, minus, offsets) -> np.ndarray:
    # One side's block intervals as sweep events in the genome's coordinate
    # space: key = 2 * position + 1 for a block start, 2 * position for a
    # block end, so a sort orders ends before starts at the same position.
    step = sizes + gaps
    before = np.cumsum(step) - step  # each block's offset from its chain's first block...
    before -= np.repeat(before[first], counts)  # ...counted from its own chain only
    # The index stores each chain's span in plus-strand coordinates; a
    # minus-strand chain's blocks run downward from the span's end.
    plus_start = np.where(
        np.repeat(minus, counts),
        np.repeat(end, counts) - before - sizes,
        np.repeat(start, counts) + before,
    )
    plus_end = plus_start + sizes
    offset = np.repeat(offsets, counts)
    keep = sizes > 0
    return np.concatenate(
        [2 * (offset + plus_start)[keep] + 1, 2 * (offset + plus_end)[keep]]
    )


def _concatenate_and_free(parts: list[np.ndarray]) -> np.ndarray:
    keys = np.concatenate(parts) if parts else np.zeros(0, np.int64)
    parts.clear()
    return keys


def _depth_summary(keys: np.ndarray, total_bases: int) -> DepthSummary:
    # Sweep the sorted events: depth is constant between consecutive events,
    # so each gap adds its length to the histogram at the depth reached.
    keys.sort()
    histogram = np.zeros(1, dtype=np.float64)
    depth_before = 0
    for start in range(0, keys.size, SWEEP_CHUNK):
        chunk = keys[start : start + SWEEP_CHUNK + 1]  # plus the next event, for the last gap
        n_events = min(SWEEP_CHUNK, keys.size - start)
        depth = depth_before + np.cumsum((chunk[:n_events] & 1) * 2 - 1)
        gaps = np.diff(chunk >> 1)
        counts = np.bincount(depth[: gaps.size], weights=gaps)
        if counts.size > histogram.size:
            counts[: histogram.size] += histogram
            histogram = counts
        else:
            histogram[: counts.size] += counts
        depth_before = int(depth[-1])
    histogram = np.rint(histogram).astype(np.int64)
    # Depth 0 is everything no block covers, including uncovered sequences.
    histogram[0] = total_bases - int(histogram[1:].sum())

    if total_bases <= 0:
        return DepthSummary(histogram, 0, 0.0, 0, 0.0)
    cumulative = np.cumsum(histogram)
    return DepthSummary(
        histogram=histogram,
        total_bases=total_bases,
        mean=float(np.dot(np.arange(histogram.size), histogram) / total_bases),
        median=int(np.searchsorted(cumulative, total_bases / 2)),
        covered_fraction=float(1 - histogram[0] / total_bases),
    )

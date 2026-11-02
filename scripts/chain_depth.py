#!/usr/bin/env python
"""Chain QC figures for a chain file on its own, without running stola.

    python scripts/chain_depth.py alignment.chain.gz --out-dir chain_qc \\
        [--ref-genome ref.fa --query-genome query.2bit]

Writes the same two figures a stola run puts in summary_figures/ (see
stola.chain_qc): chain size and score densities, and per-base aligned-block
depth across the query and reference genomes. Without the genome files,
depth covers only the sequences named in the chain file, using the sizes in
its chain headers; with them, it covers every base of each genome, as in a
stola run.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from stola.chain_qc import compute_chain_qc  # noqa: E402
from stola.genome_io import open_genome_accessor  # noqa: E402
from stola.index import build_index  # noqa: E402
from stola.plots import write_chain_qc_plots  # noqa: E402
from stola.workflow import CHAIN_DEPTH_FILENAME, CHAIN_SIZE_SCORE_FILENAME  # noqa: E402


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("chain_file", type=Path, help="Input .chain or .chain.gz file")
    parser.add_argument("--out-dir", required=True, type=Path, help="Directory for the two SVGs")
    parser.add_argument("--ref-genome", type=Path, help="Reference genome, FASTA or .2bit")
    parser.add_argument("--query-genome", type=Path, help="Query genome, FASTA or .2bit")
    args = parser.parse_args(argv)
    if (args.ref_genome is None) != (args.query_genome is None):
        parser.error("give both --ref-genome and --query-genome, or neither")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        chain_index = build_index(args.chain_file, Path(tmpdir) / "chain_index.json")
        genome_args = {}
        if args.ref_genome is not None:
            ref = open_genome_accessor(args.ref_genome, Path(tmpdir) / "ref_fasta.idx")
            query = open_genome_accessor(args.query_genome, Path(tmpdir) / "query_fasta.idx")
            genome_args = dict(
                ref_lengths=ref.sequence_lengths(),
                query_lengths=query.sequence_lengths(),
                ref_aliases=ref.aliases,
                query_aliases=query.aliases,
            )
            ref.close()
            query.close()
        qc = compute_chain_qc(chain_index, **genome_args)

    write_chain_qc_plots(
        qc, args.out_dir / CHAIN_SIZE_SCORE_FILENAME, args.out_dir / CHAIN_DEPTH_FILENAME
    )
    print(f"{qc.aligned_bases.size} chains, {int(qc.aligned_bases.sum())} aligned bases")
    for genome, depth in (("Query", qc.query_depth), ("Reference", qc.ref_depth)):
        print(
            f"{genome}: {depth.total_bases} bp, {100 * depth.covered_fraction:.1f}% covered, "
            f"mean depth {depth.mean:.2f}, median depth {depth.median}"
        )
    print(f"Wrote {CHAIN_SIZE_SCORE_FILENAME} and {CHAIN_DEPTH_FILENAME} -> {args.out_dir}")


if __name__ == "__main__":
    main()

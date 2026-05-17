from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def safe_link_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        dst.symlink_to(src.resolve())
    except OSError:
        shutil.copy2(src, dst)


def main() -> None:
    p = argparse.ArgumentParser(
        description="Prepare standardised Panaroo inputs from onebio Module B genome folders."
    )
    p.add_argument("--genomes", required=True, type=Path,
                   help="Path to work/genomes containing one subfolder per BioSample")
    p.add_argument("--out", required=True, type=Path,
                   help="Output folder for Panaroo-ready inputs")
    args = p.parse_args()

    genomes_dir = args.genomes.resolve()
    out_dir = args.out.resolve()
    input_dir = out_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    pairs_path = out_dir / "panaroo_pairs.tsv"
    missing_path = out_dir / "missing_inputs.tsv"

    pair_lines = []
    missing_lines = ["biosample_accession\treason"]

    n_ok = 0
    n_missing = 0

    for bs_dir in sorted(genomes_dir.iterdir()):
        if not bs_dir.is_dir():
            continue

        bs = bs_dir.name.strip()
        gff = bs_dir / "annotation.gff"
        fasta = bs_dir / "genome.fna"

        if not gff.exists() and not fasta.exists():
            missing_lines.append(f"{bs}\tmissing annotation.gff and genome.fna")
            n_missing += 1
            continue
        if not gff.exists():
            missing_lines.append(f"{bs}\tmissing annotation.gff")
            n_missing += 1
            continue
        if not fasta.exists():
            missing_lines.append(f"{bs}\tmissing genome.fna")
            n_missing += 1
            continue

        # Standardise filenames so Panaroo sample names become the BioSample IDs
        out_gff = input_dir / f"{bs}.gff"
        out_fa = input_dir / f"{bs}.fa"

        safe_link_or_copy(gff, out_gff)
        safe_link_or_copy(fasta, out_fa)

        pair_lines.append(f"{out_gff}\t{out_fa}")
        n_ok += 1

    pairs_path.write_text("\n".join(pair_lines) + ("\n" if pair_lines else ""), encoding="utf-8")
    missing_path.write_text("\n".join(missing_lines) + "\n", encoding="utf-8")

    print(f"[OK] Panaroo pair list: {pairs_path}")
    print(f"[OK] Standardised input dir: {input_dir}")
    print(f"[OK] Missing-input report: {missing_path}")
    print(f"[SUMMARY] usable={n_ok} missing={n_missing}")


if __name__ == "__main__":
    main()
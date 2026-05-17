from __future__ import annotations

import json
from pathlib import Path

from .io_qc import fastq_sample_stats
from .module_c import (
    _detect_run_technology,
    align_short_reads_bowtie2,
    align_long_reads_minimap2,
    collect_alignment_metrics,
)
from .utils import ensure_dir, log


def run_doctor(ref: Path, fastq: list[Path], out: Path, threads: int = 8, full: bool = False) -> Path:
    """Small mapping diagnostic for one reference and one or more FASTQs."""
    out = ensure_dir(Path(out))
    ref = Path(ref)
    fastqs = [Path(x) for x in fastq]
    qc = fastq_sample_stats(fastqs[0], max_reads=50000)
    tech = _detect_run_technology(read_qc=qc)
    align_dir = ensure_dir(out / "alignment")
    if tech == "ont":
        bam = align_long_reads_minimap2(ref, fastqs, align_dir, threads=threads)
        aligner = "minimap2"
    else:
        bam, _ = align_short_reads_bowtie2(ref, fastqs, align_dir, threads=threads)
        aligner = "bowtie2"
    mapping = collect_alignment_metrics(bam, align_dir)
    report = {"reference": str(ref), "fastq": [str(x) for x in fastqs], "technology_guess": tech, "aligner": aligner, "raw_qc": qc, "mapping": mapping}
    report_path = out / "doctor_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    log(f"Doctor report: {report_path}")
    return report_path

from pathlib import Path
import pandas as pd
from onebio_amr_atlas.manifest_bridge import checker_csv_to_manifest


def test_checker_to_manifest(tmp_path: Path):
    p = tmp_path / "checker.csv"
    p.write_text(
        "biosample_accession,assembly_accession_best,final_conservative_hit,rnaseq_available,rnaseq_run_ids,antibiotics_resistant\n"
        "SAMN000001,GCA_000000001.1,True,True,SRR000001,ampicillin\n"
    )
    manifest, ast, *_ = checker_csv_to_manifest(p, tmp_path / "out")
    m = pd.read_csv(manifest, sep="\t")
    a = pd.read_csv(ast, sep="\t")
    assert len(m) == 1
    assert m.loc[0, "run_list"] == "SRR000001"
    assert a.loc[0, "phenotype"] == "R"

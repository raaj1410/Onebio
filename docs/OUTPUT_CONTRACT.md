# Output contract

## Discovery output

The RNA-seq AMR checker emits species-specific CSVs and a combined RNA-seq-positive CSV.
The important columns are:

- `biosample_accession`
- `assembly_accession_best`
- `final_conservative_hit`
- `rnaseq_available`
- `rnaseq_run_ids`
- `ena_fastq_ftp`
- `ena_fastq_md5`
- `antibiotics_resistant`
- `antibiotics_susceptible`
- `antibiotics_intermediate`

## Manifest bridge output

`onebio-atlas prepare-manifest` creates:

- `manifest.tsv`
- `ast_long.tsv`
- `checker_selected_rows.csv`
- `manifest_qc.json`

## Module outputs

Default run folder layout:

```text
work/run1/
  genomes/
    <BioSample>/
      genome.fna
      annotation.gff
      locus_map.tsv
      amrfinder.tsv
      amr_hits.tsv
  rnaseq/
    runs/<RunID>/
      qc/
      trim/
      alignment/
      counts/
      run_metrics.json
    rnaseq_run_index.tsv
  features/
    X_genome.tsv
    X_expr_tpm.tsv
    X_fused.tsv
  atlas/
    atlas.sqlite
  kg/
    kg_index.tsv
    top10/*.dot
```

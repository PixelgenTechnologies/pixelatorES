# `edgelist_aggregate` QC fixture

QC JSON files from an nf-core/pixelator stub run whose read-processing front end
is `edgelist` + `aggregate` rather than `amplicon` + `demux` + `collapse`. Two
pools, four samples, hashed. Laid out like `qc_jsons`: one directory per stage
at the top level, with the layout PXL files at the root. This is what
`EXPERIMENT_SUMMARY` stages into `results/`, so it is what `data_folder` points
at; the `pixelator/<stage>/` nesting seen under a pipeline `--outdir` comes from
`publishDir` and never reaches the report.

Pool1 was sequenced twice, so `edgelist/` holds both `pool1.report.json` and
`pool1_reseq.report.json`. `edgelist` outputs are named after the fastq
basename, unlike every other stage, which is named after the pool or the
sample. The accompanying samplesheet keeps its duplicate rows, one pair per
fastq, since that is what produces the two reports. Each report names its pool
in `sample_id`, which is how the two are regrouped.

The layout PXL files are copies of the `qc_jsons_hashing` ones, renamed to this
fixture's samples, since that fixture is hashed across two pools in the same
shape. The run's `.parquet` outputs and the PXL files of the upstream stages are
omitted for size; tests that need those create empty placeholders in a
temporary directory. `pipeline_info/` and `logs/` are omitted because file
discovery ignores both.

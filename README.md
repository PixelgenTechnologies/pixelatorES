# Experiment Summary 

The **Experiment Summary** (ES) generates a self-contained `.html` report summarizing a multi-sample experiment, including QC metrics, statistics, and visualizations.

## Usage

The `pixelatorES.qmd` file can be moved, renamed, and run from different locations on your computer without modifications. Given a sample sheet specifying the samples in an experiment, it processes that data and generates an experiment summary report. The ES is always run from the directory that the `.qmd` file is in, and will by default generate the report in the same directory.

### Parameters

Parameters are defined in the YAML header of the Quarto file:

```YAML
params:
  sample_sheet: "../../../ES_test_data/test_samplesheet.csv"
  data_folder: "../../../ES_test_data/"
  control_markers: ["mIgG1", "mIgG2a", "mIgG2b"]
  do_harmonize: FALSE
  harmonization_vars: ["condition"]
  norm_method: "CLR"
  clustering_resolution: 1
  annotation_method: "Seurat"
  mc_cores: 1
  debug_mode: FALSE
  test_mode: FALSE
```

Users can modify parameters in the YAML header or pass them as arguments when rendering the file from the command line.

- **`sample_sheet`**: Path to the sample sheet `.csv` file. Must contain the columns `sample`, `sample_alias`, and `condition`.
- **`data_folder`**: Path to the input data folder.
- **`control_markers`**: Markers for control samples (e.g., `c('mIgG1', 'mIgG2a', 'mIgG2b')`).
- **`do_harmonize`**: Logical flag to harmonize data (`TRUE/FALSE`).
- **`harmonization_vars`**: Variables to harmonize (e.g., `'condition'`).
- **`norm_method`**: Normalization method (`'CLR'`).
- **`clustering_resolution`**: Resolution for clustering analysis.
- **`annotation_method`**: Method for cell type annotation (`'Seurat'` or `'nmf'`).
- **`mc_cores`**: Number of cores for parallel processing.
- **`debug_mode`**: Logical flag to enable debug mode (`TRUE/FALSE`).
- **`test_mode`**: Logical flag to enable test mode (`TRUE/FALSE`). 

#### Sample Sheet

The sample sheet describes your experiment and the samples in it, and sets the conditions for differential analyses.

The sample sheet should have the following columns, but may have additional columns that will be ignored:
| sample | sample_alias | condition |
|--------|--------------|-----------|
| Exp1_PBMC_unstimulated_sample01_S1 | S1 | Unstim_PBMC | 
| Exp1_PBMC_stimulated_sample02_S2 | S2 | Stim_PBMC |

`sample` is the name of the sample, and must match the sample file base name in the data folder., i.e. the name of the sample file without the .pxl extension.

`sample_alias` (optional) is the alias for the sample; a short, unique identifier for the sample that is used in the report. If not provided, `sample` will be used.

`condition` is the condition of the sample, which is used for a number of steps in the report, such as harmonization (by default), differential analyses, and visualization.

#### Data Folder

The data folder should contain the folder hierarchy of the `nf-core/pixelator` output. The ES will search for files within subfolders of the data folder to find `.pxl` files and `.json` files with data and statistics for each sample matching the `sample_name` in the sample sheet. Therefore, the only requirement of the data folder format is that it contains a subfolder with the following `nf-core/pixelator` output folders:

- `amplicon`
- `analysis`
- `collapse`
- `demux`
- `graph`
- `layout`
- `logs` (optional)
- `post_analysis` (optional)
- `[sample_name].layout.dataset.pxl`

The ES will read `.json` files to collect statistics and metadata for each sample, and will read `.pxl` files to collect the data for each sample. The ES will read one `.pxl` file per sample, prioritizing in the order `layout` > `post_analysis` > `analysis` > `graph`.

### Running the Report

To run the report, execute the following command in the terminal, passing parameters as arguments when rendering the file. For example:

```bash
quarto render pixelatorES.qmd -P sample_sheet="a_samplesheet.csv" -P data_folder="./data/"
```

To specify the output file name, use the `--output` flag. For example:

```bash
quarto render pixelatorES.qmd -P sample_sheet="a_samplesheet.csv" -P data_folder="./data/" --output name-of-report.html
```

The report can also be rendered in RStudio by opening a copy of the `.qmd` file and pressing "Render" in the top toolbar. This will generate a report in the same directory as the `.qmd` file with the same name (e.g., `pixelatorES.html`). 


#### Running the Experiment Summary in a Docker container

The ES can be run in a Docker container with `pixelatorR` pre-installed using the `run_ES_docker.sh` bash script in the `utility/` directory.
The script downloads the Docker image from ghcr (or runs one that is available locally) and runs the ES in the container.

##### Configuration

Before you can run the script, you need to follow the steps listed [here](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry).

Step 1: create a personal access token

```bash
export CR_PAT=YOUR_TOKEN
```

Step 2: login to the container registry

```bash
echo $CR_PAT | docker login ghcr.io -u USERNAME --password-stdin
```

##### Running the script

To run the ES in a Docker container, execute the following command in the terminal:

```
bash run_ES_docker.sh -d <data directory> -n <output name> -s <sample sheet>
```

If data files are located in multiple separate folders, you can specify multiple folders using `-d`:

```
bash run_ES_docker.sh -d <data directory 1> -d <data directory 2> -n <output name> -s <sample sheet>
```

If you want to run with a specific container image you can use the `-i` flag to specify the image name (e.g. `ghcr.io/pixelgentechnologies/pixelatores:pr-15`):

```
bash run_ES_docker.sh -d <data directory> -n <output name> -s <sample sheet> -i <image name>
```

You can also pass additional parameters to the ES using the `-P` flag. For example, to set the normalization method to `CLR`, you can run:

```bash
bash run_ES_docker.sh -d <data directory> -n <output name> -s <sample sheet> -P norm_method=CLR
```

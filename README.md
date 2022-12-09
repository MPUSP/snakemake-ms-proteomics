# snakemake-ms-proteomics

A pipeline for the automatic initial processing and quality control of mass spectrometry data.

- [snakemake-ms-proteomics](#snakemake-ms-proteomics)
  - [Pipeline overview](#pipeline-overview)
  - [Installation](#installation)
    - [Snakemake](#snakemake)
    - [Additional tools](#additional-tools)
    - [Fragpipe](#fragpipe)
  - [Running the pipeline](#running-the-pipeline)
    - [Input data](#input-data)
    - [Execution](#execution)
  - [Output](#output)
      - [](#)
  - [Authors](#authors)
  - [References](#references)

## Pipeline overview

The pipeline is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and processes MS data using the following steps:

1. Prepare `workflow` file (`python` script)
2. Fetch protein database from NCBI or use user-supplied fasta file (`python`, [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/))
3. Generate decoy proteins ([DecoyPyrat](https://github.com/wtsi-proteomics/DecoyPYrat))
4. Import raw files, search protein database ([fragpipe](https://fragpipe.nesvilab.org/))
5. Align feature maps using IonQuant ([fragpipe](https://fragpipe.nesvilab.org/))
6. Import quantified features, infer and quantify proteins ([R MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html))
7. Compare different biological conditions, export results ([R MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html))
8. Generate HTML report with embedded QC plots ([R markdown](https://rmarkdown.rstudio.com/))
9. Clean up temporary files after pipeline execution (`bash` script)

## Installation

### Snakemake

Step 1: Install snakemake with `micromamba` (or any another `conda` flavor). This step generates a new conda environment called `snakemake-ms-proteomics`, which will be used for all further installations.

```
micromamba create -c conda-forge -c bioconda -n snakemake-ms-proteomics snakemake
```

Step 2: Activate conda environment with snakemake

```
source /path/to/micromamba/bin/activate
conda activate snakemake-ms-proteomics
```

### Additional tools

Install [DecoyPyrat](https://github.com/wtsi-proteomics/DecoyPYrat) from `bioconda`.
This small tool can be used to generate superior decoy proteins from proteome `*.fasta` files.

```
micromamba install -c bioconda decoypyrat
```

Install NCBI datasets command line tool from `conda-forge`.

```
micromamba install -c conda-forge ncbi-datasets-cli
```

Install a fresh **R environment** with a set of custom packages instead of the default system-wide one.

```
micromamba install -c conda-forge r-essentials
```

Then open an R session and install packages from within R like this:

```
# CRAN packages
install.packages(c("tidyverse", "ggrepel", "scales", "dendextend", "ggpubr", "BiocManager"))

# Bioconductor packages
BiocManager::install("MSstats")
```

### Fragpipe

Step 1: Download the latest [fragpipe release from github](https://github.com/Nesvilab/FragPipe/releases)

- download to desired destination
- unzip archive
- on linux, test fragpipe by running in terminal:

```
cd fragpipe/bin
./fragpipe
```

This will start the GUI (graphical user interface). To finish the installation, continue with the following steps.

Step 2: Configure Fragpipe

- in the `Config` tab, download/update the required modules (MSFragger, IonQuant, Philosopher, EasyPQP)
- some of these require the user to approve the academic license agreement
- if automatic download from GUI fails, install e.g. `philospopher` manually
- download `philospopher` `*.deb` package from [github](https://github.com/Nesvilab/philosopher/releases/) and run:

```
sudo dpkg -i path/to/*.deb
```

Step 3: In order to work with **Thermo `*.raw`** files, Mono needs to be installed:

```
sudo apt install mono-devel
```

Step 4: Set python environment.

Finally, set the `python` environment in the `Config` tab to the conda environment created in the previous step, in order to fulfill all dependencies. If the GUI config tab shows complaints about missing packages, install python packages into the specified environment (e.g. `pandas`, `numpy`, `cython`):

```
micromamba install -c anaconda cython
```

Next, You can make adjustments to the pipeline (`workflow`) or to the sample sheet (`manifest`):

- the sample sheet (`manifest` in fragpipe language) defines the experimental setup, which includes e.g. paths to input files, type of experiment, group and replicate (see next section for details)
- an example for a manifest file for LFQ test data is included in the `workflows` directory
- the `workflow` defines the fragpipe modules and their parameters
- workflows can be customised for test purposes but, in the future, shall be chosen automatically by the pipeline

## Running the pipeline

### Input data

The pipeline requires the following input files:

1. mass spectrometry data, such as Thermo `*.raw` or `*.mzML` files
2. an (organism) database in `*.fasta` format _OR_ a NCBI refseq ID. Decoys (`rev_` prefix) will be added if necessary
3. a sample sheet in tab-separated format (aka `manifest` file)
4. a `workflow` file, the pipeline definition for fragpipe


The samplesheet file has the following structure with four mandatory columns and no header (example file: `test/input/config/samplesheet.tsv`). 

- `sample`: names/paths to raw files
- `condition`: experimental group, treatments
- `replicate`: replicate number, consecutively numbered. Repeating numbers (e.g. 1,2,1,2) will be treated as paired samples!
- `type`: the type of MS data, will be used to determine the workflow
- `control`: reference condition for testing differential abudandance

| sample   | condition   | replicate | type | control     |
| -------- | ----------- | --------- | ---- | ----------- |
| sample_1 | condition_1 | 1         | DDA  | condition_1 |
| sample_2 | condition_1 | 2         | DDA  | condition_1 |
| sample_3 | condition_2 | 3         | DDA  | condition_1 |
| sample_4 | condition_2 | 4         | DDA  | condition_1 |


For manual execution of the pipeline in fragpipe GUI, open the desired `manifest` and `workflow` files and set the correct paths to input files, database, and the python environment (see `Installation --> Fragpipe`).

### Execution

To run the pipeline from command line, change the working directory.

```
cd /path/to/snakemake-ms-proteomics
```

The top level folder `snakemake-ms-proteomics` contains the `snakefile`, which defines the snakemake rules that detail the different rules for execution of the pipeline.

By default, `snakemake` treats the _first_ rule as the target rule if not another rule is passed as command line argument. Therefore, the first rule defines the desired output files.

A config file, `config.yml`, is passed to snakemake which contains global or module-specific options. The global options are the paths to the test data files included in this repository (samplesheet, workflow, `*.fasta` database, and output folder). For test data they should work out of the box on linux systems. The only option that needs to be adjusted manually is the location of the `fragpipe` binary (`module options` --> `fragpipe` --> `path`).

Before running the entire pipeline, we can perform a dry run using:

```
snakemake --configfile 'config.yml' --dry-run
```

To run the complete pipeline with test files, execute the following command. The definition of the number of compute cores is mandatory.

```
snakemake --configfile 'config.yml' --cores 10
```

To supply options that override the defaults, run the pipeline like this:

```
snakemake --cores 10 \
  --config \
  samplesheet='test/input/config/samplesheet.fp-manifest' \
  database='test/input/database/database.fasta' \
  workflow='workflows/LFQ-MBR.workflow' \
  output='test/output/'
```

**Summary of mandatory arguments**

| argument    | type                   | details             | example                                                 |
| ----------- | ---------------------- | ------------------- | ------------------------------------------------------- |
| samplesheet | `*.tsv`                | tab-separated file  | `test/input/config/samplesheet.tsv`                     |
| database    | `*.fasta` OR refseq ID | plain text          | `test/input/database/database.fasta`, `GCF_000009045.1` |
| workflow    | `*.workflow`           | a fragpipe workflow | `workflows/LFQ-MBR.workflow`                            |
| output      | path                   | valid directory     | `test/output/`                                          |


## Output

The pipeline generates the following output from its modules:

#### 

<details markdown="1">
<summary>Output files</summary>

- `decoypyrat/`
  - `decoy_database.fasta`: Original `.fasta` file supplemented with randomized protein sequences.

</details>

## Authors

- The custom `snakemake`, `R`, `R markdown`, and `python` scripts were written by Michael Jahn, PhD
- Affiliation: Max-Planck-Unit for the Science of Pathogens ([MPUSP](https://www.mpusp.mpg.de/)), Berlin, Germany
- My [page on github](https://github.com/m-jahn)

## References

- Essential tools are linked in the top section of this document
- The core of this pipeline are the two external packages **fragpipe** and **MSstats**

**fragpipe references:**

1. Kong, A. T., Leprevost, F. V., Avtonomov, D. M., Mellacheruvu, D., & Nesvizhskii, A. I. (2017). _MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry–based proteomics_. Nature Methods, 14(5), 513-520.
2. da Veiga Leprevost, F., Haynes, S. E., Avtonomov, D. M., Chang, H. Y., Shanmugam, A. K., Mellacheruvu, D., Kong, A. T., & Nesvizhskii, A. I. (2020). _Philosopher: a versatile toolkit for shotgun proteomics data analysis_. Nature Methods, 17(9), 869-870.
3. Yu, F., Haynes, S. E., & Nesvizhskii, A. I. (2021). _IonQuant enables accurate and sensitive label-free quantification with FDR-controlled match-between-runs_. Molecular & Cellular Proteomics, 20.

**MSstats references:**

1. Choi M (2014). _MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments._ Bioinformatics, 30.


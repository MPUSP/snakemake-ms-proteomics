# snakemake-ms-proteomics

A pipeline for the automatic initial processing and quality control of mass spectrometry data.

- [snakemake-ms-proteomics](#snakemake-ms-proteomics)
  - [Description](#description)
  - [Installation](#installation)
    - [Snakemake](#snakemake)
    - [Additional tools](#additional-tools)
    - [Fragpipe](#fragpipe)
  - [Running the pipeline](#running-the-pipeline)
    - [Input data](#input-data)
    - [Execution](#execution)

## Description

To be added.

## Installation

### Snakemake

Step 1: Install snakemake with micromamba (if not yet installed). This step generates a new conda environment called `snakemake-ms-proteomics`.

```
micromamba create -c conda-forge -c bioconda -n snakemake-ms-proteomics snakemake
```

Step 2: Activate conda environment with snakemake

```
conda activate snakemake-ms-proteomics
```

### Additional tools

Install [DecoyPyrat](https://github.com/wtsi-proteomics/DecoyPYrat) from `bioconda`.
This small tool can be used to generate superior decoy proteins from automatically fetched NCBI genome/proteome `*.fasta` files.

```
micromamba install -c bioconda decoypyrat
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

Next, You can make adjustments to the pipeline ("workflow") or to the sample sheet:

- the sample sheet (`manifest` in fragpipe language) defines the experimental setup, which includes e.g. sample number, type, group and replicate
- the `workflow` defines the fragpipe modules and their parameters
- default manifest and workflow files for LFQ data are included in the `test/input/config` folder
- these files can be adjusted and downloaded for test purposes, but will be supplied automatically in the final pipeline

## Running the pipeline

### Input data

The pipeline requires the following _mandatory_ files:

1. mass spectrometry data, such as Thermo `*.raw` or `*.mzML` files
2. an (organism) database in fasta format. Decoys (`_rev` prefix) will be added if necessary
3. a sample sheet in tab-separated format (aka `manifest` file)
4. a `workflow` file, the pipeline definition for fragpipe


The samplesheet file has the following structure with four mandatory columns and no header (example file: `test/input/config/samplesheet.tsv`). The last column, here named `control`, defines the condition that will be used as reference for testing differential abudandance with MSstats.

| sample   | condition   | replicate | workflow | control     |
| -------- | ----------- | --------- | -------- | ----------- |
| sample_1 | condition_1 | 1         | DDA      | condition_1 |
| sample_2 | condition_1 | 2         | DDA      | condition_1 |
| sample_3 | condition_2 | 1         | DDA      | condition_1 |
| sample_4 | condition_2 | 2         | DDA      | condition_1 |


For manual execution of the pipeline in fragpipe GUI, open the desired `manifest` and `workflow` files and set the correct paths to input files, database, and the python environment (see `Installation --> Fragpipe`).

### Execution

To run the entire pipeline from command line, change the working directory.

```
cd /path/to/snakemake-ms-proteomics
```

The top level folder `snakemake-ms-proteomics` contains the `snakefile`, which defines the snakemake rules that detail the different rules for execution of the pipeline.

By default, `snakemake` treats the _first_ rule as the target rule if not another rule is passed as command line argument. Therefore, the first rule defines the desired output files.

A config file, `config.yml`, is passed to snakemake which contains global or module-specific options. The global options are the paths to the test data files included in this repository (samplesheet, workflow, (fasta) database, and output folder). They should work out of the box on linux systems. The only option that needs to be adjusted manually is the location of the `fragpipe` binary (`module options` --> `fragpipe` --> `path`).

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
  workflow='test/input/config/LFQ-MBR.workflow' \
  database='test/input/database/database.fasta' \
  output='test/output/'
```


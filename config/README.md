# snakemake-ms-proteomics

This workflow is a best-practice workflow for the automated analysis of mass spectrometry proteomics data. It currently supports automated analysis of data-dependent acquisition (DDA) data with label-free quantification. An extension by different wokflows (DIA, isotope labeling) is planned in the future. The workflow is mainly a wrapper for the excellent tools [fragpipe](https://fragpipe.nesvilab.org/) and [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html), with additional modules that supply and check the required input files, and generate reports. The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and processes MS data using the following steps:

1. Prepare `workflow` file (`python` script)
2. check user-supplied sample sheet (`python` script)
3. Fetch protein database from NCBI or use user-supplied fasta file (`python`, [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/))
4. Generate decoy proteins ([DecoyPyrat](https://github.com/wtsi-proteomics/DecoyPYrat))
5. Import raw files, search protein database ([fragpipe](https://fragpipe.nesvilab.org/))
6. Align feature maps using IonQuant ([fragpipe](https://fragpipe.nesvilab.org/))
7. Import quantified features, infer and quantify proteins ([R MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html))
8. Compare different biological conditions, export results ([R MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html))
9. Generate HTML report with embedded QC plots ([R markdown](https://rmarkdown.rstudio.com/))
10. Generate PDF report from HTML [weasyprint](https://weasyprint.org/)
11. Send out report by email (`python` script)
12. Clean up temporary files after workflow execution (`bash` script)

If you want to contribute, report issues, or suggest features, please get in touch on [github](https://github.com/MPUSP/snakemake-ms-proteomics).

## Installation

### Snakemake

Step 1: Install snakemake with `conda`, `mamba`, `micromamba` (or any another `conda` flavor). This step generates a new conda environment called `snakemake-ms-proteomics`, which will be used for all further installations.

```bash
conda create -c conda-forge -c bioconda -n snakemake-ms-proteomics snakemake
```

Step 2: Activate conda environment with snakemake

```bash
source /path/to/conda/bin/activate
conda activate snakemake-ms-proteomics
```

Alternatively, install `snakemake` using pip:

```bash
pip install snakemake
```

Or install `snakemake` globally from linux archives:

```bash
sudo apt install snakemake
```

### Fragpipe

Fragpipe is not available on `conda` or other package archives. However, to make the workflow as user-friendly as possible, the latest [fragpipe release from github](https://github.com/Nesvilab/FragPipe/releases) (currently v22.0) is automatically installed to the respective `conda` environment when using the workflow the first time. After installation, the GUI (graphical user interface) will pop up and ask to you to finish the installation by **downloading the missing modules MSFragger, IonQuant, and Philosopher**. This step is necessary to abide to license restrictions. From then on, fragpipe will run in `headless` mode through command line only.

All other dependencies for the workflow are **automatically pulled as `conda` environments** by snakemake.

## Running the workflow

### Input data

The workflow requires the following input files:

1. mass spectrometry data, such as Thermo `*.raw` or `*.mzML` files
2. an (organism) database in `*.fasta` format _OR_ a NCBI Refseq ID. Decoys (`rev_` prefix) will be added if necessary
3. a sample sheet in tab-separated format (aka `manifest` file)
4. a `workflow` file for fragpipe (see `resources` dir)

The samplesheet file has the following structure with four mandatory columns and no header (example file: `test/input/samplesheet/samplesheet.tsv`).

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

### Execution

To run the workflow from command line, change the working directory.

```bash
cd /path/to/snakemake-ms-proteomics
```

Adjust options in the default config file `config/config.yml`.
Before running the entire workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the complete workflow with test files using **`conda`**, execute the following command. The definition of the number of compute cores is mandatory.

```bash
snakemake --cores 10 --sdm conda --directory .test
```

To supply options that override the defaults, run the workflow like this:

```bash
snakemake --cores 10 --sdm conda --directory .test \
  --configfile 'config/config.yml' \
  --config \
  samplesheet='my/sample_sheet.tsv'
```

### Parameters

This table lists all **global parameters** to the workflow.

| parameter   | type                   | details             | example                                                 |
| ----------- | ---------------------- | ------------------- | ------------------------------------------------------- |
| samplesheet | `*.tsv`                | tab-separated file  | `test/input/config/samplesheet.tsv`                     |
| database    | `*.fasta` OR refseq ID | plain text          | `test/input/database/database.fasta`, `GCF_000009045.1` |
| workflow    | `*.workflow` OR string | a fragpipe workflow | `workflows/LFQ-MBR.workflow`, `from_samplesheet`        |

This table lists all **module-specific parameters** and their default values, as included in the `config.yml` file.

| module     | parameter        | default                            | details                                                          |
| ---------- | ---------------- | ---------------------------------- | ---------------------------------------------------------------- |
| decoypyrat | `cleavage_sites` | `KR`                               | amino acids residues used for decoy peptide generation           |
|            | `decoy_prefix`   | `rev`                              | decoy prefix appended to proteins names                          |
| fragpipe   | `target_dir`     | `share`                            | default path in conda env to store fragpipe                      |
|            | `executable`     | `fragpipe/bin/fragpipe`            | path to fragpipe executable                                      |
|            | `download`       | FragPipe-22.0 (see config)         | downlowd link to Fragpipe Github repo                            |
| msstats    | `logTrans`       | `2`                                | base for log fold change transformation                          |
|            | `normalization`  | `equalizeMedians`                  | normalization strategy for feature intensity, see MSstats manual |
|            | `featureSubset`  | `all`                              | which features to use for quantification                         |
|            | `summaryMethod`  | `TMP`                              | how to calculate protein from feature intensity                  |
|            | `MBimpute`       | `True`                             | Imputes missing values with Accelerated failure time model       |
| report     | `html`           | `True`                             | Generate HTLM report                                             |
|            | `pdf`            | `True`                             | Generate PDF report                                              |
| email      | `send`           | `False`                            | whether reports should send out by email                         |
|            | `port`           | `0`                                | default port for email server                                    |
|            | `smtp_server`    | `smtp.example.com`                 | smtp server address                                              |
|            | `smtp_user`      | `user`                             | smtp server user name                                            |
|            | `smtp_pw`        | `password`                         | smtp server user password                                        |
|            | `from`           | `sender@email.com`                 | sender's email address                                           |
|            | `to`             | `["receiver@email.com"]`           | receiver's email address(es), a list                             |
|            | `subject`        | `"Results MS proteomics workflow"` | subject line for email                                           |

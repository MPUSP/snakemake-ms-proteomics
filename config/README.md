## workflow overview

<!-- include logo-->
<img src="docs/images/logo.png" align="right" />

---

This workflow is a best-practice workflow for the automated analysis of mass spectrometry proteomics data. It currently supports automated analysis of data-dependent acquisition (DDA) data with label-free quantification. An extension by different workflows (DIA, isotope labeling) is planned in the future.

The workflow is mainly a wrapper for the excellent tools [fragpipe](https://fragpipe.nesvilab.org/) and [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html), with additional modules that supply and check the required input files, and generate reports. The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and processes MS data using the following steps:

1. Prepare `workflow` file (`python` script)
2. Check user-supplied sample sheet (`python` script)
3. Check user-supplied database FASTA file (`python` script)
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

## Running the workflow

### Input data

The workflow requires the following input files:

1. mass spectrometry data, such as Thermo `*.raw` or `*.mzML` files
2. an (organism) database in `*.fasta` format. Decoys (`rev_` prefix) will be added if necessary
3. a sample sheet in tab-separated format (aka `manifest` file)
4. a `workflow` file for fragpipe (see `resources` dir)

The samplesheet file has the following structure (example file: `test/input/samplesheet/samplesheet.tsv`).

- `sample`: names/paths to raw files
- `condition`: experimental group, treatments
- `replicate`: replicate number, consecutively numbered. Repeating numbers (e.g. 1, 2, 1, 2) will be treated as paired samples!
- `type`: the type of MS data, will be used to determine the workflow
- `control`: reference condition for testing differential abundance

| sample_name | raw_file                            | condition | replicate | MS_method | comparison |
| ----------- | ----------------------------------- | --------- | --------- | --------- | ---------- |
| WT_exp_1    | data/E1_220707_22C012_sample_01.raw | WT_exp    | 1         | DDA       | WT_exp     |
| WT_exp_2    | data/E1_220707_22C012_sample_02.raw | WT_exp    | 2         | DDA       | WT_exp     |
| WT_trans_1  | data/E1_220707_22C012_sample_03.raw | WT_trans  | 3         | DDA       | WT_exp     |
| WT_trans_2  | data/E1_220707_22C012_sample_04.raw | WT_trans  | 4         | DDA       | WT_exp     |

### Fragpipe

Fragpipe is not available on `conda` or other package archives. However, to make the workflow as user-friendly as possible, the latest [fragpipe release from github](https://github.com/Nesvilab/FragPipe/releases) (currently v22.0) is automatically installed to the respective `conda` environment when using the workflow the first time. After installation, the GUI (graphical user interface) will pop up and ask to you to finish the installation by **downloading the missing modules MSFragger, IonQuant, and Philosopher**. This step is necessary to abide to license restrictions. From then on, fragpipe will run in `headless` mode through command line only.

All other dependencies for the workflow are **automatically pulled as `conda` environments** by snakemake.

### Missing value imputation

Missing value imputation happens at different stages.
First, the default strategy for `fragpipe` is to use "match between runs", i.e. non-identified features in the MS1 spectra are cross-compared with other runs of the same experiment where MS2 identification is available. This strategy is based on actual quantification data and reduces the number of missing feature quantifications.

Second, `MSstats` imputes two kinds of missing values where absolutely no feature quantification is available: Missing values at random positions are removed during summarization, and missing values due to low abundance are imputed at the feature level via accelerated failure time model.

Missing value treatment can be controlled through the `MSstats` parameters `MBimpute` and others.
See the `MSstats` manual for more information.

## Output

The workflow generates the following output from its modules:

<details markdown="1">
<summary>samplesheet</summary>

- `samplesheet.tsv`: Samplesheet after checking file paths and options
- `log.txt`: Log file for this module

</details>

<details markdown="1">
<summary>workflow</summary>

- `workflow.txt`: Configuration file for `fragpipe`, determined from samplesheet.
- `log.txt`: Log file for this module

</details>

<details markdown="1">
<summary>database</summary>

- `database.fasta`: The downloaded or user-supplied `.fasta` file. In the latter case, the file is identical to the input.
- `log.txt`: Log file for this module

</details>

<details markdown="1">
<summary>decoypyrat</summary>

- `decoy_database.fasta`: Original `.fasta` file supplemented with randomized protein sequences.
- `log.txt`: Log file for this module

</details>

<details markdown="1">
<summary>fragpipe</summary>

- `[sample_name]/`: Directory containing sample specific output files for each run
- `combined_ion.tsv`: Quantification of ion intensity per peptide
- `combined_modified_peptide.tsv`: Quantification of peptide modifications
- `combined_peptide.tsv`: Quantification of peptides/features
- `combined_protein.tsv`: Quantification of proteins from petide, inferred by `fragpipe`
- `MSstats.csv`: Qunatification of petides/features, output from fragpipe served in `MSstats` friendly format
- other files such as logs, file lists, etc.
- `log.txt`: Log file for this module

</details>

<details markdown="1">
<summary>msstats</summary>

- `comparison_result.csv`: Main table with results about the comparison between different experimental conditions
- `feature_level_data.csv`: Feature-level quantification data processed by MSstats
- `model_qc.csv`: Table with data about the fitted quantification models from MSstats
- `protein_level_data.csv`: Protein-level quantification data processed by MSstats
- `uniprot.csv`: Optionally downloaded table with protein annotation from Uniprot
- `log.txt`: Log file for this module

</details>

<details markdown="1">
<summary>report</summary>

- `report.html`: Report with figures and tables
- `report.pdf`: Report with figures and tables in PDF format. Converted from HTML
- `log.txt`: Log file for this module

</details>

</details>

<details markdown="1">
<summary>email</summary>

- `log.txt`: Log file for this module

</details>

</details>

<details markdown="1">
<summary>clean_up</summary>

- `log.txt`: Log file for this module

</details>

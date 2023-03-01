#!/usr/bin/env Rscript

# load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(MSstats)
})


# Step 1: Import and process MS peptide quantifications
# -----------------------------------------------------------------------------
#
# parse MSstats specific parameters
list_msstats_config <- snakemake@params[["config_msstats"]]

# import sample sheet and add column names
df_sample_sheet <- read_tsv(
  snakemake@input[["samplesheet"]],
  show_col_types = FALSE,
  col_names = c("file_path", "condition", "replicate", "ms_type", "control")
)
write_lines(
  file = snakemake@log[["path"]],
  x = "MSSTATS: Imported sample sheet"
)


# import MSstats result table
df_msstats <- read_csv(
  snakemake@input[["table_msstats"]],
  show_col_types = FALSE,
  na = c("", "NA", "0")
)

write_lines(
  file = snakemake@log[["path"]], append = TRUE,
  x = "MSSTATS: Imported result table"
)

# convert input data frame to an MSstats experiment. This step applies:
# - log2 transformation
# - normalization (default: median)
# - quantification based on N features (default: all)
# - summary statistic: `TMP` = Tukey's median polish by default
# - impute missing values (default is `TRUE`)
result_msstats <- dataProcess(
  df_msstats,
  logTrans = list_msstats_config[["logTrans"]],
  normalization = list_msstats_config[["normalization"]],
  featureSubset = list_msstats_config[["featureSubset"]],
  summaryMethod = list_msstats_config[["summaryMethod"]],
  MBimpute = as.logical(list_msstats_config[["MBimpute"]]),
  use_log_file = TRUE,
  log_file_path = snakemake@log[["path"]],
  append = TRUE,
  verbose = FALSE
)

write_lines(
  file = snakemake@log[["path"]], append = TRUE,
  x = "MSSTATS: Processed feature input data"
)

# Step 2: Differential protein abundance
# -----------------------------------------------------------------------------
#
# Comparing conditions with `groupComparison`
df_comparison <- data.frame(
  condition = levels(result_msstats$ProteinLevelData$GROUP)
)

if ("control" %in% colnames(df_sample_sheet)) {
  df_comparison <- df_comparison %>%
    left_join(
      by = "condition",
      df_sample_sheet %>% select(condition, control) %>% distinct()
    )
} else {
  df_comparison <- df_comparison %>%
    mutate(control = condition[1])
}

df_comparison <- df_comparison %>%
  mutate(comparison = unname(Map(c, condition, control)))

mat_contrast <- MSstatsContrastMatrix(
  contrasts = df_comparison %>%
    filter(condition != control) %>%
    pull(comparison),
  conditions = df_comparison$condition
)

result_comparison <- groupComparison(
  contrast.matrix = mat_contrast,
  data = result_msstats,
  log_base = list_msstats_config$logTrans,
  use_log_file = TRUE,
  log_file_path = snakemake@log[["path"]],
  append = TRUE,
  verbose = FALSE
)

write_lines(
  file = snakemake@log[["path"]], append = TRUE,
  x = "MSSTATS: MSSTATS: Compared conditions with 'groupComparison'"
)

# Step 3: Retrieve annotation from uniprot
# -----------------------------------------------------------------------------
#
# function to retrieve organism's proteome from uniprot;
# organism ID is guessed from the first entry ID
get_uniprot_proteome <- function(id, backup = NULL) {
  server_error <- simpleError("")
  df_uniprot <- tryCatch(
    {
      uniprot_url_beacon <- paste0(
        "https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=",
        "accession%2Cid%2Corganism_name%2Cgene_oln%2Corganism_id&format=tsv",
        "&query=%28", id, "%29&size=500"
      )
      df_uniprot_beacon <- read_tsv(uniprot_url_beacon, col_types = cols())
      write_lines(
        file = snakemake@log[["path"]], append = TRUE,
        x = "MSSTATS: Retrieved protein annotation from Uniprot"
      )

      uniprot_org_id <- df_uniprot_beacon[[1, "Organism (ID)"]]
      write_lines(file = snakemake@log[["path"]], append = TRUE, paste0(
        x = "MSSTATS: guessing proteins belong to organism: ",
        df_uniprot_beacon[[1, "Organism"]], " based on first entry: ", id, "."
      ))
      uniprot_url <- paste0(
        "https://rest.uniprot.org/uniprotkb/stream?compressed=false&fields=",
        "accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_n",
        "ame%2Clength%2Cgene_oln%2Ccc_pathway%2Cgo_id%2Ccc_subcellular_loca",
        "tion%2Cxref_refseq&format=tsv&query=%28%28taxonomy_id%3A",
        uniprot_org_id, "%29%29"
      )
      df <- read_tsv(uniprot_url, col_types = cols()) %>%
        rename_with(.fn = function(x) gsub("[\\.\\;\\,\\ ]", "_", tolower(x))) %>%
        filter(
          !is.na(`gene_names_(ordered_locus)`),
          !duplicated(`gene_names_(ordered_locus)`)
        )
      write_lines(
        file = snakemake@log[["path"]], append = TRUE,
        x = paste0(
          "MSSTATS: downloaded proteome annotation for ",
          df_uniprot_beacon[[1, "Organism"]], " with ",
          length(unique(df$entry)), " unique entries."
        )
      )
      return(df)
    },
    error = function(server_error) {
      if (!is.null(backup)) {
        warning(
          paste0(
            "Uniprot server not available, ",
            "falling back to importing local Uniprot file."
          )
        )
        return(read_csv(backup, col_types = cols(), name_repair = tolower))
      } else {
        warning(
          paste0(
            "Uniprot server not available, ",
            "skipping protein annotation step."
          )
        )
        return(NULL)
      }
    }
  )
}


# retrieve proteome from uniprot
df_uniprot <- get_uniprot_proteome(
  id = as.character(result_msstats$ProteinLevelData$Protein[1])
)


# map refseq IDs to protein table
if (!is.null(df_uniprot)) {
  df_uniprot <- df_uniprot %>%
    separate(
      refseq,
      into = c("refseq", "refseq_2"),
      sep = ";", extra = "drop"
    )

  result_msstats$ProteinLevelData <- result_msstats$ProteinLevelData %>%
    mutate(Protein = as.character(Protein)) %>%
    left_join(df_uniprot, by = c("Protein" = "refseq"))

  result_comparison$ComparisonResult <- result_comparison$ComparisonResult %>%
    mutate(Protein = as.character(Protein)) %>%
    left_join(df_uniprot, by = c("Protein" = "refseq"))
} else {
  df_uniprot <- data.frame(entry = NULL, protein = NULL, refseq = NULL)
}


# Step 4: Export result tables
# -----------------------------------------------------------------------------
output_folder <- str_remove(snakemake@output[["protein_level_data"]], "protein_level_data.csv")
write_csv(result_msstats$FeatureLevelData, snakemake@output[["feature_level_data"]])
write_csv(result_msstats$ProteinLevelData, snakemake@output[["protein_level_data"]])
write_csv(result_comparison$ComparisonResult, snakemake@output[["comparison_result"]])
write_csv(result_comparison$ModelQC, snakemake@output[["model_qc"]])
write_csv(df_uniprot, snakemake@output[["uniprot"]])

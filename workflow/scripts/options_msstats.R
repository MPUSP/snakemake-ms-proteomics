options <- setClass("options", slots = list(input = "list", output = "list", params = "list", log = "list"))

snakemake <- new("options",
  input = list(
    samplesheet = "results/samplesheet/samplesheet.tsv",
    table_msstats = "results/fragpipe/msstats.csv"
  ),
  output = list(
    path = "results/msstats/test.csv"
  ),
  params = list(
    config_msstats = list(
      logTrans = 2,
      normalization = "equalizeMedians",
      featureSubset = "all",
      summaryMethod = "TMP",
      MBimpute = TRUE,
      use_log_file = FALSE
    )
  ),
  log = list(
    path = "results/msstats/log.txt"
  )
)

options <- setClass("options", slots=list(input="list", output="list", params="list"))

snakemake = new("options",
    input = list(
        samplesheet = "test/input/config/samplesheet.tsv",
        table_msstats = "test/output/fragpipe/MSstats.csv"
    ),
    output = list(
        path = "test/output/msstats/test.csv"
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
    )
)

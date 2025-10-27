extract_transcripts_from_genes <- function(oncogene_file, conversion_file, output_file) {
    oncogenes <- readLines(oncogene_file)
    conversion_df <- read.table(conversion_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
    oncogene_transcripts <- subset(conversion_df, gene_name %in% oncogenes)
    transcripts_list <- unique(oncogene_transcripts$transcript_id) #not necessary
    writeLines(transcripts_list, output_file)
}

extract_transcripts_from_genes(
  oncogene_file = "/home/Prostate_GAS/resources/oncogenes.genes",
  conversion_file = "/home/Prostate_GAS/resources/reference_data/genecode_V49_conversion_table.tsv",
  output_file = "/home/Prostate_GAS/resources/oncogenes_transcripts"
)

extract_transcripts_from_genes(
  oncogene_file = "/home/Prostate_GAS/resources/Oncosuppr.genes",
  conversion_file = "/home/Prostate_GAS/resources/reference_data/genecode_V49_conversion_table.tsv",
  output_file = "/home/Prostate_GAS/resources/oncosuppressor_transcripts"
)
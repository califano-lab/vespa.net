library(vespa)

# generate regulon
regulon<-hparacne2regulon(snakemake@input[["network"]],snakemake@input[["peptides"]],snakemake@params[["proteinlevel"]])
saveRDS(regulon, file=snakemake@output[["regulon"]])

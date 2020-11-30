library(phosphoviper)

# generate regulon
regulon<-hparacne2regulon(snakemake@input[["network"]],snakemake@input[["peptides"]])
saveRDS(regulon, file=snakemake@output[["regulon"]])

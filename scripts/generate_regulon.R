library(phosphoviper)

# generate regulon
regulon<-hparacne2regulon(snakemake@input[["network"]],snakemake@input[["peptides"]], snakemake@input[["matrix"]], confidence_threshold = 0, priors = FALSE)
saveRDS(regulon, file=snakemake@output[["regulon"]])

library(phosphoviper)

# generate regulon
regulon<-hparacne2regulon(snakemake@input[["network"]],snakemake@input[["peptides"]], snakemake@input[["matrix"]], likelihood_threshold = snakemake@params[["likelihood_threshold"]], priors = snakemake@params[["priors"]])
saveRDS(regulon, file=snakemake@output[["regulon"]])

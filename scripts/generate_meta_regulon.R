library(phosphoviper)
library(viper)

# import preprocessed data
qmx<-export2mx(readRDS(snakemake@input[["ref"]]))

# compute substrate-level VIPER matrix
if (length(snakemake@input[["substrate_regulons"]]) == 1) {
	substrate_regulons<-readRDS(snakemake@input[["substrate_regulons"]])
	vmx<-viper(qmx, substrate_regulons, minsize=5, pleiotropy = FALSE, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, targets = 10, penalty = 20, method = "adaptive"))
	vmxa<-export2mx(vmx2pv(vmx))
} else {
	vmxa<-qmx
}

# import single regulons
single_regulons<-snakemake@input[["regulons"]]
combined_regulons<-sapply(single_regulons, function(X){pruneRegulon(subset_regulon(readRDS(X), rownames(vmxa), min_size=5), 50)})

# combine and optimize regulons
meta_regulons<-optimizeRegulon(vmxa, combined_regulons, min_size=5)

print(meta_regulons)

# save meta regulons
saveRDS(meta_regulons, snakemake@output[["meta_regulons"]])

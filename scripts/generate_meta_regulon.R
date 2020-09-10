library(phosphoviper)
library(viper)

# import preprocessed data
qmx<-export2mx(readRDS(snakemake@input[["ref"]]))

# compute substrate-level VIPER matrix
if (length(snakemake@input[["substrate_regulons"]]) == 1) {
	substrate_regulons<-readRDS(snakemake@input[["substrate_regulons"]])
	qmx<-t(apply(qmx,1,function(X){X[is.na(X)]<-min(X,na.rm=TRUE);return(X)}))
	vmx<-viper(qmx, substrate_regulons, minsize=10, pleiotropy = TRUE, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, targets = 10, penalty = 10, method = "adaptive"), cores=snakemake@threads)

	vmxa<-export2mx(vmx2pv(vmx, fasta=snakemake@input[["fasta"]]))
} else {
	vmxa<-qmx
}

# import single regulons
single_regulons<-snakemake@input[["regulons"]]

if(length(single_regulons)>1) {
	# tune all regulons
	combined_regulons<-sapply(single_regulons, function(X){pruneRegulon(subset_regulon(readRDS(X), rownames(vmxa), min_size=10), 50)})

	# combine and optimize regulons
	meta_site_regulons<-optimizeRegulon(vmxa, combined_regulons, min_size=10)
} else {
	meta_site_regulons<-pruneRegulon(subset_regulon(readRDS(single_regulons), rownames(vmxa), min_size=10), 50)
}

# generate protein-level regulons
meta_protein_regulons<-optimizeRegulon(vmxa, regulator2protein(meta_site_regulons), min_size=10)

# save meta regulons
saveRDS(meta_site_regulons, snakemake@output[["meta_site_regulons"]])
saveRDS(meta_protein_regulons, snakemake@output[["meta_protein_regulons"]])

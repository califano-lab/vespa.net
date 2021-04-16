library(viper)
library(phosphoviper)

if (snakemake@params[["fill"]] == "NA") {
	fillvalues<-NA
} else {
	fillvalues<-snakemake@params[["fill"]]
}

# import preprocessed data
qmx<-export2mx(readRDS(snakemake@input[["ref"]]), fillvalues = fillvalues)

# compute substrate-level VIPER matrix
if (length(snakemake@input[["substrate_regulons"]]) == 1) {
	substrate_regulons<-readRDS(snakemake@input[["substrate_regulons"]])
	vmx<-viper(qmx, phosphoviper::pruneRegulon(phosphoviper::subsetRegulon(substrate_regulons, rownames(qmx), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]]), minsize=snakemake@params[["minimum_targets"]], pleiotropy = snakemake@params[["ct_correction"]], pleiotropyArgs = list(regulators = snakemake@params[["ct_regulators_threshold"]], shadow = snakemake@params[["ct_shadow_threshold"]], targets = snakemake@params[["ct_minimum_targets"]], penalty = snakemake@params[["ct_penalty"]], method = "adaptive"), cores=snakemake@threads)

	vmxa<-export2mx(vmx2pv(vmx, fasta=snakemake@input[["fasta"]]), fillvalues = fillvalues)
} else {
	vmxa<-qmx
}

# import single regulons
single_regulons<-snakemake@input[["regulons"]]

if(length(single_regulons)>1) {
	# tune all regulons
	combined_regulons<-sapply(single_regulons, function(X){phosphoviper::pruneRegulon(phosphoviper::subsetRegulon(readRDS(X), rownames(vmxa), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]])})

	# combine and optimize regulons
	meta_redundantsite_regulons<-optimizeRegulon(vmxa, combined_regulons, min_size=snakemake@params[["minimum_targets"]])
} else {
	meta_redundantsite_regulons<-phosphoviper::pruneRegulon(phosphoviper::subsetRegulon(readRDS(single_regulons), rownames(vmxa), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]])
}

# generate non-redundant, non-correlated site-level regulons
meta_site_regulons<-orthogonalRegulon(vmxa, meta_redundantsite_regulons, min_size=snakemake@params[["minimum_targets"]], min_size=snakemake@params[["orthogonal_cutoff"]])

# generate protein-level regulons
meta_protein_regulons<-optimizeRegulon(vmxa, regulator2protein(meta_redundantsite_regulons), min_size=snakemake@params[["minimum_targets"]])

# save meta regulons
saveRDS(meta_redundantsite_regulons, snakemake@output[["meta_redundantsite_regulons"]])
saveRDS(meta_site_regulons, snakemake@output[["meta_site_regulons"]])
saveRDS(meta_protein_regulons, snakemake@output[["meta_protein_regulons"]])

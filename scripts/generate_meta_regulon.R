library(viper)
library(phosphoviper)
library(stringr)

if (snakemake@params[["fill"]] == "NA") {
	fillvalues<-NA
} else {
	fillvalues<-snakemake@params[["fill"]]
}

# import preprocessed data
refmx<-readRDS(snakemake@input[["ref"]])
qmx<-export2mx(refmx, fillvalues = fillvalues)

# compute substrate-level VIPER matrix
if (length(snakemake@input[["substrate_regulons"]]) == 1) {
	substrate_regulons<-readRDS(snakemake@input[["substrate_regulons"]])
	vmx<-viper(qmx, phosphoviper::pruneRegulon(phosphoviper::subsetRegulon(substrate_regulons, rownames(qmx), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]]), minsize=snakemake@params[["minimum_targets"]], pleiotropy = snakemake@params[["ct_correction"]], pleiotropyArgs = list(regulators = snakemake@params[["ct_regulators_threshold"]], shadow = snakemake@params[["ct_shadow_threshold"]], targets = snakemake@params[["ct_minimum_targets"]], penalty = snakemake@params[["ct_penalty"]], method = "adaptive"), cores=snakemake@threads)

	if (str_detect(rownames(vmx)[1],":")) {
		vmxa<-export2mx(vmx2pv(vmx, fasta=NULL), fillvalues = fillvalues)
	} else {
		vmxa<-export2mx(vmx2pv(vmx, fasta=snakemake@input[["fasta"]]), fillvalues = fillvalues)
	}
} else {
	vmxa<-qmx
}

# compute VIPER signature if controls are present
if ("control" %in% colnames(refmx)) {
	control_runs<-unique(subset(refmx, control == TRUE)$run_id)
	vmxa_sig<-viperSignature(vmxa[,-control_runs], vmxa[,control_runs], per=1000)
} else {
	vmxa_sig<-vmxa
}

# import single regulons
single_regulons<-snakemake@input[["regulons"]]

if(length(single_regulons)>1) {
	# tune all regulons
	combined_regulons<-sapply(single_regulons, function(X){phosphoviper::pruneRegulon(phosphoviper::subsetRegulon(readRDS(X), rownames(vmxa), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]])})

	# combine and optimize regulons
	meta_redundantsite_regulons<-optimizeRegulon(vmxa_sig, combined_regulons)
} else {
	meta_redundantsite_regulons<-phosphoviper::pruneRegulon(phosphoviper::subsetRegulon(readRDS(single_regulons), rownames(vmxa), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]])
}

# generate non-redundant, non-correlated site-level regulons
meta_site_regulons<-orthogonalRegulon(vmxa_sig, meta_redundantsite_regulons, cutoff=snakemake@params[["orthogonal_cutoff"]])

# generate protein-level regulons
meta_protein_regulons<-optimizeRegulon(vmxa_sig, regulator2protein(meta_redundantsite_regulons))

# save meta regulons
saveRDS(meta_redundantsite_regulons, snakemake@output[["meta_redundantsite_regulons"]])
saveRDS(meta_site_regulons, snakemake@output[["meta_site_regulons"]])
saveRDS(meta_protein_regulons, snakemake@output[["meta_protein_regulons"]])

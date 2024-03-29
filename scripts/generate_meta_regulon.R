library(viper)
library(vespa)
library(stringr)

rn_mx<-function(qmx) {
  # rank based normalization of the matrix
  qmn.d1 <- t(t(apply(qmx, 2, rank, na.last = "keep"))/(colSums(!is.na(qmx)) + 1))
  
  # rank based estimation of single sample gene expression signature across the matrix
  qmn.norm <- t(apply(qmn.d1, 1, rank, na.last = "keep"))/(rowSums(!is.na(qmn.d1)) + 1)

  return(qmn.norm)
}

zscore_mx<-function(qmx) {
	qmn<-t(apply(t(qmx),2,function(X){return((X-mean(X, na.rm=TRUE))/sd(X, na.rm=TRUE))}))

	return(qmn)
}

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
	vmx<-viper(qmx, vespa::pruneRegulon(vespa::subsetRegulon(substrate_regulons, rownames(qmx), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]]), minsize=snakemake@params[["minimum_targets"]], pleiotropy = snakemake@params[["ct_correction"]], pleiotropyArgs = list(regulators = snakemake@params[["ct_regulators_threshold"]], shadow = snakemake@params[["ct_shadow_threshold"]], targets = snakemake@params[["ct_minimum_targets"]], penalty = snakemake@params[["ct_penalty"]], method = "adaptive"), cores=snakemake@threads)

	if (str_detect(rownames(vmx)[1],":")) {
		vmxa<-export2mx(vmx2pv(vmx, fasta=NULL), fillvalues = fillvalues)
	} else {
		vmxa<-export2mx(vmx2pv(vmx, fasta=snakemake@input[["fasta"]]), fillvalues = fillvalues)
	}
} else {
	vmxa<-qmx
}

# transform matrix
if (snakemake@params[["transform"]]=="rank") {
	vmxa<-rn_mx(vmxa)
} else if (snakemake@params[["transform"]]=="zscore") {
	vmxa<-zscore_mx(vmxa)
}

# compute VIPER signature if controls are present
if ("control" %in% colnames(refmx)) {
	target_runs<-unique(subset(refmx, control == FALSE)$run_id)
	control_runs<-unique(subset(refmx, control == TRUE)$run_id)

	vmxa_sig<-viperSignature(vmxa[,target_runs], vmxa[,control_runs], per=1000)
} else {
	vmxa_sig<-vmxa
}

# import single siteregulons
single_siteregulons<-snakemake@input[["siteregulons"]]

if(length(single_siteregulons)>1) {
	# tune all regulons
	combined_siteregulons<-sapply(single_siteregulons, function(X){vespa::pruneRegulon(vespa::subsetRegulon(readRDS(X), rownames(vmxa), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]])})

	# combine and optimize regulons
	meta_redundantsite_regulons<-optimizeRegulon(vmxa_sig, combined_siteregulons)
} else {
	meta_redundantsite_regulons<-vespa::pruneRegulon(vespa::subsetRegulon(readRDS(single_siteregulons), rownames(vmxa), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]])
}

# generate non-redundant, non-correlated site-level regulons
meta_site_regulons<-orthogonalRegulon(vmxa_sig, meta_redundantsite_regulons, cutoff=snakemake@params[["orthogonal_cutoff"]])

# import single proteinregulons
single_proteinregulons<-snakemake@input[["proteinregulons"]]

if(length(single_proteinregulons)>1) {
	# tune all regulons
	combined_proteinregulons<-sapply(single_proteinregulons, function(X){vespa::pruneRegulon(vespa::subsetRegulon(readRDS(X), rownames(vmxa), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]])})

	# combine and optimize regulons
	meta_protein_regulons<-optimizeRegulon(vmxa_sig, combined_proteinregulons)
} else {
	meta_protein_regulons<-vespa::pruneRegulon(vespa::subsetRegulon(readRDS(single_proteinregulons), rownames(vmxa), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]])
}

# save meta regulons
saveRDS(meta_redundantsite_regulons, snakemake@output[["meta_redundantsite_regulons"]])
saveRDS(meta_site_regulons, snakemake@output[["meta_site_regulons"]])
saveRDS(meta_protein_regulons, snakemake@output[["meta_protein_regulons"]])

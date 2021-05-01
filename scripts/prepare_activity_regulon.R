library(viper)
library(phosphoviper)
library(phosphoviper.db)
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
phospho<-readRDS(snakemake@input[["phospho"]])
proteo<-readRDS(snakemake@input[["proteo"]])
qmx<-export2mx(phospho, fillvalues = fillvalues)

# transform matrix
if (snakemake@params[["transform"]]=="rank") {
	qmx<-rn_mx(qmx)
} else if (snakemake@params[["transform"]]=="zscore") {
	qmx<-zscore_mx(qmx)
}

# import substrate regulon
meta_substrate_regulons<-readRDS(snakemake@input[["meta_substrate_regulons"]])
print(meta_substrate_regulons)

# run VIPER
vmx<-viper(qmx, phosphoviper::pruneRegulon(phosphoviper::subsetRegulon(meta_substrate_regulons, rownames(qmx), snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]]), minsize=snakemake@params[["minimum_targets"]], pleiotropy = snakemake@params[["ct_correction"]], pleiotropyArgs = list(regulators = snakemake@params[["ct_regulators_threshold"]], shadow = snakemake@params[["ct_shadow_threshold"]], targets = snakemake@params[["ct_minimum_targets"]], penalty = snakemake@params[["ct_penalty"]], method = "adaptive"), cores=snakemake@threads)

# convert to phosphoVIPER list
if (str_detect(rownames(vmx)[1],":")) {
        pvl<-vmx2pv(vmx, fasta=NULL)
} else {
        pvl<-vmx2pv(vmx, fasta=snakemake@input[["fasta"]])
}

# check if proteo-level abundances are present
if (identical(phospho, proteo)) {
	yqml<-pvl
} else {
	yqml<-rbind(proteo, pvl, fill = TRUE)
}

# export hpARACNe input files
export2hparacne(yqml, dirname(snakemake@output[["kinases"]]), kinases, phosphatases, list("hsm"=hsmpf,"pc"=pcdb,"lp"=lpadb), confidence_threshold = snakemake@params[["hsm_threshold"]], restrict_interactions=TRUE, interaction_level = "activity")

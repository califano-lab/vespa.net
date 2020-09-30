library(phosphoviper)
library(viper)

# import preprocessed data
proteo<-readRDS(snakemake@input[["proteo"]])
qmx<-export2mx(readRDS(snakemake@input[["phospho"]]))
qmx<-t(apply(qmx,1,function(X){X[is.na(X)]<-min(X,na.rm=TRUE);return(X)}))

# import dDPI substrate regulon
meta_substrate_regulons<-readRDS(snakemake@input[["meta_substrate_regulons"]])
print(meta_substrate_regulons)

# run VIPER
vmx<-viper(qmx, meta_substrate_regulons, minsize=snakemake@params[["minimum_targets"]], pleiotropy = snakemake@params[["ct_correction"]], pleiotropyArgs = list(regulators = snakemake@params[["ct_regulators_threshold"]], shadow = snakemake@params[["ct_shadow_threshold"]], targets = snakemake@params[["ct_minimum_targets"]], penalty = snakemake@params[["ct_penalty"]], method = "adaptive"), cores=snakemake@threads)

# convert to phosphoVIPER list
pvl<-vmx2pv(vmx, fasta=snakemake@input[["fasta"]])

# check if proteo-level abundances are present
if (snakemake@input[["phospho"]] == snakemake@input[["proteo"]]) {
	yqml<-pvl
} else {
	yqml<-rbind(proteo, pvl, fill = TRUE)
}

# export hpARACNe input files
export2hparacne(yqml, dirname(snakemake@output[["kinases"]]), kinases, phosphatases, hsmpf, confidence_threshold = snakemake@params[["hsm_threshold"]], interaction_level = "activity")

library(viper)
library(vespa)
library(vespa.db)
library(stringr)

if (snakemake@params[["fill"]] == "NA") {
    fillvalues<-NA
} else {
    fillvalues<-snakemake@params[["fill"]]
}

# import preprocessed data
phospho<-readRDS(snakemake@input[["phospho"]])
proteo<-readRDS(snakemake@input[["proteo"]])
qmx<-export2mx(phospho, fillvalues = fillvalues)

# import substrate regulon
meta_substrate_regulons<-readRDS(snakemake@input[["meta_substrate_regulons"]])
print(meta_substrate_regulons)

# run VIPER
vmx<-viper(qmx, vespa::pruneRegulon(vespa::subsetRegulon(meta_substrate_regulons, rownames(qmx), snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]], adaptive=snakemake@params[["adaptive"]]), minsize=snakemake@params[["minimum_targets"]], pleiotropy = snakemake@params[["ct_correction"]], pleiotropyArgs = list(regulators = snakemake@params[["ct_regulators_threshold"]], shadow = snakemake@params[["ct_shadow_threshold"]], targets = snakemake@params[["ct_minimum_targets"]], penalty = snakemake@params[["ct_penalty"]], method = "adaptive"), cores=snakemake@threads)

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

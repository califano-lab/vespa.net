library(phosphoviper)
library(viper)

# import preprocessed data
phospho<-readRDS(snakemake@input[["phospho"]])
proteo<-readRDS(snakemake@input[["proteo"]])

# import dDPI substrate regulon
meta_substrate_regulons<-readRDS(snakemake@input[["meta_substrate_regulons"]])
print(meta_substrate_regulons)

# run VIPER
vmx<-viper(export2mx(phospho), meta_substrate_regulons, minsize=5, pleiotropy = FALSE, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, targets = 10, penalty = 20, method = "adaptive"))

# convert to phosphoVIPER list
pvl<-vmx2pv(vmx)

# check if proteo-level abundances are present
if (snakemake@input[["phospho"]] == snakemake@input[["proteo"]]) {
	yqml<-pvl
} else {
	yqml<-rbind(proteo, pvl, fill = TRUE)
}

# export hpARACNe input files
export2hparacne(yqml, dirname(snakemake@output[["kinases"]]), kinases, phosphatases, hsmp, confidence_threshold = 0.05, interaction_level = "activity")

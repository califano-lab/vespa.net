library(phosphoviper)
library(viper)

# import preprocessed data
phospho<-readRDS(snakemake@input[["phospho"]])
proteo<-readRDS(snakemake@input[["proteo"]])

# import dDPI substrate regulon
ddpihsm_substrate_regulon<-readRDS(snakemake@input[["ddpihsm_substrate_regulon"]])
print(ddpihsm_substrate_regulon)

# run VIPER
vmx<-viper(export2mx(phospho), ddpihsm_substrate_regulon, minsize=5, pleiotropy = FALSE, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, targets = 10, penalty = 20, method = "adaptive"))

# convert to phosphoVIPER list
pvl<-vmx2pv(vmx)

# check if proteo-level abundances are present
if (snakemake@input[["phospho"]] == snakemake@input[["proteo"]]) {
	yqml<-pvl
} else {
	yqml<-rbind(proteo, pvl, fill = TRUE)
}

# export hpARACNe input files
export2hparacne(yqml, dirname(snakemake@output[["kinases"]]), kinases, phosphatases, hsmpf)

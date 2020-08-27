library(phosphoviper)
library(viper)

# import preprocessed data
phospho<-readRDS(snakemake@input[["phospho"]])
proteo<-readRDS(snakemake@input[["proteo"]])

# import dDPI substrate regulon
ddpi_substrate_regulon<-readRDS(snakemake@input[["ddpi_substrate_regulon"]])
print(ddpi_substrate_regulon)
ddpi_substrate_regulon<-pruneRegulon(subset_regulon(ddpi_substrate_regulon, rownames(export2mx(phospho)), min_size=5), 50)

# import HSM/D substrate regulon
hsm_substrate_regulon<-readRDS(snakemake@input[["hsm_substrate_regulon"]])
print(hsm_substrate_regulon)
hsm_substrate_regulon<-pruneRegulon(subset_regulon(hsm_substrate_regulon, rownames(export2mx(phospho)), min_size=5), 50)

# combine and optimize regulon
combined_substrate_regulon<-list("ddpi_substrate_regulon"=ddpi_substrate_regulon, "hsm_substrate_regulon"=hsm_substrate_regulon)
combined_substrate_regulon_optimized<-optimizeRegulon(export2mx(phospho), combined_substrate_regulon, min_size=5)

# run VIPER
vmx<-viper(export2mx(phospho), combined_substrate_regulon_optimized, minsize=5, pleiotropy = FALSE, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, targets = 10, penalty = 20, method = "adaptive"))

# convert to phosphoVIPER list
pvl<-vmx2pv(vmx)

# export hpARACNe input files
yqml<-rbind(proteo, pvl, fill = TRUE)
export2hparacne(yqml, dirname(snakemake@output[["kinases"]]), kinases, phosphatases, hsmpf, confidence_threshold = 0)

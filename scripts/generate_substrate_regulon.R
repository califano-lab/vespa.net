library(phosphoviper)
library(viper)

# import preprocessed data
qmx<-export2mx(readRDS(snakemake@input[["ref"]]))

# import dDPI substrate regulon
ddpi_substrate_regulon<-readRDS(snakemake@input[["ddpi_substrate_regulon"]])
print(ddpi_substrate_regulon)
ddpi_substrate_regulon<-pruneRegulon(subset_regulon(ddpi_substrate_regulon, rownames(qmx), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]])

# import HSM/D substrate regulon
hsm_substrate_regulon<-readRDS(snakemake@input[["hsm_substrate_regulon"]])
print(hsm_substrate_regulon)
hsm_substrate_regulon<-pruneRegulon(subset_regulon(hsm_substrate_regulon, rownames(qmx), min_size=snakemake@params[["minimum_targets"]]), snakemake@params[["maximum_targets"]])

# combine and optimize regulon
ddpihsm_substrate_regulon<-list("ddpi_substrate_regulon"=ddpi_substrate_regulon, "hsm_substrate_regulon"=hsm_substrate_regulon)
ddpihsm_substrate_regulon_optimized<-optimizeRegulon(qmx, ddpihsm_substrate_regulon, min_size=snakemake@params[["minimum_targets"]])

# save meta regulon
saveRDS(ddpihsm_substrate_regulon_optimized, snakemake@output[["ddpihsm_substrate_regulon"]])

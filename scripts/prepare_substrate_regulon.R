library(vespa)
library(vespa.db)

# import preprocessed data
phospho<-readRDS(snakemake@input[["phospho"]])
proteo<-readRDS(snakemake@input[["proteo"]])

# check if proteo-level abundances are present
if (identical(phospho, proteo)) {
	qml<-phospho
} else {
	qml<-rbind(proteo, phospho, fill=TRUE)
}

# check if target sites need to be restricted
if (snakemake@params[["restrict_peptides"]]) {
	target_sites<-unique(readRDS(snakemake@input[["ref"]])$site_id)
} else {
	target_sites<-NULL
}

# export hpARACNe input files
export2hparacne(qml, dirname(snakemake@output[["kinases"]]), kinases, phosphatases, list("hsm"=hsmpf,"pc"=pcdb,"lp"=lpsdb), target_sites, confidence_threshold = snakemake@params[["hsm_threshold"]], restrict_interactions=TRUE, interaction_level = "substrate")

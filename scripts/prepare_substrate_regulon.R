library(phosphoviper)

# import preprocessed data
phospho<-readRDS(snakemake@input[["phospho"]])
proteo<-readRDS(snakemake@input[["proteo"]])

# filter phosphopeptides
ref<-readRDS(snakemake@input[["ref"]])
phospho<-subset(phospho, site_id %in% unique(ref$site_id))

qml<-rbind(proteo, phospho, fill=TRUE)

# export hpARACNe input files
export2hparacne(qml, dirname(snakemake@output[["kinases"]]), kinases, phosphatases, hsmd, confidence_threshold = 0)

library(data.table)

peptides<-as.data.frame(fread(snakemake@input[["peptides"]]))
sitenetwork<-as.data.frame(fread(list.files(dirname(snakemake@input[["iteration"]]),"bootstrapNetwork", full.names=TRUE)))

proteinnetwork<-sitenetwork
proteinnetwork$Regulator<-peptides$protein_id[match(proteinnetwork$Regulator, peptides$peptide_id)]

write.table(sitenetwork, file=snakemake@output[["sitenetwork"]], quote=FALSE, row.names=FALSE, sep="\t")
write.table(proteinnetwork, file=snakemake@output[["proteinnetwork"]], quote=FALSE, row.names=FALSE, sep="\t")
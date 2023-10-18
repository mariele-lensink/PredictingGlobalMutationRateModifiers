library(data.table)

gff<-fread("data/TAIR10_GFF3_genes.gff")
gff<-gff[V3=="gene" & V1 %in% c("Chr1","Chr2","Chr3","Chr4","Chr5"),]
gff<-gff[,c(1,4,5)]

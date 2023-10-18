library(data.table)

gene_index<-fread("gene100.csv")
gene_index_subset<-gene_index[1:5,1:2]

#I read in a tair file, its filtered for lines with "gene" to distinguish coding vs noncoding (room for improvement here)
#It creates a logical vector, 'isgene' that will determine the probability of that location being determined as "under selection"

  #gene_vector<-vector(mode="logical",length = max(gene_inds))
lynchdt<-data.table(bp = seq(1,max(gene_index_subset)))
lynchdt[, isgene := bp %in% unlist(Map(`:`, gene_index_subset$V1, gene_index_subset$V2))]
#Inference of the Distribution of Selection Coefficients for New Nonsynonymous Mutations Using Large Samples
#juat picked a random shape and rate from the above paper for the sake of getting it working
lynchdt[, dfeMean := ifelse(selection, rgamma(.N, shape = 0.215, rate = 567.1),0)]

                        



lynchdt<-data.table("isgene"= logical(),"isselected"=logical(),"s"=numeric())
lynchdt<-data.table()
lynchdt$isgene<-populate_genecol(gene_index_subset)

dt$isselected<-populate_Scol(isgene)


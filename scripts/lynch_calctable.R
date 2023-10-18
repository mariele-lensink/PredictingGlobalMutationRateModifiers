library(data.table)
library(dplyr)

#read in the necessary input files
gene_index<-fread("GFF_wCHROM.txt")
gene_index_subset<-gene_index[1:50,]
#make a column indicating every basepair (to expand from gene ranges in input file)
lynchdt<-data.table(bp = seq(1,max(gene_index_subset)))
#read in normalized chipseq table
chipseq<-fread("data/normalized_depth_chipseq.csv")

#####Make pere BP Genome Table#####
#make a column T/F of basepair falls within gene range from input file
lynchdt[, isgene := bp %in% unlist(Map(`:`, gene_index_subset$start, gene_index_subset$stop))]
# Create a new column "selection" based on "isgene" column
lynchdt[, selection := ifelse(isgene, runif(.N) <= 0.7, runif(.N) <= 0.3)]
#Inference of the Distribution of Selection Coefficients for New Nonsynonymous Mutations Using Large Samples
#juat picked a random shape and rate from the above paper for the sake of getting it working
lynchdt[, dfeMean := ifelse(selection, rgamma(.N, shape = 0.215, rate = 567.1),0)]
#adding enrichment by basepair from chip seq
lynchdt[,enrichment:=chipseq$normalized_depth]
#adding l column, based on enrighment, should be modified
lynchdt[,l:= ifelse(enrichment >0,TRUE, FALSE)]

#####Variables for equation#####
#write out on paper, math theorem 
Ne<-300000
#why you chose each variable, from a paper
u<-0.000000001
du<-10 #consider as relative repair
s<-as.numeric(lynchdt[selection == TRUE,mean(dfeMean)])##make this include and L
lsegment<-as.numeric(lynchdt[l==T,.N])#check,make sure if contingent on affecting l
#proportion of l segment that is under selection
p<-as.numeric(lynchdt[selection == T & l ==T,.N]/lsegment)

ranknum<-Ne*u*du*s*lsegment*p


#ranknum critical threshold greater than 1
#find parameter space where drift overpowers selection

#compare distribution of raw dfe values in s column to what would be expected under "randomness"
#randomize L vector and shuffle it
#show density distrubtion of this effect (global dfe)
#when theres targeted repair toward d regions, we predict that global mutation rate will increase
#mess with u, can it drift upwards 
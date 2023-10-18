chipfiles <- c("chipseqdata/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
               "chipseqdata/GSM4668651_Col0_rep2_H3K4me1_ChIP_unique_reads.bedGraph.gz"
)
infiles<-c("data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz", 
           "data/GSM4668652_Col0_rep2_H3K4me1_Input_unique_reads.bedGraph.gz")

bg<-fread("chipseqdata/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz", sep = "\t", header = FALSE)

colnames(bg)<-c("CHROM","START","STOP","DEPTH")
bg$window<-rep(1:1000, length.out = nrow(bg))

w=1
for(w in 1:1000){
  t<-Sys.time()
  out<-rbindlist(apply(bg[window==w], 1, function(row){
    chr=row["CHROM"]
    start=as.numeric(row["START"])+1
    stop=as.numeric(row["STOP"])
    depth=as.numeric(row["DEPTH"])
    long<-data.table(CHROM=chr, POS=start:stop, DEPTH=depth)
   
    data.table(CHROM=chr, POS=start:stop, DEPTH=depth)
    
  }))
  fwrite(out, paste0("~/Desktop/",w,"_long_bg_rep1_H3K4me1.txt"))
}


library(seqinr)

genome<-read.fasta("~/Downloads/TAIR10_chr_all.fas.gz")

long_table<-rbindlist(lapply(names(genome), function(x){
  CHROM=gsub("Chr","", x)
  Length=length(genome[[x]])
  data.table(CHROM, POS=1:Length)
  
}))



final<-merge(long_table, long_bg, by=c("CHROM","POS"), all=T)
final$DEPTH[is.na(final$DEPTH)]<-0


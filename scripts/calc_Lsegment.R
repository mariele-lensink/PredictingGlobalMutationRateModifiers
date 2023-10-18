library(data.table)
library(parallel)
library(R.utils)
library(ggplot2)
library(doParallel)
library(foreach)

calculate_average_depth_perbp <- function(replicate_files) {

  #read in each replicate file
  rep_values<-lapply(replicate_files, function(bedfile) {
    cl <- makeCluster(detectCores()-2)
    bg<-fread(bedfile, sep = "\t", header = FALSE)
    colnames(bg)<-c("CHROM","START","STOP","DEPTH")
    #get rid of mt and pt DNA and chromosome
    bg<-bg[CHROM %in% c(1:5),]
    pb <- txtProgressBar(min = 0, max = nrow(bg), style = 3)
    
    # A helper function to update the progress bar
    .progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    registerDoParallel(cl)
    #parallelized row expansion
    list_of_tables <- foreach(row=split(bg,1:nrow(bg)), .packages = 'data.table') %dopar% {
      chr <- row$CHROM
      start <- as.numeric(row$START) + 1
      stop <- as.numeric(row$STOP)
      depth <- as.numeric(row$DEPTH)
      data.table(CHROM=chr, POS=start:stop, DEPTH=depth)
    }
    close(pb)
    dt<-rbindlist(list_of_tables)
    stopCluster(cl)
    return(dt)
  })
  # Combine all data.tables into a single data.table
  combined_dt <- rbindlist(rep_values)
  # Calculate the average depth for each "bp"
  average_depth_dt <- combined_dt[, .(avg_depth = mean(DEPTH)), by = .(CHROM,POS)]
  return(average_depth_dt)
}

chipfiles <- c("chipseqdata/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
               "chipseqdata/GSM4668651_Col0_rep2_H3K4me1_ChIP_unique_reads.bedGraph.gz")

infiles<-c("data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz", 
           "data/GSM4668652_Col0_rep2_H3K4me1_Input_unique_reads.bedGraph.gz")

avg_depth_chip<-calculate_average_depth_perbp(chipfiles)
fwrite(avg_depth_chip,"data/avg_chipseq_depth_wchrom.csv")
avg_depth_input<-calculate_average_depth_perbp(infiles)
#fwrite(avg_depth_input,"data/avg_input_depth.csv")


chiptotal<-sum(avg_depth_chip$avg_depth)
inputtotal<-sum(avg_depth_input$avg_depth)

merged_dt<-merge(avg_depth_chip,avg_depth_input,by="bp",suffixes = c("_chip","_input"))
merged_dt[,normalized_depth:= log2(avg_depth_chip/chiptotal)-log2(avg_depth_input/inputtotal)]

normalchip<-merged_dt[,c(1,4)]
fwrite(normalchip,"data/normalized_depth_chipseq_wCHROM.csv")






chipfiles <- c("chipseqdata/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz",
               "chipseqdata/GSM4668651_Col0_rep2_H3K4me1_ChIP_unique_reads.bedGraph.gz")

infiles<-c("chipseqdata/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz", 
           "chipseqdata/GSM4668652_Col0_rep2_H3K4me1_Input_unique_reads.bedGraph.gz")

calculate_average_depth_perbp <- function(replicate_files) {
  
  results <- list()
  
  for(file in replicate_files){
    message("Processing file: ", file)
    # Read a chunk of the file without loading the whole file into memory
    bg <- fread(file, sep = "\t", header = FALSE)
    colnames(bg) <- c("CHROM","START","STOP","DEPTH")
    
    # Get rid of unwanted chromosomes
    bg <- bg[CHROM %in% c(1:5),]
    
    # Define the chunk size
    chunk_size <- 10000 # Adjust as per your machine's memory 
    num_chunks <- ceiling(nrow(bg) / chunk_size)
    
    # Split the data into chunks
    chunks <- split(bg, rep(1:num_chunks, each=chunk_size, length.out=nrow(bg)))
    
    list_of_tables <- list()
    # Initialize progress bar
    pb <- txtProgressBar(min = 0, max = num_chunks, style = 3)
    
    for(i in seq_along(chunks)){
      chunk <- chunks[[i]]
      chunk_tables <- lapply(1:nrow(chunk), function(j) {
        chr <- chunk$CHROM[j]
        start <- as.numeric(chunk$START[j]) + 1
        stop <- as.numeric(chunk$STOP[j])
        depth <- as.numeric(chunk$DEPTH[j])
        data.table(CHROM=chr, POS=start:stop, DEPTH=depth)
      })
      list_of_tables[[length(list_of_tables) + 1]] <- rbindlist(chunk_tables)
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    # Close progress bar
    close(pb)
    results[[length(results) + 1]] <- rbindlist(list_of_tables)
  }
  message("Combining Replicates")
  combined_dt <- rbindlist(results)
  
  # Calculate the average depth for each "bp"
  message("Calculate Average")
  average_depth_dt <- combined_dt[, .(avg_depth = mean(DEPTH)), by = .(CHROM,POS)]
  return(average_depth_dt)
}

#get the average between replicates and expand the files to per bp
avg_depth_chip<-calculate_average_depth_perbp(chipfiles)
avg_depth_input<-calculate_average_depth_perbp(infiles)

#calculate the sum total depth of each average file
chiptotal<-sum(avg_depth_chip$avg_depth)
inputtotal<-sum(avg_depth_input$avg_depth)

#merge the files and get a normalized depth for both
merged_dt<-merge(avg_depth_chip,avg_depth_input,by=c('CHROM','POS'),suffixes = c("_chip","_input"))
merged_dt[,normalized_depth:= log2(avg_depth_chip/chiptotal)-log2(avg_depth_input/inputtotal)]

normalchip<-merged_dt[,c(1:2,4)]
fwrite(normalchip,"data/normalized_depth_chipseq_wCHROM.csv")

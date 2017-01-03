library(optparse)
library(R.utils)
source("src/mcf.R")

#parse options
opt <- get_options_mcf()
save(opt,file = "opt.RData")

levels <- opt$levels
min_region_size <- opt$min_region_size
max_region_size <- opt$max_region_size
min_start_stop_score <- opt$min_start_stop_score
min_mean_score <- opt$min_mean_score
min_median_score <- opt$min_median_score
min_score <- opt$min_score
output.file <- opt$out
data_file <- opt$data_file
flood_level_file <- opt$flood_level_file

rgn.file <- paste(output.file,"bed", sep = ".")

if (length(flood_level_file) > 0) {
  levels <- as.vector(as.matrix(read.table(flood_level_file)))
  cat(paste("Number of flood levels is ", length(levels) - 1, "\n", sep = ""))
}

blocks <- MCFGetRegionBlocks(data_file)

append <- FALSE
append.rgn <- FALSE
for (i in 1: length(blocks)) {
  print(paste("Block number ", i, " out of ", length(blocks), sep = ""))
  block <- blocks[[i]]

  scores <- as.numeric(block[,4])
  start <- as.numeric(block[,2])
  chrom <- as.numeric(block[1,1])
  #just make sure the input gives scores for consecutive bases
  stopifnot(max(as.numeric(block[,2])) - min(as.numeric(block[,2])) + 1 == length(scores))

  result <- MCF(scores,start,chrom,levels, min_region_size, max_region_size,
             min_start_stop_score,
             min_mean_score = min_mean_score,
             min_median_score = min_median_score,
             min_score = min_score)
  block <- cbind(block,result$mcf)
  names(block)[6] <- "mcf"
  bed  <- result$bed

  #write output
  if (append) {
    write.table(block, file = output.file, sep = '\t', row.names = FALSE,
                append=TRUE,col.names=FALSE, quote = FALSE)
  } else {
    write.table(block, file = output.file, sep = '\t', row.names = FALSE,
                append=FALSE,col.names=TRUE, quote = FALSE)
    append <- TRUE
  }

  if (append.rgn) {
    if (!is.null(bed)) {
      write.table(bed, file = rgn.file, sep = '\t', row.names = FALSE,
                append=TRUE,col.names=FALSE, quote = FALSE)
    }
  } else {
    if (!is.null(bed)) {
      write.table(bed, file = rgn.file, sep = '\t', row.names = FALSE,
                append=FALSE,col.names=TRUE, quote = FALSE)
      append.rgn <- TRUE
    }
  }
}


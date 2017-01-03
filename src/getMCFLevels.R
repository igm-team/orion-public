library(optparse)
library(R.utils)
source("src/mcf.R")

#parse options
opt <- get_options_mcf()
save(opt,file = "opt.RData")
levels <- opt$levels

output.file <- opt$out
data_file <- opt$data_file

all.scores <- NULL
blocks <- MCFGetRegionBlocks(data_file)
for (i in 1: length(blocks)) {
  print(paste("Block number ", i, " out of ", length(blocks), sep = ""))
  block <- blocks[[i]]
  scores <- as.numeric(block[,4])
  all.scores <- c(all.scores,scores[!is.na(scores)])
}
f.levels <- MCFGetFloodLevels(all.scores,levels)
write(f.levels, file = output.file,ncolumns = 1)

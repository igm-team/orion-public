get_options_mcf <- function() {
  option_list <- list(
    make_option(c("-d", "--orion-raw-data"), type="character",  dest='data_file', default=NULL,
                help="Tab-delimited file(with header) with five columns.", metavar="string"),
    make_option(c("-l", "--levels"), type="integer",  dest='levels', default=100,
                help="Default number of flood levels", metavar="integer"),
    make_option(c("--min_region_size"), type="integer",  dest='min_region_size', default=10,
                help="Minumum region size", metavar="integer"),
    make_option(c("--max_region_size"), type="integer",  dest='max_region_size', default=1000,
                help="Maximum region size", metavar="integer"),
    make_option(c("-o", "--out"), type="character",  dest='out', default=NULL,
                help="output file name", metavar="string"),

    #optional parameers
    make_option(c("--min_start_stop_score"), type="double",  dest='min_start_stop_score', default=NULL,
                help="Minimum start/stop score", metavar="double"),
    make_option(c("--min_mean_score"), type="double",  dest='min_mean_score', default=NULL,
                help="Minimum mean score of a region", metavar="double"),
    make_option(c("--min_median_score"), type="double",  dest='min_median_score', default=NULL,
                help="Minimum median score of a region", metavar="double"),
    make_option(c("--min_score"), type="double",  dest='min_score', default=NULL,
                help="Minimum score of a region", metavar="double"),
    make_option(c("--flood_level_file"), type="character",  dest='flood_level_file', default=NULL,
                help="Precomputed", metavar="string")
  )

  opt <- parse_args(OptionParser(option_list=option_list))

  if (is.null(opt$d)) {
    stop("Input raw score file is missing")
  }
  if (is.null(opt$o)) {
    stop("output file is missing")
  }
  return(opt)
}

MCF <- function(scores,start,chrom,levels, min_region_size, max_region_size,
                min_start_stop_score = NULL,min_mean_score = NULL,
                min_median_score = NULL, min_score = NULL) {

  if (length(levels) == 1) {   #generate levels is there is no precomputed levels
     floodlevels <- MCFGetFloodLevels(scores, levels)
     #write(floodlevels, file = "levels.txt")
  } else { #use precomputed level
    floodlevels <- levels
    print("using precomputed flood levels")
  }

  scores[is.na(scores)] <- floodlevels[1]
  indicators <- vector(length = length(scores))

  #reorder the flooding levels
  floodlevels <- sort(floodlevels,decreasing = TRUE)
  for (i in 1: (length(floodlevels)-1)) {
    current_flood_level <- floodlevels[i]
    curent_flood_status <- scores >= current_flood_level
    RLE <- rle(scores >= current_flood_level)
    ends <- c(0,cumsum(RLE$lengths)) #ends of intervals

    regions_under_size_control <- which(RLE$values & RLE$lengths>=min_region_size & RLE$lengths<=max_region_size)
    if (length(regions_under_size_control) > 0) {
       regions.start <- ends[regions_under_size_control] +1
       regions.end <- ends[regions_under_size_control+1]
       #sanity check, to be removed after testing runs
       stopifnot(all(regions.end - regions.start + 1 - RLE$lengths[regions_under_size_control] == 0))

       #to be implemented in C eventually
       for (j in 1: length(regions.start)) {

         if (!is.null(min_start_stop_score)) {
           if (scores[regions.start[j]] < min_start_stop_score || regions.end[j]  < min_start_stop_score) {
             next
           }
         }

         if (!is.null(min_mean_score)) {
           mean_score <- mean(scores[regions.start[j]:regions.end[j]])
           if (mean_score < min_mean_score) {
             next
           }
         }

         if (!is.null(min_median_score)) {
           median_score <- median(scores[regions.start[j]:regions.end[j]])
           if (median_score < min_median_score) {
             next
           }
         }

         if (!is.null(min_score)) {
           m_score <- min(scores[regions.start[j]:regions.end[j]])
           if (m_score < min_score) {
             next
           }
         }

         indicators[regions.start[j]:regions.end[j]] <- TRUE
       }
    }
  }
  result <- list()

  RLE <- rle(indicators)
  ends <- c(0,cumsum(RLE$lengths)) #ends of intervals
  regions_true <- which(RLE$values)
  if (length(regions_true) > 0) {
    START <- start[ends[regions_true] + 1]
    END <- start[ends[regions_true+1]] +1 #the right most base is not included
    bed <- as.data.frame(cbind(chrom,START,END))
    names(bed) <- c("#CHROM", "START", "END")
    result$bed <- bed
  } else {
    result$bed <- NULL
  }
  result$mcf <- indicators
  result
}

#read a bed file that have all orion scores
MCFGetRegionBlocks <- function(bed.file, blocksize = 1000000) {

  #get the total number of lines from the bed file
  total.lines <- countLines(bed.file)
  total.lines <- total.lines - 1 #remove header line

  con <- file(bed.file)
  open(con)

  header <- scan(con, what = character(),nlines =1, quiet = TRUE)
  stopifnot(length(header)==5) #make sure there are 5 columns in each row

  linecount <- 0

  blocks <- list()
  blockID <- 0

  #step 1: read raw blocks
  while (linecount < total.lines) {
    #read in data in blocks
    block <- scan(con, what = list("character","integer","integer","numeric","numeric"),sep = "\t",
                  nmax = blocksize, na.strings = "NA",quiet = TRUE)
    block <- as.data.frame.list(block,stringsAsFactors = FALSE)
    colnames(block) <- header
    linecount <- linecount + dim(block)[1]

    blockID <- blockID + 1
    blocks[[blockID]] <- block
    #print(blockID)
  }
  close(con)

  #step 2: addjust blocks
  if (length(blocks) > 1) {
    for (i in 1:(length(blocks)-1)) {
      if (!is.na(blocks[[i]][dim(blocks[[i]])[1],4])) { #a candidate block that needs to be fixed
        RLE <- rle(is.na(blocks[[i]][,4]))
        ends <- c(0,cumsum(RLE$lengths)) #ends of intervals
        last.interval.start <- ends[length(ends) -1] +1
        last.interval.end <- ends[length(ends)]
        blocks[[i+1]] <- rbind(blocks[[i]][last.interval.start:last.interval.end,],blocks[[i+1]])
        blocks[[i]] <- blocks[[i]][1:(last.interval.start-1),]
      }
    }
  }
  blocks
}

#give the number of levels, and considering not knowing the underline disctribtuion,
#the best guess is to use percentle
MCFGetFloodLevels <- function(scores, levels) {
  valid.scores <- scores[!is.na(scores)]
  flood_levels <- unique(quantile(valid.scores,(0:(levels-1)) / levels))
  #100 is just a number here, in fact, any possitive number will do
  flood_levels <- c(min(flood_levels) - 100,flood_levels)
  flood_levels
}

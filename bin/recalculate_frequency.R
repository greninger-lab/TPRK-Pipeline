# Recalculates relative frequency given a pre-filtered file (i.e. allreads_filtered.csv).
# This step is done in the tprk_pipeline.py. recalculate_frequency.R will replace the file.

suppressMessages(suppressWarnings(require("optparse")))
suppressMessages(suppressWarnings(require(vegan)))

option_list <- list(
  make_option(c("-d", "--directory"), type = "character", default = NULL, help = "Specify working directory", metavar = "character"),
  make_option(c("-f", "--filename"),
    type = "character", default = NULL, help = "Specify file to recalculate frequency.",
    metavar = "character"
  ),
  make_option(c("-m", "--metadata"), type = "character", default = NULL, help = "Specify metadata", metavar = "character"),
  make_option(c("-i", "--illumina_check"),
    type = "character", default = FALSE, help = "Specify if these files are only PacBio.",
    metavar = "character", action = "store_true"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

path <- "./"
filename <- (opt$filename)

#####

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding line of path
## but remember to recomment the lines if you want to run the script automatically in the pipeline.
## path refers to the folder your metadata.csv and sequencing files (.fastq) are.
## filename refers to the pre-filtered file (i.e. allreads_filtered.csv) that you want to recalculate frequencies for.

# path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/testing2"
# filename <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/testing2/allreads_filtered.csv"

## This script can also be run from the command line.
## Usage: rscript \path\to\recalculate_frequency.R -d [path] -f [file_name]

#####
allreads_filtered <- read.table(filename, sep = ",", header = TRUE)
# Grabs the actual number of samples.
numsamples <- (length(colnames(allreads_filtered)) - 2) / 4
metadata <- read.table(opt$metadata, sep = ",", header = TRUE)
sample_names <- c(as.character(metadata$SampleName))

allreads_filtered1 <- allreads_filtered

# Allreads_filtered.csv has both PacBio and Illumina data, so must run through twice.
if (grepl("allreads_filtered.csv", filename)) {
  numsamples <- numsamples * 2
}

# Loops through samples and recalculates frequency for each sample (column) for each region.
for (sample in c(1:(numsamples))) {
  freqcol <- (sample * 2) + 1
  countcol <- (sample * 2) + 2
  for (region in unique(allreads_filtered1$Region)) {
    allreads_filtered1[which(allreads_filtered1[, 1] == region), freqcol] <-
      allreads_filtered1[which(allreads_filtered1[, 1] == region), countcol] / sum(allreads_filtered1[which(allreads_filtered1[, 1] == region), countcol], na.rm = TRUE) * 100
  }
  allreads_filtered1[, freqcol] <- trimws(format(allreads_filtered1[, freqcol], digits = 4, nsmall = 4))
}

# Rewrites the file with the newly calculated relative frequencies.
write.csv(allreads_filtered1, file = filename, row.names = FALSE, quote = FALSE)

# Also makes Amin's requested CSV
if (grepl("allreads_filtered.csv", filename)) {
  amin_csv <- t(allreads_filtered1)
  variable_region_lists <- split(allreads_filtered1, f = allreads_filtered1$Region)
  for (i in 1:length(variable_region_lists)) {
    amin_csv <- t(variable_region_lists[[i]])
    colnames(amin_csv) <- NULL
    to_remove <- "Region"
    for (row in 1:nrow(amin_csv)) {
      row_name <- rownames(amin_csv)[row]
      if (grepl("_RelativeFreq", row_name)) {
        to_remove <- c(to_remove, row_name)
      }
    }
    removedamin_csv <- amin_csv[!(row.names(amin_csv) %in% to_remove), ]
    path <- gsub("[\r\n]", "", path)
    dir.create(file.path(path, "Diversity"), showWarnings = FALSE)
    write.table(removedamin_csv, file = paste(path, "/Diversity/V", i, ".csv", sep = ""), sep = ",", quote = FALSE, col.names = FALSE)
  }
}
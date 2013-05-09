library(stringr)
library(GenomicRanges)
# Setup parameters
#
# Version of Human Genome references. Current version is hg19. (options: 18 or 19)
hg_v <- 19


# Setting working directory
setwd("cnvdata")
# regex for hg version
hg19 <- paste0(".+_hg",hg_v,".*")
# regex for TCGA barcode
pattern <- "*TCGA-([A-Z0-9]{2})-([A-Z0-9]{4})-((0[1-9]|[1-2][0-9])[A-Z])-((0[1-9]|[1-9][0-9])[DGHRTWX])-([A-Z0-9]{4})-([A-Z0-9]{2})*"

files <- list.files(pattern=paste0(pattern, hg19), recursive=TRUE)

# Merge all files into a dataframe
print(paste0("Total files found: ", length(files)))
df <- NULL
for (i in seq(files)) {
    print(paste0("Reading file: ", i, "/", length(files)))
    f <- files[i]

    if (is.null(df)) {
        df <- read.table(f, sep ="\t", header=TRUE)
    } else {
        temp <- read.table(f, sep ='\t', header=TRUE)
        df <- rbind(df, temp)
    }
}

# Create a collection of genomic features from the imported data
gr <- GRanges(seqnames = df$chromosome,
              ranges = IRanges(start = df$start, end = df$stop),
              GC = df$seg.mean)
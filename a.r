library(stringr)
library(GenomicRanges)
# Setup parameters
#
# Version of Human Genome references. Current version is hg19. (options: 18 or 19)
hg_v <- 19

# Setting working directory
# setwd("cnvdata")

# regex for hg version with cnv and nocnv
regex.all <- paste0(".+hg",hg_v,".*")
# regex for hg version with cnv (exclude files that contain "nocnv_hg" on its name)
regex.cnv <- paste0(".+analysis.hg",hg_v,".*")
# regex for hg version nocnv (include only files that contain "nocnv_hg" on its name)
regex.nocnv <- paste0(".+analysis.nocnv_hg",hg_v,".*")
# regex for TCGA barcode
pattern <- "*(TCGA)-([A-Z0-9]{2})-([A-Z0-9]{4})-(0[1-9]|[1-2][0-9])([A-Z])-(0[1-9]|[1-9][0-9])([DGHRTWX])-([A-Z0-9]{4})-([A-Z0-9]{2})*"

# merge
source("merge.R")
# all segments
files <- list.files(pattern=paste0(pattern, regex.all), recursive=TRUE)
seg.all <- mergeFiles(files)
save(seg.all, file="data/lgg.all.rda")
# cnv
files <- list.files(pattern=paste0(pattern, regex.cnv), recursive=TRUE)
seg.cnv <- mergeFiles(files)
save(seg.cnv, file="data/lgg.cnv.rda")
# nocnv
files <- list.files(pattern=paste0(pattern, regex.nocnv), recursive=TRUE)
seg.nocnv <- mergeFiles(files)
save(seg.nocnv, file="data/lgg.nocnv.rda")

# # remove chromosome X and Y
# df.c <- subset(df, ((df$"chromosome"!="X")&(df$"chromosome"!="Y")))
# save(df.c, file="data/noxy.lgg.cnv.rda")

# Manifest file
manifest <- data.frame(stringsAsFactors=FALSE)

# by file
# barcodes <- c()
# for (i in seq(files)) {
#     f <- files[i]

#     barcodes <- append(barcodes,c(str_match(f, pattern)[,1])) 
# }

# by array with barcode
barcodes <- unique(seg.all$barcode)

# generate Manifest
source("barcode.r")
manifest <- getMetadata(barcodes)
write.table(manifest, file="MANIFEST.txt", row.names=F)
save(manifest, file="data/manifest.rda")

# choose sample type (NORMAL|TUMOR|UNKNOWN) as interest query
tumor <- levels(droplevels(subset(manifest, manifest$"sample.type"=="TUMOR")$"barcode"))
normal <- levels(droplevels(subset(manifest, manifest$"sample.type"=="NORMAL")$"barcode"))

# interest subset of all
seg.all.tumor <- seg.all[seg.all$barcode %in% tumor,]
save(seg.all.tumor, file="data/tumor.lgg.all.rda")
seg.all.normal <- seg.all[seg.all$barcode %in% normal,]
save(seg.all.normal, file="data/normal.lgg.all.rda")
# interest subset of cnv
seg.cnv.tumor <- seg.cnv[seg.cnv$barcode %in% tumor,]
save(seg.cnv.tumor, file="data/tumor.lgg.cnv.rda")
seg.cnv.normal <- seg.cnv[seg.cnv$barcode %in% normal,]
save(seg.cnv.normal, file="data/normal.lgg.cnv.rda")
# interest subset of nocnv
seg.nocnv.tumor <- seg.nocnv[seg.nocnv$barcode %in% tumor,]
save(seg.nocnv.tumor, file="data/tumor.lgg.nocnv.rda")
seg.nocnv.normal <- seg.nocnv[seg.nocnv$barcode %in% normal,]
save(seg.nocnv.normal, file="data/normal.lgg.nocnv.rda")

# histogram of all
tiff("histogram/all.tiff")
hist(seg.all$"seg.mean", breaks=100, main="Histogram of Segment Mean (CNV and NOCNV)", xlab="Segment Mean")
graphics.off()
tiff("histogram/all.tumor.tiff")
hist(seg.all.tumor$"seg.mean", breaks=100, main="Histogram of Segment Mean (CNV and NOCNV) - TUMOR", xlab="Segment Mean")
graphics.off()
tiff("histogram/all.normal.tiff")
hist(seg.all.normal$"seg.mean", breaks=100, main="Histogram of Segment Mean (CNV and NOCNV) - NORMAL", xlab="Segment Mean")
graphics.off()
# histogram of cnv
tiff("histogram/cnv.tiff")
hist(seg.cnv$"seg.mean", breaks=100, main="Histogram of Segment Mean (CNV)", xlab="Segment Mean")
graphics.off()
tiff("histogram/cnv.tumor.tiff")
hist(seg.cnv.tumor$"seg.mean", breaks=100, main="Histogram of Segment Mean (CNV) - TUMOR", xlab="Segment Mean")
graphics.off()
tiff("histogram/cnv.normal.tiff")
hist(seg.cnv.normal$"seg.mean", breaks=100, main="Histogram of Segment Mean (CNV) - NORMAL", xlab="Segment Mean")
graphics.off()
# histogram of nocnv
tiff("histogram/nocnv.tiff")
hist(seg.nocnv$"seg.mean", breaks=100, main="Histogram of Segment Mean (NOCNV)", xlab="Segment Mean")
graphics.off()
tiff("histogram/nocnv.tumor.tiff")
hist(seg.nocnv.tumor$"seg.mean", breaks=100, main="Histogram of Segment Mean (NOCNV) - TUMOR", xlab="Segment Mean")
graphics.off()
tiff("histogram/nocnv.normal.tiff")
hist(seg.nocnv.normal$"seg.mean", breaks=100, main="Histogram of Segment Mean (NOCNV) - NORMAL", xlab="Segment Mean")
graphics.off()


# Create a collection of genomic features from the imported data
gr <- GRanges(seqnames = df$chromosome,
              ranges = IRanges(start = df$start, end = df$stop),
              GC = df$seg.mean)

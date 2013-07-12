library(iCluster)
library(GenomicRanges)

load("data/manifest.rda")
load("data/noxy.tumor.lgg.cnv.rda")

df.chr01 <- subset(df.q, chromosome==1)

fit=iCluster(as.list(df.chr01$seg.mean), k=4, lambda=0.2, max.iter=50)
plot.iCluster(fit=fit, label=rownames(datasets[[2]]))

gr <- GRanges(seqnames = df.chr01$chromosome,
              ranges = IRanges(start = df.chr01$start, end = df.chr01$stop, names=df.chr01$barcode),
              num.mark = df.chr01$num.mark,
              seg.mean = df.chr01$seg.mean)

CNA.object <- CNA(cbind(coriell$Coriell.05296),coriell$Chromosome,coriell$Position,data.type="logratio",sampleid="c05296")
CNA.object <- CNA(cbind(coriell$Coriell.05296),coriell$Chromosome,coriell$Position,data.type="logratio",sampleid="c05296")



source("http://bioconductor.org/biocLite.R")
biocLite("biomvRCNS")
library(biomvRCNS)
data('coriell', package='biomvRCNS')

# http://www.cbioportal.org/public-portal/index.do
# Select Cancer Study: Brain Lower Grade Glioma
# Putative copy-number alterations from GISTIC data
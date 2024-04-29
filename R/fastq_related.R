library(Biostrings)
library(GenomicAlignments)
library(Rsamtools)

raw_fa <- list.files("../rawdata", pattern = "\\.fasta(\\.gz)?", full.names = TRUE)

for (i in raw_fa){
    fa <- readDNAStringSet(i, format = "fasta")
    avg <- ceiling(mean(fa@ranges@width))
    std <- ceiling(sd(fa@ranges@width))
    max_ <- max(fa@ranges@width)
    min_ <- min(fa@ranges@width)
    cat(i, ":", avg, "\U00B1", std, paste0("(", min_, ", ", max_, ")"), "\n")
    duplicated_reads <- names(fa)[which(duplicated(names(fa)))]
    cat("Duplicated reads: ", duplicated_reads, "\n\n")
}

intergenic <- readDNAStringSet("D:/bcst/YAMADA/jklai/Araport11/Araport11_intergenic_20220504.fa")
w <- width(intergenic)
avg <- ceiling(mean(w))
std <- ceiling(sd(w))
max_ <- max(w)
min_ <- min(w)
cat("Intergenic region : ", avg, "\U00B1", std)

###########################################################################################
# Check BAM
###########################################################################################
bam_file <- "../BAM_mm2/hq_transcripts-4.bam"

bam <- readGAlignments(bam_file)
seqlevels(bam)

which_ <- GRanges(c(
    "Chr1:1-20",
    "Chr2:10-100",
    "Chr2:50-150"
))
what_ <- c("rname", "qwidth", "seq")
param <- ScanBamParam(which = which_, what = what_)

bam <- scanBam(bam_file, param = param)
bam <- scanBam(file = bam_file)
names(bam[[1]])

bam_ <- bam[[1]]
bam_$mrnm

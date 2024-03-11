source("./load_packages.R")
load_packages(BiocGenerics, GenomicRanges, GenomicFeatures, TxDb.Athaliana.BioMart.plantsmart51, tidyverse)

get_intergenic <- function(
        fasta_dir = "D:/bcst/YAMADA/jklai/Araport11/Araport11_intergenic_20220504.fa",
        output_gtf = TRUE,
        output_dir = NULL
){
    fa <- readLines(fasta_dir)
    fa <- fa[unlist(gregexpr(">", fa)) == 1]
    
    df0 <- data.frame(
        seqname = gsub(pattern = ".*\\|\\W(.+):.*", replacement = "\\1", x = fa),
        source = "Araport11",
        feature = "intergenic",
        start = gsub(pattern = ".*Chr\\w:(.+)-(.+)\\WFORWARD.*", replacement = "\\1", x = fa),
        end = gsub(pattern = ".*Chr\\w:(.+)-(.+)\\WFOR.*", replacement = "\\2", x = fa),
        score = ".",
        strand = "*",
        frame = ".",
        attribute = gsub(pattern = ">(.+) \\|.*", replacement = "\\1", x = fa)
    )
    
    res <- df0 %>% 
        mutate(across(c(start, end), as.numeric)) %>% 
        # mutate(attribute = paste0("\"", attribute, "\"")) %>% 
        mutate(attribute = paste0("transcript_id \"", attribute, "\"; gene_id \"", attribute, "\";"))
    
    if (output_gtf & !is.null(output_dir)){
        write.table(res, file = output_dir, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    }else{
        return(res)
    }
    
}


# genic <- GenomicFeatures::genes(TxDb.Athaliana.BioMart.plantsmart51)
# gene_id <- genic$gene_id
# genic <- GenomicRanges::reduce(genic, ignore.strand = TRUE)
# intergenic <- GenomicRanges::gaps(genic)
# intergenic <- intergenic[BiocGenerics::strand(intergenic) == "*"]
# 
# 
# df0 <- as.data.frame(intergenic) %>% 
#     mutate(
#         seqnames = paste0("Chr", seqnames),
#         seqnames = case_when(seqnames == "ChrMt" ~ "ChrM",
#                              seqnames == "ChrPt" ~ "ChrC",
#                              TRUE ~ seqnames)
#     )
# 
# gtf <- data.frame(
#     seqname = df0$seqnames,
#     source = "TxDb.Athaliana.BioMart.plantsmart51",
#     feature = "intergenic",
#     start = df0$start,
#     end = df0$end,
#     score = ".",
#     strand = df0$strand,
#     frame = ".",
#     attribute = ""
# )
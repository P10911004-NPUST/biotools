read_gtf <- function(gtf_file){
    if (!require(readr)) install.packages("readr")
    library(readr)
    gtf <- read_tsv(gtf_file, comment = "#", col_names = FALSE)
    colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
    return(gtf)
}

get_novel_isoform <- function(gtf_file, feature = "all", known_isoform_attr = "ref_gene_id"){
    if (!require(dplyr)) install.packages("dplyr")
    library(dplyr)
    
    novel_isoform <- read_gtf(gtf_file) %>% 
        dplyr::filter(!str_detect(attributes, !!known_isoform_attr))
    
    if (tolower(feature) != "all"){
        novel_isoform <- novel_isoform %>% 
            dplyr::filter(feature == !!feature)
    }
    
    return(novel_isoform)
}

write_gtf <- function(x, file){
    write_tsv(x, file, col_names = FALSE, quote = "none", escape = "none")
    
}
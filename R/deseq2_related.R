suppressMessages({
    if (!require(readr)) install.packages("readr")
    library(readr)
    if (!require(stringr)) install.packages("stringr")
    library(stringr)
    if (!require(dplyr)) install.packages("dplyr")
    library(dplyr)
    if (!require(DESeq2)) install.packages("DESeq2")
    library(DESeq2)
    # options(scipen = 100)
    
    source("https://github.com/P10911004-NPUST/biotools/blob/main/R/arabidopsis.R?raw=true")
})


get_norm_count <- function(dds_object, returnData = TRUE, ...){
    dds <- DESeq2::estimateSizeFactors(dds_object)
    dds_counts <- DESeq2::counts(dds, normalized = TRUE)
    return(dds_counts)
}


get_result <- function(
        object = NULL,
        contrast = NULL,
        pvalueCutoff = 0.05,
        lfcCutOff = 0,
        entrez_id = FALSE
){
    res <- DESeq2::results(
        object = object,
        contrast = contrast,
        alpha = pvalueCutoff
    ) %>% 
        as.data.frame() %>% 
        tidyr::drop_na() %>%
        tibble::rownames_to_column("gene_id") %>%
        dplyr::filter(padj < pvalueCutoff) %>% 
        dplyr::mutate(
            reg_pos = case_when(
                log2FoldChange > abs(lfcCutOff) ~ "up",
                log2FoldChange < -abs(lfcCutOff) ~ "down",
            )
        ) %>%
        tidyr::drop_na()
    
    if (entrez_id) res <- res %>% dplyr::mutate(entrez_id = tair2entrez(gene_id))
    
    return(res)
}















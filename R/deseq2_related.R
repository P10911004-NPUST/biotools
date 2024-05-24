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
})


get_norm_count <- function(dds_object, returnData = TRUE, ...){
    dds <- DESeq2::estimateSizeFactors(dds_object)
    dds_counts <- DESeq2::counts(dds, normalized = TRUE)
    return(dds_counts)
}


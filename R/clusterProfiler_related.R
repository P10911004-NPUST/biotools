suppressMessages({
    if (!require(dplyr)) BiocManager::install("dplyr")
    if (!require(magrittr)) BiocManager::install("magrittr")
    if (!require(BiocManager)) BiocManager::install("BiocManager")
    if (!require(GO.db)) BiocManager::install("GO.db")
    if (!require(org.At.tair.db)) BiocManager::install("org.At.tair.db")
    if (!require(AnnotationDbi)) BiocManager::install("AnnotationDbi")
    if (!require(clusterProfiler)) BiocManager::install("clusterProfiler")
    if (!require(biomartr)) BiocManager::install("biomartr")
    
    library(dplyr)
    library(magrittr)
    library(GO.db)
    library(org.At.tair.db)
    library(AnnotationDbi)
    library(clusterProfiler)
    # options(scipen = 100)  # large value will prevent scientific notation
})


tair2entrez <- function(tair_id, output_data_types = c("vector", "list", "data.frame")){
    data_type <- match.arg(output_data_types)
    df0 <- clusterProfiler::bitr(
        geneID = tair_id,
        OrgDb = org.At.tair.db,
        fromType = "TAIR",
        toType = "ENTREZID",
        drop = FALSE
    )
    
    # Output a vector
    if (data_type == "vector") return(df0 %>% dplyr::pull(ENTREZID, TAIR))
    
    # Output a list
    if (data_type == "vector") return(as.list(df0))
    
    # Output a data frame
    if (data_type == "data.frame") return(df0)
}




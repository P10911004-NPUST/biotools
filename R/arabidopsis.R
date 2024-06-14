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

entrez2tair <- function(entrez_id, output_data_types = c("vector", "list", "data.frame")){
    data_type <- match.arg(output_data_types)
    
    df0 <- clusterProfiler::bitr(
        geneID = entrez_id,
        OrgDb = org.At.tair.db,
        fromType = "ENTREZID",
        toType = "TAIR",
        drop = FALSE
    )
    
    # Output a vector
    if (data_type == "vector") return(df0 %>% dplyr::pull(TAIR, ENTREZID))
    
    # Output a list
    if (data_type == "vector") return(as.list(df0))
    
    # Output a data frame
    if (data_type == "data.frame") return(df0)
}


# ==== Over Representation Analysis (ORA) ====
tair_enrichGO <- function(
        gene_id, 
        ontology = c("ALL", "BP", "CC", "MF"), 
        qvalueCutoff = 0.05,  # default: 0.2
        pvalueCutoff = 0.01,  # default: 0.05
        simplify = TRUE,
        return_as_dataframe = TRUE
){
    df0 <- "NA"
    
    res <- clusterProfiler::enrichGO(
        gene = gene_id,
        # universe = universe_tair_id,
        OrgDb = org.At.tair.db,
        keyType = "TAIR",
        ont = match.arg(ontology),  # c("ALL", "BP", "CC", "MF")
        qvalueCutoff = qvalueCutoff,  # default: 0.2
        pvalueCutoff = pvalueCutoff,  # default: 0.05
        pAdjustMethod = "BH",  # default: "BH"
        readable = FALSE,  # TRUE: show Gene Symbol, FALSE: show TAIR_id
        pool = FALSE  # FALSE: Show ontology separately
    )
    
    if (nrow(res@result > 0)){
        if (simplify) res <- res %>% clusterProfiler::simplify()
        
        df0 <- res@result %>% 
            mutate(
                bg_ratio = Count / as.numeric(str_split_i(BgRatio, "/", 2)),
                gene_ratio = Count / as.numeric(str_split_i(GeneRatio, "/", 2)),
                fold_enrich = gene_ratio / bg_ratio,
                rich_factor = Count / as.numeric(sub("/\\d+", "", BgRatio)),
            )
    }
    
    ifelse(return_as_dataframe, return(df0), return(res))
}


# ==== Gene Set Enrichment Analysis (GSEA) ====
tair_gseGO <- function(
        gene_named_list, 
        minGSSize = 5,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        eps = 0,
        simplify = TRUE
){
    gene_named_list <- sort(gene_named_list, decreasing = TRUE)
    
    res <- clusterProfiler::gseGO(
        geneList = gene_named_list,
        OrgDb = org.At.tair.db,
        keyType = "ENTREZID",
        ont = "ALL",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff,  # set to 1 in order to obtain all possible terms 
        pAdjustMethod = "BH",
        eps = eps,
        verbose = FALSE,  # FALSE: Show ontology separately
        seed = TRUE
    )
    
    if (nrow(res) < 1) return(NULL)
    
    if (simplify) res <- res %>% clusterProfiler::simplify()
    
    return(as.data.frame(res))
}


# ==== KEGG (KEGG) ====
tair_gseKEGG <- function(
        gene_named_list,
        pvalueCutoff = 0.05,
        minGSSize = 3,
        maxGSSize = 500
){
    gene_named_list <- sort(gene_named_list, decreasing = TRUE)
    
    res <- clusterProfiler::gseKEGG(
        geneList = gene_named_list,
        organism = "ath",
        use_internal_data = FALSE,
        pvalueCutoff = pvalueCutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        verbose = FALSE
    )
    
    return(res@result)
}


# ==== Compare cluster ====
tair_compare_cluster <- function(
        gene_named_list, 
        data = NULL,
        apply_method = c("ora", "gsea", "kegg"), 
        qvalueCutoff = 0.05,  # default: 0.2
        pvalueCutoff = 0.01,  # default: 0.05
        minGSSize = 5,
        maxGSSize = 500,
        eps = 0,
        simplify = TRUE
){
    apply_method <- match.arg(apply_method)
    
    if (tolower(apply_method) == "ora"){
        res <- clusterProfiler::compareCluster(
            geneClusters = entrez_id ~ group,
            OrgDb = org.At.tair.db,
            fun = enrichGO,
            data = cluster_df
        )
    }
    
    if (tolower(apply_method) == "gsea"){
        res <- clusterProfiler::gseGO(
            geneList = gene_named_list,
            OrgDb = org.At.tair.db,
            keyType = "ENTREZID",
            ont = "ALL",
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            pvalueCutoff = pvalueCutoff,  # set to 1 in order to obtain all possible terms 
            pAdjustMethod = "BH",
            eps = eps,
            verbose = FALSE,  # FALSE: Show ontology separately
            seed = TRUE
        )
    }
    
    if (tolower(apply_method) == "kegg"){
        res <- NULL
        print("Not yet")
    }
    
    if (nrow(res) < 1) return(NULL)
    
    if (simplify) res %>% clusterProfiler::simplify()
    
    res <- as.data.frame(res)
    
    return(res)
}

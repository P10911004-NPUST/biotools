suppressMessages({
    if (!require(dplyr)) install.packages("dplyr")
    if (!require(magrittr)) install.packages("magrittr")
    if (!require(BiocManager)) install.packages("BiocManager")
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
    
    windowsFonts(
        palatino = windowsFont("Palatino Linotype"),
        times = windowsFont("Times New Roman"),
        arial = windowsFont("Arial")
    )
    
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
        qvalueCutoff = 0.2,  # default: 0.2       # consider using 0.05
        pvalueCutoff = 0.05,  # default: 0.05     # consider using 0.01
        simplify = TRUE,
        return_as_dataframe = TRUE,
        ...
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
        pool = FALSE,  # FALSE: Show ontology separately
        ...
    )
    
    if (nrow(res@result > 0)){
        if (simplify) res <- res %>% clusterProfiler::simplify()
        
        df0 <- res@result %>% 
            dplyr::mutate(
                bg_ratio = Count / as.numeric(str_split_i(BgRatio, "/", 2)),
                gene_ratio = Count / as.numeric(str_split_i(GeneRatio, "/", 2)),
                fold_enrich = gene_ratio / bg_ratio,
                rich_factor = Count / as.numeric(sub("/\\d+", "", BgRatio)),
            ) %>% 
            dplyr::arrange(ONTOLOGY, desc(rich_factor))
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


# GO dotplot ====
go_dotplot <- function(
        data, 
        x, 
        y, 
        group_var = NULL, 
        size_var, 
        color_var,
        shape_var = NULL,
        text_wrap = 25,
        text_size = 22,
        legend_text_size = 20
){
    x <- deparse(substitute(x))
    y <- deparse(substitute(y))
    group_var <- deparse(substitute(group_var))
    size_var <- deparse(substitute(size_var))
    color_var <- deparse(substitute(color_var))
    shape_var <- deparse(substitute(shape_var))
    
    df0 <- data.frame(
        x = data[[x]],
        y = data[[y]],
        size_var = data[[size_var]],
        color_var = data[[color_var]]
    )
    
    p1 <- ggplot(
        data = df0, 
        mapping = aes(
            x = x, 
            y = fct_reorder(y, x)
        )
    ) +
        theme_bw() +
        labs(
            x = sym(x),
            y = sym(y),
            size = sym(size_var),
            color = sym(color_var),
        ) +
        geom_point(aes(size = size_var, color = color_var)) +
        geom_segment(
            mapping = aes(xend = 0, yend = y, color = color_var, linewidth = size_var), 
            show.legend = FALSE
        ) +
        scale_y_discrete(labels = function(x) str_wrap(x, text_wrap)) +
        scale_color_gradientn(
            colours = c("#a8dadc", "#457b9d", "#1d3557"),
            trans = "log10",
            guide = guide_colorbar(
                reverse = TRUE,
                order = 1,
                theme = theme(
                    legend.key.size = grid::unit(1, "cm")
                )
            ),
        ) +
        scale_size_continuous(
            range = c(5, 12),
            guide = guide_legend(
                order = 2,
                theme = theme(
                    legend.key.size = grid::unit(0.7, "cm")
                )
            )
        ) +
        scale_shape(
            guide = guide_legend(
                order = 3,
                override.aes = list(size = 5)
            )
        ) +
        theme(
            text = element_text(family = "arial", face = "bold", size = text_size),
            axis.title.x = element_markdown(margin = ggplot2::margin(t = 7)),
            axis.text.x.bottom = element_text(size = text_size),
            axis.title.y = element_blank(),
            axis.text.y = element_text(family = "arial", face = "bold", size = text_size),
            legend.title = element_markdown(family = "arial", face = "bold", size = legend_text_size),
            legend.text = element_text(family = "arial", face = "plain", size = legend_text_size - 1),
            legend.background = element_blank(),
            legend.position = "right"
        )
    
    if (group_var != "NULL"){
        p1$data$group_var <- data[[group_var]]
        p1 <- p1 +
            facet_wrap(~ group_var, nrow = 1, scales = "free")
    }
    
    if (shape_var != "NULL"){
        p1$layers[[1]] <- NULL
        p1$data$shape_var <- data[[shape_var]]
        p1 <- p1 +
            geom_point(aes(size = size_var, color = color_var, shape = shape_var)) +
            labs(
                x = sym(x),
                y = sym(y),
                size = sym(size_var),
                color = sym(color_var),
                shape = sym(shape_var)
            )
    }
    
    return(p1)
}

# go_dotplot(df0_up_set1, rich_factor, Description, gene_ratio, p.adjust)


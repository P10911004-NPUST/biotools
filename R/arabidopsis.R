suppressMessages({
    if (!require(dplyr)) install.packages("dplyr")
    if (!require(ggplot2)) install.packages("ggplot2")
    if (!require(magrittr)) install.packages("magrittr")
    if (!require(BiocManager)) install.packages("BiocManager")
    if (!require(ggvenn)) install.packages("ggvenn")
    if (!require(ggtext)) install.packages("ggtext")
    if (!require(GO.db)) BiocManager::install("GO.db")
    if (!require(org.At.tair.db)) BiocManager::install("org.At.tair.db")
    if (!require(AnnotationDbi)) BiocManager::install("AnnotationDbi")
    if (!require(clusterProfiler)) BiocManager::install("clusterProfiler")
    if (!require(biomartr)) BiocManager::install("biomartr")
    if (!require(KEGGREST)) BiocManager::install("KEGGREST")
    if (!require(pathview)) BiocManager::install("pathview")
    
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


compound2kegg <- function(compound_name, compound_mass = NULL) {
    cn <- gsub("N-|N,|\U03B1|\U03B2|\U03B3|\U03B4", "", compound_name)  # alpha, beta, gamma, delta
    # cn <- gsub("\U03B1", "", compound_name)  # alpha
    # cn <- gsub("\U03B2", "", cn)  # beta
    # cn <- gsub("\U03B3", "", cn)  # gamma
    # cn <- gsub("\U03B4", "", cn)  # delta
    cn <- gsub("/|:|-|,|<|>|\\(|\\)", " ", cn)
    cn <- strsplit(cn, "\\s")[[1]]
    cn <- cn[!cn %in% c("", " ")]
    cn <- cn[!base::duplicated(cn)]
    
    ret_id <- KEGGREST::keggFind("compound", cn)
    ret_id <- names(ret_id)
    ret_mass <- ret_id  # if ret_id has more than one elements, will be changed as below
    
    if (length(ret_id) > 1 & !is.null(compound_mass)) {
        compound_mass <- round(compound_mass, 2)
        ret_mass <- KEGGREST::keggFind("compound", compound_mass, option = "exact_mass")
        ret_mass <- names(ret_mass)
    }
    
    ret <- intersect(ret_id, ret_mass)
    ret <- vapply(ret, function(x) gsub("cpd:", "", x), character(1), USE.NAMES = FALSE)
    ret <- paste(ret, collapse = "/")
    if (ret == "") ret <- NA_character_
    
    return(ret)
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



# Venn plot ====
venn_label_position <- list(
    x = c(
        "A" = -0.85, "B" = 0, "C" = 0.85, 
        "D" = -1.1, "E" = 1.1, 
        "F" = -1.5, "G" = 1.5, 
        "H" = -0.5, "I" = 0.5, 
        "J" = -0.85, "K" = 0, "L" = 0.85, 
        "M" = -0.35, "N" = 0.35, "O" = 0
    ),
    # A     B     C     D     E     F     G     
    y = c(
        "A" = 0.95, "B" = 0.65, "C" = 0.95, 
        "D" = 0.55, "E" = 0.55, 
        "F" = 0.40, "G" = 0.40, 
        "H" = 0.10, "I" = 0.10, 
        "J" = -0.60, "K" = -0.45, "L" = -0.60, 
        "M" = -0.90, "N" = -0.90, "O" = -1.17
    )
)

ggvenn4 <- function(
        data = list(), 
        
        show_set_names = TRUE,
        set_names = c(),
        set_names_text_size = 6,
        
        show_labels = TRUE,
        label_text_size = 5,
        
        text_size = 5,
        stroke_size = NA,
        ...
){
    v1 <- data[[1]]
    v2 <- data[[2]]
    v3 <- data[[3]]
    v4 <- data[[4]]
    # layer 1 = polygon
    # layer 2 = polygon stroke
    # layer 3 = set names
    # layer 4 = set data
    set_elements <- list(
        set_A = v2 %>% setdiff(v3) %>% setdiff(v1) %>% setdiff(v4),
        set_B = v2 %>% intersect(v3) %>% setdiff(v1) %>% setdiff(v4),
        set_C = v3 %>% setdiff(v2) %>% setdiff(v1) %>% setdiff(v4),
        set_D = v2 %>% intersect(v1) %>% setdiff(v3) %>% setdiff(v4),
        set_E = v3 %>% intersect(v4) %>% setdiff(v2) %>% setdiff(v1),
        set_F = v1 %>% setdiff(v2) %>% setdiff(v3) %>% setdiff(v4),
        set_G = v4 %>% setdiff(v3) %>% setdiff(v2) %>% setdiff(v1),
        set_H = v1 %>% intersect(v2) %>% intersect(v3) %>% setdiff(v4),
        set_I = v2 %>% intersect(v3) %>% intersect(v4) %>% setdiff(v1),
        set_J = v1 %>% intersect(v3) %>% setdiff(v2) %>% setdiff(v4),
        set_K = v1 %>% intersect(v2) %>% intersect(v3) %>% intersect(v4),
        set_L = v4 %>% intersect(v2) %>% setdiff(v3) %>% setdiff(v1),
        set_M = v1 %>% intersect(v4) %>% intersect(v3) %>% setdiff(v2),
        set_N = v4 %>% intersect(v2) %>% intersect(v1) %>% setdiff(v3),
        set_O = v1 %>% intersect(v4) %>% setdiff(v2) %>% setdiff(v3)
    )
    
    p0 <- ggvenn(
        data = data,
        text_size = text_size,
        stroke_size = stroke_size,
        ...
    ) +
        theme(
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            text = element_text(size = 22, face = "bold", family = "sans")
        )
    
    p0[["layers"]][[3]] <- NULL
    
    if (show_set_names){
        if (length(set_names) == 0) set_names <- names(data)
        p0 <- p0 +
            annotate(
                geom = "richtext",
                x = c(-1.5, -1, 1, 1.5),
                y = c(-1, 1.2, 1.2, -1),
                label = set_names,
                size = set_names_text_size,
                fontface = "bold",
                angle = c(310, 0, 0, 50),
                fill = NA,
                label.color = NA
            )
    }
    
    if (show_labels){
        p0 <- p0 +
            annotate(
                geom = "text",
                x = venn_label_position$x,
                y = venn_label_position$y,
                label = LETTERS[1:length(venn_label_position$x)],
                size = label_text_size,
                fontface = "bold"
            )
    }
    
    return(list(figure = p0, set_elements = set_elements))
}


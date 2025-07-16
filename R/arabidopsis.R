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


tair2entrez <- function(tair_id, drop_NA = FALSE) {
    df0 <- clusterProfiler::bitr(
        geneID = tair_id,
        OrgDb = org.At.tair.db,
        fromType = "TAIR",
        toType = "ENTREZID",
        drop = drop_NA
    )
    ret <- setNames(df0[["ENTREZID"]], df0[["TAIR"]])
    return(ret)
}

entrez2tair <- function(tair_id, drop_NA = FALSE) {
    df0 <- clusterProfiler::bitr(
        geneID = tair_id,
        OrgDb = org.At.tair.db,
        fromType = "ENTREZID",
        toType = "TAIR",
        drop = drop_NA
    )
    ret <- setNames(df0[["TAIR"]], df0[["ENTREZID"]])
    return(ret)
}

# entrez2tair <- function(entrez_id, output_data_types = c("vector", "list", "data.frame")){
#     data_type <- match.arg(output_data_types)
#     
#     df0 <- clusterProfiler::bitr(
#         geneID = entrez_id,
#         OrgDb = org.At.tair.db,
#         fromType = "ENTREZID",
#         toType = "TAIR",
#         drop = FALSE
#     )
#     
#     # Output a vector
#     if (data_type == "vector") return(df0 %>% dplyr::pull(TAIR, ENTREZID))
#     
#     # Output a list
#     if (data_type == "vector") return(as.list(df0))
#     
#     # Output a data frame
#     if (data_type == "data.frame") return(df0)
# }


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
        pvalueCutoff = 0.05,  # default: 0.05     # consider using 0.01
        qvalueCutoff = 0.2,  # default: 0.2       # consider using 0.05
        simplify = TRUE,
        return_as_dataframe = TRUE,
        ...
){
    
    ontology <- match.arg(ontology)  # c("ALL", "BP", "CC", "MF")
    
    df0 <- NULL
    res <- NULL
    
    res <- clusterProfiler::enrichGO(
        gene = gene_id,
        # universe = universe_tair_id,
        OrgDb = org.At.tair.db,
        keyType = "TAIR",
        ont = ontology,  
        qvalueCutoff = qvalueCutoff,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = "BH",
        readable = FALSE,  # set to TRUE if wanna show Gene Symbol, otherwise show TAIR_id
        pool = FALSE,  # set to FALSE to show the ontology categories separately
        ...
    )
    
    # if res is NULL object, then nrow() will raise error
    # therefore, `||` is necessary, DO NOT change to `|`
    if (!(is.null(res) || nrow(res@result) == 0)){ 
        if (simplify) res <- clusterProfiler::simplify(res)
        
        df0 <- res@result
        if (ontology != "ALL") df0$ONTOLOGY <- ontology
        
        df0 <- df0 %>%
            dplyr::mutate(
                bg_ratio = parse_ratio(BgRatio),
                gene_ratio = parse_ratio(GeneRatio),
                gene_id = geneID
            ) %>%
            dplyr::relocate(ID, .after = last_col()) %>%
            dplyr::relocate(geneID, .after = ONTOLOGY) %>%
            dplyr::arrange(
                ONTOLOGY, 
                dplyr::desc(RichFactor), 
                dplyr::desc(FoldEnrichment),
                p.adjust,
                dplyr::desc(gene_ratio)
            ) %>% 
            dplyr::mutate(
                rich_factor = round(RichFactor, 2),
                fold_enrich = round(FoldEnrichment, 2)
            ) %>% 
            dplyr::select(
                ONTOLOGY, gene_id, Description, Count, 
                rich_factor, fold_enrich, gene_ratio, bg_ratio,
                p.adjust, qvalue, ID
            )
    }
    
    ifelse(return_as_dataframe, return(df0), return(res))
}


# ==== Gene Set Enrichment Analysis (GSEA) ====
tair_gseGO <- function(
        gene_named_vector, # gene_id as names, log2FoldChange as values
        minGSSize = 3,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        # eps = 0,
        simplify = TRUE
){
    gene_named_vector <- sort(gene_named_vector, decreasing = TRUE)
    
    res <- clusterProfiler::gseGO(
        geneList = gene_named_vector,
        OrgDb = org.At.tair.db,
        keyType = "ENTREZID",
        ont = "ALL",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff,  # set to 1 in order to obtain all possible terms 
        pAdjustMethod = "BH",
        by = "fgsea",
        # eps = eps,
        # verbose = FALSE,  # FALSE: Show ontology separately
        seed = TRUE
    )
    
    if (nrow(res) < 1) return(NULL)
    
    if (simplify) 
    {
        res <- res %>% 
            clusterProfiler::simplify() %>% 
            as.data.frame()
    }
    
    .to_gene_id <- function(entrez_id)
    {
        gene_id <- entrez_id %>% 
            str_split_1("/") %>% 
            entrez2tair()
        gene_id <- unname(gene_id[complete.cases(gene_id)])
        gene_id <- paste(gene_id, collapse = "/")
        return(gene_id)
    }
    
    suppressMessages(
        res$gene_id <- vapply(res$core_enrichment, .to_gene_id, character(1))
    )
    
    return(res)
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
GO_dotplot <- function(
        go_df, 
        title = "",
        subtitle = "",
        text_size = 24,
        text_wrap = 17,
        text_angle = 0,
        scale_x_limits = NULL,
        scale_x_breaks = seq(0, scale_x_limits[2], length.out = 5),
        save_dir = NULL,
        dpi = 660,
        height = 19,
        width = 27
){
    if (is.null(scale_x_limits))
        scale_x_limits <- as.integer(c(0, ceiling(max(go_df[["fold_enrich"]]))))
    
    if (is.null(scale_x_breaks))
        scale_x_breaks <- seq(0, scale_x_limits[2], length.out = 5)
    
    p1 <- go_df %>% 
        ggplot(
            mapping = aes(x = fold_enrich, y = fct_reorder(Description, fold_enrich))
        ) +
        theme_bw() +
        labs(
            title = title,
            subtitle = subtitle,
            x = "Fold enrichment",
            size = "Count",
            color = "<i>p<sub>adj</sub></i>"
        ) +
        geom_point(aes(size = Count, color = p.adjust)) +
        geom_segment(
            mapping = aes(xend = min(scale_x_limits), 
                          yend = Description, 
                          color = p.adjust, 
                          linewidth = Count), 
            show.legend = FALSE
        ) +
        scale_x_continuous(limits = scale_x_limits, breaks = scale_x_breaks) +
        scale_y_discrete(labels = function(x) str_wrap(x, text_wrap)) +
        scale_color_gradientn(
            colours = c("#FFB000", "#FE6100", "#DC267F", "#785EF0"),
            trans = "log10",
            guide = guide_colorbar(
                reverse = FALSE,
                order = 1,
                title.theme = element_markdown(margin = ggplot2::margin(b = 15)),
                theme = theme(legend.key.size = grid::unit(1, "cm"))
            )
        ) +
        scale_size_continuous(
            range = c(7, 15),
            guide = guide_legend(
                order = 2,
                theme = theme(legend.key.size = grid::unit(0.7, "cm"))
            )
        ) +
        scale_linewidth_continuous(
            range = c(1, 4),
            guide = guide_legend(
                order = 2,
                theme = theme(legend.key.size = grid::unit(0.7, "cm"))
            )
        ) +
        theme(
            text = element_text(family = "sans", face = "bold", size = text_size),
            plot.title = element_markdown(),
            plot.subtitle = element_markdown(),
            axis.title.x = element_markdown(margin = ggplot2::margin(t = 7)),
            axis.text.x.bottom = element_text(size = text_size - 2, angle = text_angle),
            axis.title.y = element_blank(),
            axis.text.y = element_text(family = "sans", face = "bold", size = text_size),
            legend.title = element_markdown(
                family = "sans", 
                face = "bold", 
                size = text_size - 4,
                padding = ggplot2::margin(b = 0)
            ),
            legend.text = element_text(family = "sans", face = "plain", size = text_size - 4),
            legend.background = element_blank(),
            legend.position = "right"
        )
    
    if (!is.null(save_dir)){
        ggsave(
            plot = p1, 
            filename = basename(save_dir), 
            path = dirname(save_dir),
            device = "jpeg",
            units = "cm",
            dpi = dpi,
            height = height,
            width = width 
        )
    }
    
    return(p1)
}


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

ggvenn3 <- function(
        data = list(), 
        
        show_set_names = TRUE,
        set_names = c(),
        set_names_text_size = 6,
        set_names_text_angle = c(0, 0, 0),
        
        show_labels = TRUE,
        label_text_size = 4.5,
        
        text_size = 4,
        stroke_size = NA,
        
        fill_color = c("#56B4E9", "#F0E442", "#CC79A7"),
        fill_alpha = 0.33,
        ...
){
    v1 <- data[[1]]
    v2 <- data[[2]]
    v3 <- data[[3]]
    # layer 1: polygon
    # layer 2: polygon stroke
    # layer 3: set names
    # layer 4: set data
    set_elements <- list(
        set_A = v1 %>% setdiff(v2) %>% setdiff(v3),
        set_B = v1 %>% intersect(v2) %>% setdiff(v3),
        set_C = v2 %>% setdiff(v1) %>% setdiff(v3),
        set_D = v1 %>% intersect(v3) %>% setdiff(v2),
        set_E = v1 %>% intersect(v2) %>% intersect(v3),
        set_F = v2 %>% intersect(v3) %>% setdiff(v1),
        set_G = v3 %>% setdiff(v1) %>% setdiff(v2)
    )
    
    p0 <- ggvenn::ggvenn(
        data = data,
        text_size = text_size,
        stroke_size = stroke_size,
        fill_color = fill_color,
        fill_alpha = fill_alpha,
        ...
    ) +
        theme(
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            text = element_text(size = 22, face = "bold", family = "sans")
        )
    
    p0[["layers"]][[3]] <- NULL
    
    if (show_set_names)
    {
        if (length(set_names) == 0) 
            set_names <- names(data)
        
        p0 <- p0 +
            ggplot2::annotate(
                geom = "richtext",
                x = c(-0.75, 0.75, 0),
                y = c(1.8, 1.8, -1.8),
                label = set_names,
                size = set_names_text_size,
                fontface = "bold",
                angle = set_names_text_angle,
                fill = NA,
                label.color = NA
            )
    }
    
    if (show_labels)
    {
        p0 <- p0 +
            ggplot2::annotate(
                geom = "text",
                # x = venn_label_position$x,
                x = c(
                    "A" = -0.75, "B" = 0, "C" = 0.75, 
                    "D" = -0.85, "E" = 0, "F" = 0.85, "G" = 0
                ),
                # y = venn_label_position$y,
                y = c(
                    "A" = 1.5, "B" = 1.2, "C" = 1.5, 
                    "D" = -0.3, "E" = -0.1, "F" = -0.3, "G" = -1.4
                ),
                label = LETTERS[1:7],  # A to G
                size = label_text_size,
                fontface = "bold"
            )
    }
    
    return(list(figure = p0, set_elements = set_elements))
}

ggvenn4 <- function(
        data = list(), 
        
        show_set_names = TRUE,
        set_names = c(),
        set_names_text_size = 6,
        
        show_labels = TRUE,
        label_text_size = 4.5,
        
        text_size = 4,
        stroke_size = NA,
        
        fill_color = c("#CC79A7", "#009E73", "#E69F00", "#56B4E9"),
        # fill_color = c("#882255", "#117733", "#DDCC77", "#88CCEE"),
        fill_alpha = 0.3,
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
    
    p0 <- ggvenn::ggvenn(
        data = data,
        text_size = text_size,
        stroke_size = stroke_size,
        fill_color = fill_color,
        fill_alpha = fill_alpha,
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
            ggplot2::annotate(
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
            ggplot2::annotate(
                geom = "text",
                # x = venn_label_position$x,
                x = c(
                    "A" = -0.85, "B" = 0, "C" = 0.85, 
                    "D" = -1.1, "E" = 1.1, 
                    "F" = -1.5, "G" = 1.5, 
                    "H" = -0.5, "I" = 0.5, 
                    "J" = -0.85, "K" = 0, "L" = 0.85, 
                    "M" = -0.35, "N" = 0.35, "O" = 0
                ),
                # y = venn_label_position$y,
                y = c(
                    "A" = 0.95, "B" = 0.65, "C" = 0.95, 
                    "D" = 0.55, "E" = 0.55, 
                    "F" = 0.40, "G" = 0.40, 
                    "H" = 0.10, "I" = 0.10, 
                    "J" = -0.60, "K" = -0.45, "L" = -0.60, 
                    "M" = -0.90, "N" = -0.90, "O" = -1.17
                ),
                label = LETTERS[1:15],  # A to O
                size = label_text_size,
                fontface = "bold"
            )
    }
    
    return(list(figure = p0, set_elements = set_elements))
}


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
    # library(biomartr) 
})


# NOT RUN {
if (FALSE){
    tair_id <- "AT1G08830"
    go_id <- c("GO:0070482", "GO:0071453", "GO:0071456")
}
# }

parse_eval <- function(v){
    sapply(v, function(.x) eval(parse(text = .x)))
}

# ---- Retrieve GO data based on TAIR GeneID -----------------------------------------
tair_desc <- function(tair_id){
    stopifnot(inherits(tair_id, "character"))
    
    tair_ <- suppressMessages({
        AnnotationDbi::select(
            x = org.At.tair.db,
            keys = tair_id,
            keytype = "TAIR",
            columns = c(
                "SYMBOL", 
                "GENENAME", 
                "GOALL", 
                "ARACYC",  # exp. superoxide radicals degradation
                "ARACYCENZYME",  # exp: superoxide dismutase
                "ENZYME"  # exp: 1.15.1.1
            )
        ) %>% 
            dplyr::select(-EVIDENCEALL, -ONTOLOGYALL) %>% 
            dplyr::group_by(TAIR) %>% 
            dplyr::summarise_all(function(x) paste(unique(x), collapse = "/"))
    })
    
    return(tair_)
}


# ---- Retrieve GO data based on TAIR GeneID -----------------------------------------
go_desc <- function(go_id){
    stopifnot(inherits(go_id, "character"))
    
    tair_all_GO <- AnnotationDbi::keys(org.At.tair.db, "GOALL")
    avail_GO <- go_id[go_id %in% tair_all_GO]
    
    go_desc_ <- suppressMessages({
        AnnotationDbi::select(
            x = org.At.tair.db,
            keys = avail_GO,
            keytype = "GOALL",
            columns = "TAIR"
        ) %>% 
            dplyr::select(-EVIDENCEALL, -ONTOLOGYALL) %>% 
            dplyr::group_by(GOALL) %>%
            dplyr::summarise_all(function(x) paste(unique(x), collapse = "/"))
    })
    
    go_ <- suppressMessages({
        AnnotationDbi::select(
            x = GO.db,
            keys = avail_GO,
            keytype = "GOID",
            columns = AnnotationDbi::columns(GO.db)
        ) %>% 
            dplyr::rename(GOALL = GOID) %>% 
            dplyr::left_join(go_desc_, by = "GOALL")
    })
    
    return(go_)
}


# ---- Retrieve GO data based on TAIR GeneID -----------------------------------------
tair_enrichGO <- function(
        tair_id, 
        universe_tair_id,
        ontology = "ALL",
        showCategory = 10,
        returnData = FALSE
){
    eg_ <- clusterProfiler::enrichGO(
        gene = tair_id,
        universe = universe_tair_id,
        OrgDb = org.At.tair.db,
        keyType = "TAIR",
        ont = ontology,  # c("ALL", "BP", "CC", "MF")
        qvalueCutoff = 0.2,  # default: 0.2
        pvalueCutoff = 0.05,  # default: 0.05
        pAdjustMethod = "BH",  # default: "ALL"
        readable = FALSE,  # TRUE: show Gene Symbol, FALSE: show TAIR_id
        pool = FALSE  # FALSE: Show ontology separately
    ) %>% 
        clusterProfiler::simplify()
    
    eg_df_ <- as.data.frame(eg_@result) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate_at(vars(GeneRatio, BgRatio), parse_eval) %>% 
        tidyr::drop_na() %>% 
        dplyr::mutate(fold.enrich = as.numeric(GeneRatio) / as.numeric(BgRatio)) %>% 
        arrange(desc(GeneRatio))
    
    eg_df_$Description <- factor(
        x = eg_df_$Description, 
        levels = eg_df_ %>% dplyr::arrange(GeneRatio) %>% .$Description
    )
    
    if (!"ONTOLOGY" %in% colnames(eg_df_)) eg_df_$ONTOLOGY <- rep(ontology, nrow(eg_df_))
    
    if (nrow(eg_df_) > showCategory) eg_df_ <- eg_df_[1:showCategory, ]
    # eg@result$fold.enrich <- eg_df$fold.enrich
    # 
    # eg_dotplot <- clusterProfiler::dotplot(
    #     object = eg,
    #     title = "DEGs",
    #     showCategory = 10,
    #     x = "GeneRatio",
    #     size = "fold.enrich",
    #     orderBy = "x"
    # )
    
    eg_dotplot_ <- ggplot(
        data = eg_df_, 
        mapping = aes(Description, GeneRatio, color = p.adjust, size = fold.enrich, shape = ONTOLOGY)
    ) +
        geom_point() +
        theme_bw() +
        scale_x_discrete(labels = scales::label_wrap(25)) +
        scale_color_continuous(type = "gradient") +
        coord_flip()
    
    ifelse(returnData, return(eg_df_), return(eg_dotplot_))
}




# ---- Show description of the GO EVIDENCE CODE -----------------------------------------
showGOEvidenceCodes <- function(detailed = FALSE){
    if (!detailed){
        print(c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                "ISS", "ISO", "ISA", "ISM", "IGC", "IBA",
                "IBD", "IKR", "IRD", "RCA", "TAS", "NAS", 
                "IC",  "ND",  "IEA", "NR" )
        )
    }else{
        cat("
    GO Term Evidence Code
    
    #> Experimental Evidence Codes
    EXP: Inferred from Experiment
    IDA: Inferred from Direct Assay
    IPI: Inferred from Physical Interaction
    IMP: Inferred from Mutant Phenotype
    IGI: Inferred from Genetic Interaction
    IEP: Inferred from Expression Pattern
    
    #> Computational Analysis Evidence Codes
    ISS: Inferred from Sequence or Structural Similarity
    ISO: Inferred from Sequence Orthology
    ISA: Inferred from Sequence Alignment
    ISM: Inferred from Sequence Model
    IGC: Inferred from Genomic Context
    IBA: Inferred from Biological aspect of Ancestor
    IBD: Inferred from Biological aspect of Descendant
    IKR: Inferred from Key Residues
    IRD: Inferred from Rapid Divergence
    RCA: inferred from Reviewed Computational Analysis
    
    #> Author Statement Evidence Codes
    TAS: Traceable Author Statement
    NAS: Non-traceable Author Statement
    #> Curator Statement Evidence Codes
    IC: Inferred by Curator
    ND: No biological Data available
    Automatically-assigned Evidence Codes
    IEA: Inferred from Electronic Annotation
    
    #> Obsolete Evidence Codes
    NR: Not Recorded
          ")
    }
}



suppressMessages({
    if (!require(dplyr)) BiocManager::install("dplyr")
    if (!require(magrittr)) BiocManager::install("magrittr")
    if (!require(ggplot2)) BiocManager::install("ggplot2")
})


plotFPKM <- function(dds, gene_id, intgroup, returnData = FALSE){
    gene_id <- sym(gene_id)
    gene_id_str <- deparse(substitute(gene_id))
    
    fpkm_mat <- DESeq2::fpkm(dds)
    fpkm_mat <- fpkm_mat[rownames(fpkm_mat) == gene_id_str, , drop = FALSE] %>% 
        t() %>% 
        as.data.frame() %>% 
        cbind(., dds@colData) %>% 
        dplyr::select(all_of(c(gene_id_str, intgroup)))
    
    fpkm_mat$intgroup <- fpkm_mat %>% 
        dplyr::select(intgroup) %>% 
        apply(., 1, function(x) paste(x, collapse = "_"))
    
    if (returnData){
        return(fpkm_mat)
    }else{
        gp <- ggplot(fpkm_mat, aes(intgroup, {{gene_id}})) +
            labs(title = gene_id_str, x = element_blank(), y = "FPKM") +
            geom_point(position = position_jitter(width = 0.1), alpha = .7) +
            geom_boxplot(outlier.shape = NA, fill = "transparent") +
            theme_bw() +
            theme(
                axis.title.x.bottom = element_blank()
            )
        return(gp)
    }
}


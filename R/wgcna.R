# Weighted Gene Correlation Network Analysis

if (!require(ps)) install.packages("ps", quiet = TRUE)
if (!require(tidyverse)) install.packages("tidyverse", quiet = TRUE)
if (!require(WGCNA)) BiocManager::install("WGCNA", ask = FALSE, update = FALSE)
if (!require(flashClust)) BiocManager::install("flashClust", ask = FALSE, update = FALSE)

.outliers_filtering <- function(expr_mat){
    res <- list()
    gsg <- WGCNA::goodSamplesGenes(expr_mat)
    
    if (any(!gsg$goodGenes)) {
        res[["outlier_genes"]] <- colnames(expr_mat)[!gsg$goodGenes]
        expr_mat <- expr_mat[gsg$goodGenes]
    }
    
    if (any(!gsg$goodSamples)) {
        res[["outlier_samples"]] <- rownames(expr_mat)[!gsg$goodSamples]
        expr_mat <- expr_mat[gsg$goodSamples, ]
    }
    
    return(res)
}

.show_sample_clustering <- function(expr_mat, dist_method, hclust_method){
    sample_tree <- hclust(
        d = dist(expr_mat, method = dist_method),
        method = hclust_method
    )
    par(cex = 0.6);
    par(mar = c(0,4,2,0));
    plot(
        x = sample_tree,
        main = sprintf("Sample clusters (%s + %s)", dist_method, hclust_method),
        sub = "",
        xlab = "",
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 2
    )
    return(sample_tree)
}

wgcna <- function(
        expr_mat, 
        filter_outliers = TRUE,
        show_sample_clustering = FALSE,
        dist_method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
        hclust_method = c("complete", "average", "ward", "single", "mcquitty", "median", "centroid"),
        power_n = NULL,
        allow_multithreads = TRUE,
        ...
){
    cor <- WGCNA::cor
    res <- list()
    dist_method <- match.arg(dist_method)
    hclust_method <- match.arg(hclust_method)
    
    # Multi-threads ====
    use_cores <- 1
    if (allow_multithreads){
        use_cores <- ps::ps_cpu_count() - 1
        WGCNA::allowWGCNAThreads(nThreads = use_cores)
    }
    
    # expr_mat: (rownames: gene_id, colnames: sample_id) ====
    expr_mat[is.na(expr_mat)] <- 0
    expr_mat <- t(expr_mat[rowSums(expr_mat) > 0, ])
    # df0 <- as.data.frame(expr_mat)
    res[["data"]] <- expr_mat
    
    # Outliers filtering ====
    if (filter_outliers) {
        res <- c(res, .outliers_filtering(expr_mat))
    }
    
    # Samples clustering ====
    if (show_sample_clustering) {
        res[["sample_tree"]] <- .show_sample_clustering(expr_mat, dist_method, hclust_method)
    }
    
    
    # Soft-picked threshold ====
    # spt <- pickSoftThreshold(df0)
    spt <- pickSoftThreshold(expr_mat)
    res[["spt_fit"]] <- spt$fitIndices
    if (is.null(power_n)) power_n <- spt$powerEstimate
    if (is.na(power_n)) power_n <- which.max(spt$fitIndices$SFT.R.sq)
    tmp_plot <- spt$fitIndices %>% 
        select(Power, SFT.R.sq, mean.k.) %>% 
        pivot_longer(cols = c(SFT.R.sq, mean.k.), names_to = "parameters", values_to = "val") %>% 
        ggplot(aes(Power, val)) +
        theme_bw() +
        facet_wrap( ~ parameters, scales = "free_y") +
        geom_text(aes(label = Power)) +
        geom_vline(xintercept = power_n, linetype = "dashed", color = "maroon")
    show(tmp_plot)
    res[["power_estimate"]] <- power_n
    
    network <- WGCNA::blockwiseModules(
        datExpr = expr_mat,
        power = power_n,
        # networkType = "signed",
        # deepSplit = 2,
        # pamRespectsDendro = FALSE,
        # detectCutHeight = 0.75,
        # minModuleSize = 30,
        # maxBlockSize = 4000,
        # reassignThreshold = 0,
        # mergeCutHeight = 0.25,
        # numericLabels = TRUE,
        nThreads = use_cores,
        verbose = 3,
        ...
    )
    res[["network"]] <- network
    
    # Modules genes ====
    module_df <- data.frame(
        gene_id = names(network$colors),
        modules = network$colors
    )
    res[["modules"]] <- module_df
    
    # Eigen-genes of each modules ====
    eigen_genes <- moduleEigengenes(expr_mat, network$colors)$eigengenes
    eigen_genes <- orderMEs(eigen_genes)
    colnames(eigen_genes) <- gsub("ME", "", colnames(eigen_genes))
    res[["eigen_genes"]] <- eigen_genes
    
    # # Adjacency matrix of the weighted network ====
    # adjacency <- adjacency(expr_mat, power = power_n)
    # 
    # # Topological Overlap Matrix (TOM) ====
    # TOM_similarity <- TOMsimilarity(adjacency)
    # TOM_dissimilarity <- 1 - TOM_similarity
    # 
    # # Hierarchical Clustering Analysis ====
    # gene_tree <- hclust(as.dist(TOM_dissimilarity), method = hclust_method)
    # gene_modules <- cutreeDynamic(
    #     dendro = gene_tree,
    #     distM = TOM_dissimilarity,
    #     deepSplit = 2,
    #     pamRespectsDendro = FALSE,
    #     minClusterSize = 30
    # )
    # # gene_module_colors <- labels2colors(gene_modules)
    # 
    # res[["gene_modules"]] <- gene_modules
    
    # Eigengenes identification from the modules ====
    # module_eigen_genes <- moduleEigengenes(expr_mat, colors = labels2colors(gene_modules))
    # var_exp <- list_c(module_eigen_genes$varExplained[, , drop = TRUE])
    # PVE <- cumsum(var_exp / sum(var_exp))
    # res[["module_eigen_genes"]] <- module_eigen_genes
    # res[["PVE"]] <- PVE
    cor <- stats::cor
    return(res)
}












count2fpkm <- function(count_mat, gene_length, returnTPM = FALSE){
    stopifnot(is.vector(gene_length))
    stopifnot(inherits(names(gene_length), "character"))
    stopifnot(identical(rownames(count_mat), names(gene_length)))
    
    norm_cmat <- t(t(count_mat) / (colSums(count_mat) / 1e6))
    fpkm_mat <- norm_cmat / (gene_length / 1000)
    
    if (returnTPM){
        fpkm_mat <- t(t(fpkm_mat) / colSums(fpkm_mat)) * 1e6
    }
    
    return(fpkm_mat)
}


count2tpm <- function(count_mat, gene_length){
    stopifnot(is.vector(gene_length))
    stopifnot(inherits(names(gene_length), "character"))
    stopifnot(identical(rownames(count_mat), names(gene_length)))
    
    norm_mat <- count_mat / (gene_length / 1000)
    tpm_mat <- t(t(norm_mat) / (colSums(norm_mat) / 1e6))
    
    return(tpm_mat)
}


fpkm2tpm <- function(fpkm_mat){
    tpm_mat <- t(t(fpkm_mat) / colSums(fpkm_mat)) * 1e6
    return(tpm_mat)
}


tpm2fpkm <- function(tpm_mat){
    next
}
    
    
    
    
    
    
    
suppressMessages({
    if (!require(readr)) install.packages("readr")
    library(readr)
    if (!require(stringr)) install.packages("stringr")
    library(stringr)
})


read_gtf <- function(gtf_file){
    gtf <- read_tsv(gtf_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
    colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
    return(gtf)
}


write_gtf <- function(x, file){
    write_tsv(x, file, col_names = FALSE, quote = "none", escape = "none")
    
}


extract_attribute_feature <- function(x, feature_name = NULL){
    if (is.null(feature_name)){
        cat("\n'feature_name' should not be null\n")
    }
    
    if (is.vector(x)){
        feat <- gsub(paste0('.*', feature_name, ' "([^"]+)".*'), "\\1", x)
        feat[sapply(feat, function(x) grepl(";", x))] <- NA_character_
    }
    
    if (is.data.frame(x)){
        feat <- x[ , tail(colnames(x), 1), drop = TRUE]
        feat <- gsub(paste0('.*', feature_name, '\\s"([^"]+)".*'), "\\1", feat)
        feat[sapply(feat, function(x) grepl(";", x))] <- NA_character_
    }
    
    return(feat)
}


gff2gtf <- function(gff_file){
    
    col_names <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
    
    if (inherits(gff_file, "character")){
        gff_file <- read_tsv(gff_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
    }
    
    if (identical(colnames(gff_file), col_names)){
        colnames(gff_file) <- col_names
    }
    
    gff <- gff_file %>% 
        mutate(attributes = str_replace_all(attributes, "=", " \"")) %>% 
        mutate(attributes = str_replace_all(attributes, ";", "\"; "))
    
    return(gff)
}





suppressMessages({
    if (!require(readr)) install.packages("readr")
    if (!require(stringr)) install.packages("stringr")
    if (!require(dplyr)) install.packages("dplyr")
})


gtf_colnames <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

read_gtf <- function(gtf_file){
    gtf <- readr::read_tsv(gtf_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
    colnames(gtf) <- gtf_colnames
    gtf$start <- as.numeric(format(gtf$start, scientific = FALSE))
    gtf$end <- as.numeric(format(gtf$end, scientific = FALSE))
    gtf$attributes <- as.character(gtf$attributes)
    return(gtf)
}


write_gtf <- function(x, file){
    readr::write_tsv(x, file, col_names = FALSE, quote = "none", escape = "none")
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
    
    if (inherits(gff_file, "character")){
        gff_file <- readr::read_tsv(gff_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
    }
    
    if (identical(colnames(gff_file), gtf_colnames)){
        colnames(gff_file) <- gtf_colnames
    }
    
    gff <- gff_file %>% 
        dplyr::mutate(attributes = stringr::str_replace_all(attributes, "=", " \"")) %>% 
        dplyr::mutate(attributes = stringr::str_replace_all(attributes, ";", "\"; ")) %>% 
        dplyr::mutate(attributes = paste0(attributes, "\";"))
    
    return(gff)
}





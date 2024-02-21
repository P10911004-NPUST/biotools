# options(dplyr.summarise.inform = FALSE)

# # Load self-defined functions
# sf <- new.env()
# sapply(list.files("./utils", pattern = ".R", full.names = TRUE), 
#        function(x){source(x, local = sf, chdir = TRUE)})

.library_chr_only <- function(x) library(x, character.only = TRUE)

load_packages <- function(...){
    pkgs_list <- sapply(match.call(), as.character)[-1]
    installed_pkgs <- rownames(installed.packages())
    direct_load_pkgs <- intersect(pkgs_list, installed_pkgs)
    
    sapply(direct_load_pkgs, .library_chr_only)
    pkgs_list <- setdiff(pkgs_list, direct_load_pkgs)
    
    if (length(pkgs_list) != 0){
        if (!require(utils)) install.packages("utils", quiet = TRUE)
        if (!require(BiocManager, quietly = TRUE)) install.packages("BiocManager", quiet = TRUE)
        
        cran_pkgs <- rownames(utils::available.packages(quiet = TRUE))
        bioc_pkgs <- setdiff(BiocManager::available(), cran_pkgs)
        
        download_from_cran <- pkgs_list[pkgs_list %in% cran_pkgs]
        download_from_bioc <- pkgs_list[pkgs_list %in% bioc_pkgs]
        
        download_from_other <- setdiff(pkgs_list, union(download_from_cran, download_from_bioc))
        if (length(download_from_other) == 0) download_from_other <- ""
        
        sapply(download_from_cran, function(x) install.packages(x, dependencies = TRUE, quiet = TRUE))
        sapply(download_from_bioc, function(x) BiocManager::install(x, update = FALSE, ask = FALSE))
        
        sapply(setdiff(pkgs_list, download_from_other), .library_chr_only) 
        
        if (length(download_from_other) == 0) download_from_other <- 0
        message(sprintf("Packages not found in CRAN and Bioconductor: %s", paste(download_from_other, collapse = ", ")))
    }
}


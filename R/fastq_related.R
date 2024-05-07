suppressMessages({
    source("./load_packages.R")
    source("./os_related.R")
    load_packages(Biostrings, GenomicAlignments, Rsamtools, ShortRead, Rfastp)
})

extract_barcodes <- function(filepath, format = "fastq"){
    tryCatch(
        expr = {
            if (format == "fastq"){
                input_seq <- readFastq(filepath)
            }else{
                input_seq <- readFasta(filepath)
            }
            ids <- as.character(input_seq@id)
            
            rm(input_seq)
            gc()
            
            ids <- sapply(strsplit(ids, ":"), function(x) tail(x, 1))
            ids <- unique(ids)
        }, 
        error = function(err){
            print(err)
        },
        finally = {
            return(ids)
        }
    )
}


.seqkit_extract_barcode <- function(inFile, outFile){
    # Install seqkit: https://bioinf.shenwei.me/seqkit/download/
    .command <- paste(
        # seqkit command --> extract header from fastq
        "seqkit fx2tab", "--threads", use_cores, "--name", inFile,
        "|",
        # separate the header by ':' and extract the last column
        "awk -F ':' '{print $NF}'",
        "|",
        # extract unique values
        "cut -d\\n -f1 | sort | uniq",
        ">",
        outFile,
        sep = " "
    )
    cat(.command)
    system(.command)
}


readLines_iteratively <- function(filepath, n = 1){
    con <- file(filepath, "r")
    while ( TRUE ) {
        line <- readLines(con, n = n)
        if ( length(line) == 0 ) break
        # if (TRUE){
        #     line <- line[grepl("@", line)]
        #     print(line)
        # }
    }
    close(con)
}

get_seq_from_granges <- function(genome, gr_obj){
    for (i in 1:length(gr_obj)){
        gr <- gr_obj[i]
        if (i == 1){
            res <- genome[gr]
        }else{
            res <- append(res, genome[gr])
        }
    }
    names(res) <- names(gr_obj)
    return(res)
}


fasta_or_fastq <- function(file_dir){
    fa_pattern <- paste(paste0(c("\\.fa", "\\.fasta"), "(\\.gz)?$"), collapse = "|")
    fq_pattern <- paste(paste0(c("\\.fq", "\\.fastq"), "(\\.gz)?$"), collapse = "|")
    if (missingArg(file_dir)) stop("Please input the file directory")
    stopifnot(grepl(fa_pattern, file_dir) | grepl(fq_pattern, file_dir))
    if (grepl(fa_pattern, file_dir)) return("fasta")
    if (grepl(fq_pattern, file_dir)) return("fastq")
}


# intergenic <- readDNAStringSet("D:/bcst/YAMADA/jklai/Araport11/Araport11_intergenic_20220504.fa")
# w <- width(intergenic)
# avg <- ceiling(mean(w))
# std <- ceiling(sd(w))
# max_ <- max(w)
# min_ <- min(w)
# cat("Intergenic region : ", avg, "\U00B1", std)


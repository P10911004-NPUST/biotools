rm(list = ls()); gc()

## NOT RUN {
### Execute in terminal
if (FALSE){"
scp ./Rsubread_miRNA.R holst:/bcst/YAMADA/jklai/rcode/
ssh holst
module load slurm R
cd /bcst/YAMADA/jklai
srun --exclusive --mincpus=40 --mem 400000 Rscript ./rcode/Rsubread_miRNA.R
"}
## }

if (!require(BiocManager)) install.packages("BiocManager")
if (!require(benchmarkme)) install.packages("benchmarkme")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(Rsubread)) BiocManager::install("Rsubread")
library(benchmarkme)
library(tidyverse)
library(Rsubread)

use_cores <- floor(get_cpu()$no_of_cores * 0.99)
# use_ram <- get_ram() * 0.99
sprintf("Use / Available CPUs: %s / %s", use_cores, get_cpu()$no_of_cores)
# sprintf("Use / Available RAM: %.2f GB / %.2f GB", use_ram / 1e9, get_ram() / 1e9)


# ---- Project directory --------------------------------------------------------------------------
project_dir <- "D:/bcst/YAMADA/jklai/project/Project_nx079_20230330_S1_Masashi_Yamada"


# ----- Reference genome (FASTA file) ------------------------------------------------
ref_genome <- list.files("D:/bcst/YAMADA/jklai/Araport11", pattern = "\\.fa$", full.names = TRUE)


# ----- Samples sequences (folders containing FASTQ files) ----------------------------------
sample_list_dir <- "D:/bcst/YAMADA/jklai/rawdata/Project_nx079_20230330_S1_Masashi_Yamada"
sample_list <- list.files(sample_list_dir, full.names = TRUE)


# ---- Annotation (GTF file) ------------------------------------------------------------------
anno_gtf <- list.files("/bcst/YAMADA/jklai/Araport11", pattern = "\\.gtf$", full.names = TRUE)

####################################################################################################
# ---- Index FASTA --------------------------------------------------------------------------
idx_folder <- file.path(dirname(ref_genome), "index")
if (!dir.exists(idx_folder)) dir.create(idx_folder, recursive = TRUE)

if (length(list.files(idx_folder, pattern = "ref_idx")) == 0){
    print(paste0("Building index for reference genome in ", idx_folder))
    Rsubread::buildindex(
        basename = file.path(idx_folder, "ref_idx"),
        reference = ref_genome,
        # memory = use_ram,
        indexSplit = TRUE
    )
}else{
    print("Indexed reference genome already exists...")
}


# ---- BAM folder ------------------------------------------------
BAM_folder <- file.path(project_dir, "BAM")
if (!dir.exists(BAM_folder)) dir.create(BAM_folder, recursive = TRUE)


# ---- Mapping miRNA ------------------------------------------------------------------
for (i in sample_list){
    FASTQ_file <- list.files(i, pattern = "\\.fastq\\.gz", full.names = TRUE)
    
    stopifnot(length(FASTQ_file) == 2)
    
    BAM_file <- paste0(file.path(BAM_folder, basename(i)), ".bam")
    
    if (!file.exists(BAM_file)){
        align(
            index = file.path(idx_folder, "ref_idx"),
            readfile1 = FASTQ_file[1],
            readfile2 = FASTQ_file[2],
            type = 0,  # RNA-seq data
            PE_orientation = "fr",  # Forward-Reverse
            output_file = BAM_file,
            useAnnotation = TRUE,
            annot.ext = anno_gtf,
            isGTF = TRUE,
            output_format = "BAM",
            sortReadsByCoordinates = TRUE,    # sorting and indexing (.bai)
            TH1 = 7,  # only >= 22 bp can be detected (at least 22 perfectly matched bases)
            nsubreads = 35,
            maxMismatches = 3,
            indels = 0,
            unique = FALSE,
            nBestLocations = 10,
            nthreads = use_cores
        )
    }else{
        print(paste0("BAM file already exists: ", basename(BAM_file)))
    }
}


# ---- Count miRNA ---------------------------------------------------------------------
BAM_list <- list.files(BAM_folder, pattern = "\\.bam$", full.names = TRUE)

FC_folder <- file.path(project_dir, "FC")
if (!dir.exists(FC_folder)) dir.create(FC_folder, recursive = TRUE)

fc <- featureCounts(
    # BAM files
    files = BAM_list,
    
    # Annotation GTF file
    annot.ext = anno_gtf,
    isGTFAnnotationFile = TRUE,
    
    # Set feature type as miRNA
    GTF.featureType = "miRNA_primary_transcript",
    
    # Paired-end reads
    countReadPairs = TRUE,
    isPairedEnd = TRUE,
    
    # Strand-specific (0: unstranded, 1: stranded, 2: reversely stranded)
    strandSpecific = 1,
    
    # Number of threads
    nthreads = use_cores,
    tmpDir = FC_folder
)

fc_df <- as.data.frame(fc$counts) %>%
    tibble::rownames_to_column(var = "GeneID") %>%
    dplyr::left_join(fc$annotation, ., by = "GeneID")

write.table(
    x = fc_df, 
    file = file.path(FC_folder, "count_matrix_miRNA.txt"), 
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)



# ---- Start counting ---------------------------------------------------------------
# paste(
#     "featureCounts",
#     
#     # Paired-end reads
#     paste0("-p --countReadPairs"),
#     
#     # Annotation GTF file
#     paste0("-a ", anno_gtf),
#     
#     # Set feature type as exon
#     paste0("-t exon"),
#     
#     # Specify the attribute type used to group features
#     paste0("-g gene_id"),
#     
#     # Strand-specific (0: unstranded, 1: stranded, 2: reversely stranded)
#     paste0("-s 1"),
#     
#     # Number of threads
#     paste0("-T ", use_cores),
#     
#     # Output feature counts table
#     # paste0("-o ", file.path(FC_folder, basename(i)), ".txt"),
#     paste0("-o ", file.path(FC_folder, "count_matrix.txt")),
#     
#     # Input the sorted-indexed BAM file
#     paste(BAM_list, collapse = " "),
#     
#     sep = " "
# ) |>
#     system()


rm(list = ls()); gc()

# NOT RUN
# Execute in terminal
if (FALSE){"
scp ./hisat2.R holst:/bcst/YAMADA/jklai/rcode/
ssh holst
module load slurm R hisat2 samtools
cd /bcst/YAMADA/jklai
srun --exclusive --mincpus=40 --mem 400000 Rscript rcode/hisat2.R
"}

if (!require(stringr)) install.packages("stringr")
if (!require(benchmarkme)) install.packages("benchmarkme")
library(stringr)
library(benchmarkme)

use_cores <- floor(get_cpu()$no_of_cores * 0.99)
use_ram <- get_ram() * 0.99
sprintf("Use / Available CPUs: %s / %s", use_cores, get_cpu()$no_of_cores)
sprintf("Use / Available RAM: %.2f GB / %.2f GB", use_ram / 1e9, get_ram() / 1e9)

# ----- Reference genome (FASTA file) ------------------------------------------------
ref_genome <- list.files("/bcst/YAMADA/jklai/Araport11", pattern = "\\.fa$", full.names = TRUE)

# ----- Dataset ----------------------------------------------------------------------
sample_list_dir <- "/bcst/YAMADA/jklai/rawdata/Project_nx053_20220922_S1_Masashi_Yamada"
sample_list <- list.files(sample_list_dir, full.names = TRUE)

# ----- Output BAM directory ---------------------------------------------------------
project_dir <- gsub("rawdata", "project", sample_list_dir)
BAM_folder <- file.path(project_dir, "BAM")
if (!dir.exists(BAM_folder)) dir.create(BAM_folder, recursive = TRUE)

idx_folder <- file.path(dirname(ref_genome), "index")


# ----- Build index for the reference genome ---------------------------------------------
if (!dir.exists(idx_folder)) dir.create(idx_folder, recursive = TRUE)
if (length(list.files(idx_folder, pattern = "ref_idx")) == 0){
    print(paste0("Building index for reference genome in ", idx_folder))
    system(
        paste(
            "hisat2-build",
            ref_genome,
            file.path(idx_folder, "ref_idx"),
            sep = " "
        )
    )
}else{
    print("Indexed reference genome already exists...")
}

# Mapping ------------------------------------------------------------------
for (i in sample_list){
    # Use fastq files after QC or use the original fastq?
    fastq_folder <- list.files(i, pattern = "\\.fastq", full.names = TRUE) %>% 
        str_subset(pattern = "trimmed", negate = TRUE)
    
    # Check whether there are 2 fastq files
    if (length(fastq_folder) != 2){
        print(paste0("Skip: ", i))
        next
    }
    
    BAM_file <- paste0(file.path(BAM_folder, basename(i)), ".bam")
    
    if (file.exists(BAM_file)){
        # Skip mapping if the SAM file already exists
        print(paste0(BAM_file, " already exists..."))
    }else{
        print(paste0("Mapping: ", i))
        paste(
            "hisat2",
            
            # Input the path and prefix of the indexed fasta
            paste0("-x ", file.path(idx_folder, "ref_idx")),
            
            # Input the first fastq file (in case of paired-end)
            paste0("-1 ", fastq_folder[1]),
            
            # Input the second fastq file (in case of paired-end)
            paste0("-2 ", fastq_folder[2]),
            
            # Input the unpaired fastq file (in case of single-end)
            ## paste0("-U ", fastq01),
            
            # lower down score threshold for reads validation
            paste0("--score-min L,0.0,-1"),
            
            # Forward-Reverse reads
            paste0("--rna-strandness FR"),
            
            # Number of threads
            paste0("--threads ", use_cores),
            
            # Generate statistic reports
            paste0("2> ", gsub("\\.bam", "\\.log", BAM_file)),
            
            # Output a SAM file and pipe to the next conversion command
            "|",
            
            # Convert SAM to BAM
            paste0("samtools view -b"),
            
            # Number of threads
            paste0("--threads ", use_cores),
            
            # Input SAM from the previous pipe operator
            "-",
            
            # output a BAM file and pipe to the next sorting command
            "|",
            
            # Sort the BAM
            paste0("samtools sort"),
            
            # Number of threads
            paste0("-@ ", use_cores),
            
            # Input BAM file from the previous pipe operator
            "-",
            
            # Output sorted-BAM filename
            paste0("-o ", BAM_file),
            
            # Execute on background
            # "&",
            
            sep = " "
        ) |>
            system()
    }
}


# Indexing the sorted-BAM ------------------------------------------------------
BAM_list <- list.files(BAM_folder, pattern = "\\.bam$", full.names = TRUE)
for (i in BAM_list){
    print(paste0("Indexing: ", i))
    
    paste(
        "samtools index",
        
        # Number of threads
        paste0("-@ ", use_cores),
        
        # Input sorted-BAM filename
        i,
        
        # Execute on background
        # "&",
        
        sep = " "
    ) |>
        system()
}

print("Finished!")


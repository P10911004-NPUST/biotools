rm(list = ls()); gc()

# NOT RUN
# Execute in terminal
if (FALSE){"
scp ./subread.R holst:/bcst/YAMADA/jklai/rcode/
ssh holst
module load slurm R subread
cd /bcst/YAMADA/jklai
srun --nodes=1 --mincpus=40 --mem=400G --job-name=subread Rscript ./rcode/subread.R
"}

if (!require(benchmarkme)) install.packages("benchmarkme")
library(benchmarkme)

use_cores <- floor(get_cpu()$no_of_cores * 0.99)
use_ram <- get_ram() * 0.99
sprintf("Use / Available CPUs: %s / %s", use_cores, get_cpu()$no_of_cores)
sprintf("Use / Available RAM: %.2f GB / %.2f GB", use_ram / 1e9, get_ram() / 1e9)

# Gene annotation (GTF file) ---------------------------------------------------
anno_gtf <- list.files("/bcst/YAMADA/jklai/Araport11", pattern = "\\.gtf$", full.names = TRUE)

# sorted- and indexed-BAM files ------------------------------------------------
project_dir <- "/bcst/YAMADA/jklai/project/Project_nx053_20220922_S1_Masashi_Yamada"
BAM_folder <- file.path(project_dir, "BAM")
BAM_list <- list.files(BAM_folder, pattern = "\\.bam$", full.names = TRUE)

# Output directory of feature counts results -----------------------------------
FC_folder <- file.path(project_dir, "FC")
if (!dir.exists(FC_folder)) dir.create(FC_folder, recursive = TRUE)


# ---------- Start counting ---------------------------------------------------------------

paste(
    "featureCounts",
    
    # Paired-end reads
    paste0("-p --countReadPairs"),
    
    # Annotation GTF file
    paste0("-a ", anno_gtf),
    
    # Set feature type as exon
    paste0("-t exon"),
    
    # Specify the attribute type used to group features
    paste0("-g gene_id"),
    
    # Strand-specific (0: unstranded, 1: stranded, 2: reversely stranded)
    paste0("-s 1"),
    
    # Number of threads
    paste0("-T ", use_cores),
    
    # Output feature counts table
    # paste0("-o ", file.path(FC_folder, basename(i)), ".txt"),
    paste0("-o ", file.path(FC_folder, "count_matrix.txt")),
    
    # Input the sorted-indexed BAM file
    paste(BAM_list, collapse = " "),
    
    sep = " "
) |>
    system()


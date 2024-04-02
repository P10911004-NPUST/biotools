suppressMessages({
  if (!require(BiocManager)) install.packages("BiocManager", quiet = TRUE)
  if (!require(Biostrings)) BiocManager::install("Biostrings", ask = FALSE, update = FALSE)
})

percent_identity <- function(x, ref){
  pa <- Biostrings::pairwiseAlignment(
    pattern = x,
    subject = ref
  )
  perc_ident <- Biostrings::pid(pa)
  return(perc_ident)
}
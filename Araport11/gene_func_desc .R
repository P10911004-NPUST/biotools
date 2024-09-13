# source("https://github.com/P10911004-NPUST/biotools/blob/main/R/isoquant_related.R?raw=TRUE")
# source("https://github.com/P10911004-NPUST/biotools/blob/main/R/gtf_related.R?raw=TRUE")
source("../R/gtf_related.R", chdir = TRUE)

features <- c(
    "mRNA",
    "exon",
    "five_prime_UTR",
    "three_prime_UTR",
    "miRNA_primary_transcript",
    "uORF",
    "pseudogenic_exon"
)

gtf <- read_gtf("./Araport11_GTF_genes_transposons.20240331.gtf") %>% 
    filter(feature %in% features) %>% 
    mutate(
        transcript_id = extract_attribute_feature(attributes, "transcript_id"),
        gene_id = extract_attribute_feature(attributes, "gene_id")
    ) %>% 
    select(gene_id, transcript_id, seqname, start, end)

gene_alias <- read_tsv(
    file = "./gene_aliases_20230407.txt",
    show_col_types = FALSE
) %>% 
    rename(gene_id = locus_name) %>% 
    summarise(
        across(c(symbol, full_name), function(x) paste(x, collapse = "; ")), 
        .by = gene_id
    )

func_desc <- read_tsv(
    file = "./Araport11_functional_descriptions_20230407.txt", 
    show_col_types = FALSE
) %>% 
    rename(transcript_id = name) %>% 
    mutate(
        Curator_summary = case_when(
            is.null(Curator_summary) | Curator_summary == "NULL" ~ Computational_description,
            TRUE ~ Curator_summary
        )
    ) %>% 
    # mutate_if(is.character, function(x) str_replace(x, "<|>|<.*>", "")) %>%
    select(transcript_id, short_description, Curator_summary, Computational_description) %>%
    rename(
        curator_summary = Curator_summary,
        computational_description = Computational_description
    )


df0 <- gtf %>% 
    left_join(gene_alias, by = "gene_id") %>% 
    left_join(func_desc, by = "transcript_id") %>% 
    distinct(gene_id, .keep_all = TRUE) %>% 
    select(-transcript_id) %>% 
    mutate(across(everything(), function(x) str_replace_all(x, "<.*>", ""))) %>% 
    mutate(across(everything(), function(x) str_replace_all(x, "&#945;", "alpha"))) %>% 
    mutate(across(everything(), function(x) str_replace_all(x, "&#946;", "beta"))) %>% 
    mutate(across(everything(), function(x) str_replace_all(x, "&#947;", "gamma"))) %>% 
    mutate(across(everything(), function(x) str_replace_all(x, "&#948;", "delta"))) %>% 
    mutate(across(everything(), function(x) str_replace_all(x, "&#952;", "theta")))

write.csv(df0, "gene_func_desc.csv", row.names = FALSE)



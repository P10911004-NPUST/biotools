if (!require(readr)) install.packages("readr")

read_assign <- function(x){
    df0 <- readr::read_tsv(
        file = x, 
        comment = "#", 
        col_names = FALSE,
        show_col_types = FALSE
    )
    colnames(df0) <- c(
        "read_id", "chr", "strand", "isoform_id", "gene_id", 
        "assignment_type", "assignment_events", "exons", "additional"
    )
    return(df0)
}

read_sqanti <- function(x){
    df0 <- readr::read_tsv(
        file = x,
        comment = "#",
        col_names = TRUE,
        show_col_types = FALSE
    )
    colnames(df0) <- c(
        "isoform", "chrom", "strand", "length", "exons",              
        "structural_category", "associated_gene", "associated_transcript", "ref_length", "ref_exons",
         "diff_to_TSS", "diff_to_TTS", "diff_to_gene_TSS", "diff_to_gene_TTS", "subcategory",
         "RTS_stage", "all_canonical", "min_sample_cov", "min_cov", "min_cov_pos",
         "sd_cov", "FL", "n_indels", "n_indels_junc", "bite",       
         "iso_exp", "gene_exp", "ratio_exp", "FSM_class", "coding",       
         "ORF_length", "CDS_length", "CDS_start", "CDS_end", "CDS_genomic_start",
         "CDS_genomic_end", "predicted_NMD", "perc_A_downstream_TTS", "seq_A_downstream_TTS", "dist_to_CAGE_peak",
         "within_CAGE_peak", "dist_to_polyA_site", "within_polyA_site", "polyA_motif", "polyA_dist",
         "polyA_motif_found", "ORF_seq", "ratio_TSS" 
    )
    return(df0)
}

assignment_events <- function(
        events = c("all", "consistent", "artifact", "intron", "inconsistent", "alt.tss.tes")
){
    event <- match.arg(events)
    
    list0 <- list(
        consistent = c(
            "none", 
            "\\.",
            "undefined",
            "mono_exon_match",
            "fsm",
            "ism_5",
            "ism_3",
            "ism_internal",
            "mono_exonic",
            "tss_match",
            "tes_match"
        ),
        artifact = c(
            "intron_shift",
            "exon_misalignment",
            "fake_terminal_exon_5",
            "fake_terminal_exon_3",
            "terminal_exon_misalignment_5",
            "terminal_exon_misalignment_3",
            "exon_elongation_5",
            "exon_elongation_3",
            "fake_micro_intron_retention"
        ),
        intron = c(
            "intron_retention",
            "unspliced_intron_retention",
            "incomplete_intron_retention_5",
            "incomplete_intron_retention_3"
        ),
        inconsistent = c(
            "major_exon_elongation_5",
            "major_exon_elongation_3",
            "extra_intron_5",
            "extra_intron_3",
            "extra_intron",
            "alt_donor_site",
            "alt_acceptor_site",
            "intron_migration",
            "intron_alternation",
            "mutually_exclusive_exons",
            "exon_skipping",
            "exon_merge",
            "exon_gain",
            "exon_detach",
            "terminal_exon_shift",
            "alternative_structure"
        ),
        alt.tss.tes = c(
            "alternative_polya_site",
            "internal_polya_site",
            "correct_polya_site",
            "aligned_polya_tail",
            "alternative_tss"
        )
    )
    
    ifelse(
        identical(event, "all"), 
        return(unname(unlist(list0))), 
        return(list0[[event]])
    )
}



# private
bed_to_feat <- function(data_bed) {

    result <- rphast::feat(seqname = data_bed$V1,
                   start = data_bed$V2,
                   end = data_bed$V3,
                   score = data_bed$V5,
                   strand = data_bed$V6,
                   attribute = data_bed$V4)

    return(result)
}



# --------------------------------------------------------------------
# public functions
# --------------------------------------------------------------------

conserved_elements_computation <- function(alignment_id, clade, chrs,
                                           required_feats_func) {



    chr <- chrs[1]

    # debug:
    # chr <- 22
    # alignment_id <- '100_way'
    # clade <- 'sarcopterygii'
    # required_feats_func <- required_species_features_sarcopterygii_100way

    config <- load_config()
    filtering_config <- config$conservation$filtering[[alignment_id]]
    sequence_names <- filtering_config[[clade]]

    align <- load_multiz_alignment(alignment_id, chr, sequence_names)

    neutral_model <- load_neutral_model(alignment_id)

    elements <- rphast::phastCons(align, neutral_model, expected.length=45,
                          target.coverage = 0.3, rho = 0.3, viterbi = TRUE)

    required_feats <- required_feats_func(align)

    conserved_elements <- rphast::coverage.feat(
        elements$most.conserved, required_feats,
        or = FALSE, get.feats = TRUE)

    conserved_elements$seqname <- paste0('chr', chr)

    # save to bed file

    out_file_name <- conserved_mostConserved_file_name(chr)
    out_file_path <- save_conserved_elements(elements = conserved_elements,
                        alignment_id = alignment_id, clade = clade,
                        output_file_name = out_file_name)

    rm()
    gc()


    # TODO
    if(is.na(chrs)) {

        #puede venir un vector de chrs

    } else {




    }


    #debug:
    return(out_file_path)
}



# conserved_path_1 <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/conservation/mammals/chr22_mostConserved.bed'
# conserved_path_2 <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/conservation/sarcopterygii/chr22_mostConserved.bed'
conserved_elements_in_common <- function(conserved_path_1, conserved_path_2,
                                         alignment_id,
                                         clade_1, clade_2, chr) {

    bed_1 <- read.delim(conserved_path_1, sep = '\t', header = FALSE)
    bed_2 <- read.delim(conserved_path_2, sep = '\t', header = FALSE)

    feat_1 <- bed_to_feat(bed_1)
    feat_2 <- bed_to_feat(bed_2)

    regions_in_common <- rphast::coverage.feat(feat_1, feat_2,
                                          or = FALSE, get.feats = TRUE)

    out_file_name <- conserved_mostConserved_in_common_file_name(chr)
    out_file_path <- save_conserved_elements(elements = regions_in_common,
                         alignment_id = alignment_id,
                         clade = paste(clade_1, clade_2, sep = "_"),
                         output_file_name = out_file_name)

    return(out_file_path)

}


# dado chr
# 1. load alignment
# 2. load neutral model
# 3. compute phastCons mostConserved
# 4. intersect with req features
# 5. save bed files

# revisar que las especies que figuran en configuracion son las que estan en el texto del paper
# ok para 100 way





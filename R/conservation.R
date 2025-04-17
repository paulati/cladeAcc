# private

conserved_elements_computation_chr <- function(alignment_id, clade, chr,
                                    sequence_names, required_feats_func) {

    align <- load_multiz_alignment(alignment_id, chr, sequence_names)

    neutral_model <- load_neutral_model(alignment_id, step = "cons")

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
                                             alignment_id = alignment_id,
                                             clade = clade,
                                             output_file_name = out_file_name)


    return(out_file_path)
}




# --------------------------------------------------------------------
# public functions
# --------------------------------------------------------------------


#' Compute conservation for clade species
#' @param alignment_id char One of these possible values (100_way, 77_way)
#' @param clade to complete
#' @param chrs description
#' @param required_feats_func delegate to function
#' @export
conserved_elements_computation <- function(alignment_id, clade, chrs = NA,
                                           required_feats_func) {

    if(anyNA(chrs)) {

        # TODO: process all
        out_file_paths <- NA

    } else {

        out_file_paths <- lapply(chrs, function(chr) {

            # debug:
            # chr <- 17
            # alignment_id <- '100_way'
            # clade <- 'sarcopterygii'
            # required_feats_func <- required_species_features_sarcopterygii_100way

            config <- load_config()
            filtering_config <- config$conservation$filtering[[alignment_id]]
            sequence_names <- filtering_config[[clade]]

            conserved_elements_file_path <-
                conserved_elements_computation_chr(alignment_id, clade, chr,
                                                   sequence_names, required_feats_func)

            rm()
            gc()

            return(conserved_elements_file_path)

        })

        names(out_file_paths) <- chrs

    }

    return(out_file_paths)
}



# conserved_path_1 <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/conservation/mammals/chr22_mostConserved.bed'
# conserved_path_2 <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/conservation/sarcopterygii/chr22_mostConserved.bed'


#' Compute intersection between elements sets
#' @param alignment_id char One of these possible values (100_way, 77_way)
#' @param clade_1 to complete
#' @param clade_2 to complete
#' @param chr description
#' @param conserved_path_1 description
#' @param conserved_path_2 description
#' @export
conserved_elements_in_common <- function(alignment_id,
                                         clade_1, clade_2, chr,
                                         conserved_path_1 = NA,
                                         conserved_path_2 = NA) {

    if(! is.na(conserved_path_1)) {
        bed_1 <- read.delim(conserved_path_1, sep = '\t', header = FALSE)
        feat_1 <- bed_to_feat(bed_1)
    } else {
        feat_1 <- load_conserved_elements(alignment_id, clade_1, chr)
    }

    if(! is.na(conserved_path_2)) {
        bed_2 <- read.delim(conserved_path_2, sep = '\t', header = FALSE)
        feat_2 <- bed_to_feat(bed_2)
    } else {
        feat_2 <- load_conserved_elements(alignment_id, clade_2, chr)
    }

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






# private functions

split_candidate_regions <- function(alignment_id, clade, chr, feat_length) {


    base_path <- acc_candidate_elements_path(alignment_id, clade, feat_length)
    if(! dir.exists(base_path)) {
        dir.create(base_path, recursive = TRUE)
    }
    file_name <- paste0("chr", chr, "_", clade, "_feat_",  feat_length, ".csv")
    file_path <- file.path(base_path, file_name)

    if(file.exists(file_path)) {

        data <- read.delim(file_path, sep = '\t', header = TRUE)
        result <- bed_to_feat(data)

    } else {

        candidate_regions <- load_conserved_elements(alignment_id,
                                                     clade, chr,
                                                     in_common = TRUE)

        result <- rphast::split.feat(candidate_regions, f = feat_length,
                                     drop = TRUE, pointer.only = FALSE)

        write.table(result, file_path, sep="\t", col.names = TRUE,
                    row.names = FALSE, quote = FALSE)

    }

    return(result)

}



# branch <- branch_label
compute_observed_phyloP <- function(align, neutral_model_labeled, branch,
                            candidate_regions_feat) {

    phyloP_obs_result <- rphast::phyloP(neutral_model_labeled, msa = align,
                                mode = "ACC",
                                features = candidate_regions_feat,
                                branches = branch)

    # sort by score
    result <- phyloP_obs_result[order(-phyloP_obs_result$score), ]

    return(result)
}

# TODO: cambiar este nombre process_acc_obs_phyloP_computation
# candidate_regions_clade values: 'mammals_sarcopterygii' or 'aves_sarcopterygii'
# candidate_regions_clade <- 'mammals_sarcopterygii'
process_acc_obs_phyloP_computation <- function(alignment_id, alignment, acc_clade,
                                               candidate_regions_clade,
                                               chr, branch_label,
                                               neutral_model_labeled,
                                               feat_lengths = c(25, 50)) {



    obs_phyloP_paths <- lapply(feat_lengths, function(x) {

            candidate_regions_feat <- split_candidate_regions(alignment_id,
                                                        candidate_regions_clade,
                                                        chr,
                                                        feat_length = x)

            elements <- compute_observed_phyloP(alignment,
                                        neutral_model_labeled,
                                        branch_label,
                                        candidate_regions_feat)

            out_file_name <- accelerated_observed_phyloP_file_name(chr, x)
            # TODO: ver si tiene sentido guardar en carpetas separadas 25 y 50
            # si no tienen sentido, sacar x de la llamada a save_accelerated_elements
            elements_file_path <- save_obs_phyloP_elements(elements,
                                                       alignment_id, acc_clade,
                                                       x, out_file_name)




            rm()
            gc()

            return(elements_file_path)

        })

    names(obs_phyloP_paths) <- paste0('len_', feat_lengths)

    return(obs_phyloP_paths)

}


compute_accelerated_elements_chr <- function(alignment_id, acc_clade,
                                         candidate_regions_clade,
                                         chr, sequence_names, branch_label,
                                         neutral_model_labeled,
                                         feat_lengths = c(25, 50)) {

    alignment <- load_multiz_alignment(alignment_id, chr, sequence_names)

    obs_phyloP_result_paths <- process_acc_obs_phyloP_computation(
        alignment_id, alignment, acc_clade,
        candidate_regions_clade,
        chr, branch_label,
        neutral_model_labeled,
        feat_lengths = feat_lengths)

    nrep <- 100 # TODO: cambiar a 100 mil

    # non_parametric_stats based on chr, ...
    non_parametric_stats_paths  <- calculate_non_parametric_stats(
        alignment_id, alignment,
        neutral_model_labeled,
        acc_clade, branch_label,
        candidate_regions_clade,
        chr, feat_lengths, nrep)


    # non_parametric_significance pval FDR
    non_parametric_significance_paths <- calculate_non_parametric_significance(
        obs_phyloP_result_paths, non_parametric_stats_paths,
        alignment_id, acc_clade, chr)




}

# ------------------------

# public functions

# alignment_id <- alignment_id_mam
# acc_clade <- acc_clade_mam
# candidate_regions_clade <- candidate_regions_clade_mam
# chrs <- chrs_mam
# branch_label <- branch_label_mam
# label_neutral_model_func <- label_neutral_model_func_mam
# feat_lengths <- c(25, 50)

#' Compute accelerated elements ...
#' @param alignment_id char One of these possible values (100_way, 77_way)
#' @param acc_clade to complete
#' @param candidate_regions_clade to complete
#' @param chrs description
#' @param branch_label description
#' @param label_neutral_model_func description
#' @param feat_lengths description
#' @export
compute_accelerated_elements <- function(alignment_id, acc_clade,
                                             candidate_regions_clade,
                                             chrs, branch_label,
                                             label_neutral_model_func,
                                             feat_lengths = c(25, 50)) {

    # same neutral model for every chr

    neutral_model <- load_neutral_model(alignment_id, step = "acc")
    neutral_model_labeled <- label_neutral_model_func(neutral_model)
    # write neutral model as debug info
    save_labeled_neutral_model(alignment_id, neutral_model_labeled)


    config <- load_config()
    filtering_config <- config$acceleration$filtering[[alignment_id]]
    sequence_names <- filtering_config[[acc_clade]]

    # alignent for chr
    acc_elements_paths <- lapply(chrs, function(chr) {

        compute_accelerated_elements_chr(alignment_id, acc_clade,
                                                candidate_regions_clade,
                                                chr, sequence_names, branch_label,
                                                neutral_model_labeled,
                                                feat_lengths = c(25, 50))

    })

    names(acc_elements_paths) <- chrs

    return(acc_elements_paths)

}

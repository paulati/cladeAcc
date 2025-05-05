# TODO: check NA de local base path - DONE
# para todas las carpetas checkear si existe el dir y crearlo (lo hago en input output)

# ----------------------------------------------------------
# tmp dir

pkg_data_tmp_base_path <- function() {
    result <- tools::R_user_dir("cladeAcc", which = "data")
    return(result)
}

user_local_data_base_path <- function() {

    if(exists("user_base_path")) {
        if(!dir.exists(user_base_path)) {
            dir.create(user_base_path, showWarnings = FALSE, recursive = TRUE)
        }
        result <- user_base_path
    } else {
        result <- NA
    }
    return(result)
}

# ----------------------------------------------------------
# multiple alignment

multiz_alignment_relative_path <- function(alignment_id) {
    result <- file.path(alignment_id, 'raw')
    return(result)
}

multiz_alignment_paths <- function(alignment_id, chr) {

    file_name <- paste0('chr', chr, '.maf')
    storage_relative_path <- multiz_alignment_relative_path(alignment_id)

    local_storage_base_path <- user_local_data_base_path()
    if(is.na(local_storage_base_path)) {

        result_path_base <- pkg_data_tmp_base_path()
    } else  {

        result_path_base <- local_storage_base_path
    }

    result_path <- file.path(result_path_base,
                             storage_relative_path,
                             file_name)

    result_gz_path <- paste0(result_path, '.gz')

    result <- list('path' = result_path,
                   'path_gz' = result_gz_path)

    return(result)

}

fasta_alignment_base_paths <- function(alignment_id, clade, feat_length, chr) {

    storage_relative_path <- fasta_alignment_relative_path(alignment_id,
                                                           clade, feat_length)


    local_storage_base_path <- user_local_data_base_path()

    if(is.na(local_storage_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()
    } else {
        storage_base_path <- local_storage_base_path
    }

    base_path <- file.path(storage_base_path,
                               storage_relative_path$path,
                               paste0('chr', chr))


    base_path_gz <- file.path(storage_base_path,
                                  storage_relative_path$path_gz,
                                  paste0('chr', chr))

    result <- list(
        'path' = base_path,
        'path_gz' = base_path_gz
    )

    return(result)

}

# alignment_id <- '77_way'
# feat_length <- 25
# clade <- 'aves'
# chr <- 24
# start <- 10
# end <- 35
fasta_alignment_paths <- function(alignment_id, clade, feat_length,
                                  chr, start, end) {

    file_name <- paste0('acc_', feat_length, '_chr', chr, '_', start, '_',
                        end, '.fasta')

    # ejemplo
    # /home/rstudio/disco_tmp/data/acc_maf_ingroup/50/chr2/chr2_acc_50_15656775_15656825.maf

    base_path <- fasta_alignment_base_paths(alignment_id, clade, feat_length,
                                            chr)

    fasta_alignment_path <- file.path(base_path$path, file_name)
    fasta_alignment_gz_path <- file.path(base_path$path_gz,
                                         paste0(file_name, '.gz'))

    result <- list('path' = fasta_alignment_path,
                   'path_gz' = fasta_alignment_gz_path)

    return(result)

}


fasta_alignment_relative_path <- function(alignment_id, clade, len) {
    out_base <- file.path(alignment_id, 'output', 'acceleration', clade)
    out <- file.path(out_base, 'fasta_alignments', len)
    out_gz <- file.path(out_base, 'fasta_alignments_gz', len)
    result <- list(
        'path' = out,
        'path_gz' = out_gz
    )
    return(result)
}


# ----------------------------------------------------------
# consensus alignment / sequences


fasta_alignment_consensus_sequence_paths <- function(alignment_id, clade,
                                                     feat_length, chr) {

    file_name <- paste0('chr', chr, '_consensus_', feat_length,  '.txt')
    file_name_gz <- paste0(file_name,  '.gz')

    storage_relative_path <-
        fasta_alignment_consensus_sequence_relative_path(alignment_id, clade,
                                                         feat_length)

    # ejemplo
    # outgroup_consensus_file_name_head: "chr"
    # outgroup_consensus_file_name_tail: "_25_consensus.txt"

    local_storage_base_path <- user_local_data_base_path()

    if(is.na(local_storage_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()
    } else {
        storage_base_path <- local_storage_base_path
    }

    file_path <- file.path(storage_base_path,
                           storage_relative_path$path,
                           file_name)
    file_gz_path <- file.path(storage_base_path,
                              storage_relative_path$path_gz,
                              file_name_gz)

    result <- list('path' = file_path,
                   'path_gz' = file_gz_path)

    return(result)


}

fasta_alignment_consensus_sequence_relative_path <- function(alignment_id,
                                                             clade, len) {
    out_base <- file.path(alignment_id, 'output', 'acceleration', clade,
                          'fasta_alignments', len)
    out <- file.path(out_base, 'consensus_sequence')

    out_gz <- file.path(out_base, 'consensus_sequence_gz')

    result <- list(
        'path' = out,
        'path_gz' = out_gz
    )

    return(result)
}

# ----------------------------------------------------------
# scoring
# acc raw scoring, acc filtered scoring, cons raw scoring, raw filtered scoring

acc_raw_scoring_relative_path <- function(alignment_id, clade, len) {

    out <- file.path(alignment_id, 'output', 'acceleration', clade,
                        'score_raw', len)

    out_gz <- file.path(alignment_id, 'output', 'acceleration', clade,
                        'score_raw_gz', len)

    result <- list(
        'path' = out,
        'path_gz' = out_gz
    )

    return(result)

}


# score_elements_filtered_base_path <- function(alignment_id, clade,
#                                              feat_length) {
custom_filtering_base_path <- function(alignment_id, clade) {


    local_base_path <- user_local_data_base_path()

    if(is.na(local_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()
    } else {
        storage_base_path <- local_base_path
    }

    result <- file.path(storage_base_path, alignment_id,
                            'output', 'custom_filtering', clade)

    return(result)
}


custom_filtering_relative_path <- function(alignment_id, clade, len) {

    out_base <- file.path(alignment_id, 'output', 'custom_filtering',
                          clade, 'acceleration')
    out <- file.path(out_base, 'score_filtered', len)
    out_gz <- file.path(out_base, 'score_filtered_gz', len)

    result <- list(
        'path' = out,
        'path_gz' = out_gz
    )

    return(result)
}


acc_raw_scoring_paths <- function(alignment_id, clade, feat_length, chr){

    file_name <- paste0('chr', chr, '_score_', feat_length, '.csv')
    file_name_gz <- paste0(file_name, '.gz')
    storage_relative_path <- acc_raw_scoring_relative_path(alignment_id,
                                                           clade, feat_length)

    local_storage_base_path <- user_local_data_base_path()

    if(is.na(local_storage_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()
    } else {
        storage_base_path <- local_storage_base_path
    }

    file_path <- file.path(storage_base_path,
                               storage_relative_path$path,
                               file_name)
    file_gz_path <- file.path(storage_base_path,
                                  storage_relative_path$path_gz,
                                  file_name_gz)

    result <- list('path' = file_path,
                   'path_gz' = file_gz_path)

    return(result)
}

acc_filtered_scoring_path <- function(alignment_id, clade, feat_length, chr){

    out_file_name <- paste0('chr', chr, '_score_', feat_length,
                            '_filtered_norm.csv')
    out_file_name_gz <- paste0(out_file_name, '.gz')

    relative_path <- custom_filtering_relative_path(alignment_id, clade,
                                                    feat_length)

    local_base_path <- user_local_data_base_path()

    if(is.na(local_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()
    } else {
        storage_base_path <- local_base_path
    }

    out_base_path <- file.path(storage_base_path, relative_path$path)
    out_base_path_gz <- file.path(storage_base_path, relative_path$path_gz)

    out_file_path <- file.path(out_base_path, out_file_name)
    out_file_path_gz <- file.path(out_base_path_gz, out_file_name_gz)


    result <- list('path' = out_file_path,
                   'path_gz' = out_file_path_gz)

    return(result)

}



# ----------------------------------------------------------
# neutral model

neutral_model_relative_path <- function(alignment_id) {
    result <- file.path(alignment_id, 'neutral_model')
    return(result)
}

# step in c("cons", "acc")
neutral_model_paths <- function(alignment_id, step = "cons") {

    config <- load_config()
    if(step == "cons") {
        step_config <- config$conservation
    } else if(step == "acc") {
        step_config <- config$acceleration
    } else {
        print ("error: step values should be 'cons' or 'acc'")
        return()
    }

    file_name <- step_config$neutral_model[[alignment_id]]$file_name
    storage_relative_path <- neutral_model_relative_path(alignment_id)

    local_storage_base_path <- user_local_data_base_path()

    if(is.na(local_storage_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()

    } else {
        storage_base_path <- local_storage_base_path
    }

    result <- file.path(storage_base_path, storage_relative_path, file_name)

    return(result)

}


# ----------------------------------------------------------
# conservation

conserved_elements_path <- function(alignment_id, clade) {

    out_base_tail <- file.path(alignment_id, 'output', 'conservation', clade)

    local_base_path <- user_local_data_base_path()

    if(is.na(local_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()
    } else {
        storage_base_path <- local_base_path
    }

    result  <- file.path(storage_base_path, out_base_tail)

    return(result)
}

conserved_elements_filtered_path <- function(alignment_id, clade){

    base_path <- conserved_elements_path(alignment_id, clade)
    result <- file.path(base_path, 'filtered_elements')
    return(result)
}

conserved_mostConserved_file_name <- function(chr) {
    result <- glue::glue("chr{chr}_mostConserved.bed")
    return(result)
}

conserved_mostConserved_in_common_file_name <- function(chr) {
    result <- glue::glue("chr{chr}_mostConserved_in_common.bed")
    return(result)
}

conserved_filtered_file_name <- function(chr) {
    result <- glue::glue("chr{chr}_mostConserved_in_common_filtered.bed")
    return(result)
}

# ----------------------------------------------------------
# candidate elements

candidate_elements_raw_scoring_path <- function(alignment_id, clade,
                                                feat_length, chr) {

    base_path <- custom_filtering_base_path(alignment_id, clade)
    out_base_path <- file.path(base_path, 'candidate_elements', 'score_raw',
                               feat_length)
    file_name <- paste0('chr', chr, '_score_', feat_length, '.csv')
    result <- file.path(out_base_path, file_name)
    return(result)
}



# ----------------------------------------------------------
# acceleration


accelerated_elements_path <- function(alignment_id, clade) {

    local_base_path <- user_local_data_base_path()

    common_acc_base_path <- file.path(alignment_id, 'output',
                                      'acceleration', clade)

    if(is.na(local_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()
    } else {
        storage_base_path <- local_base_path
    }

    acc_base_path <- file.path(storage_base_path, common_acc_base_path)

    result <- list(
        observed_phyloP = file.path(acc_base_path, 'observed_phyloP'),
        non_parametric_phyloP_distribution =
            file.path(acc_base_path, 'non_parametric_distribution_phyloP'),
        non_parametric_phyloP =
            file.path(acc_base_path, 'non_parametric_phyloP'))

    return(result)

}

acc_candidate_elements_path <- function(alignment_id, clade, len) {

    common_acc_base_path <- file.path(alignment_id, 'output',
                                      'acceleration', clade, 'candidates', len)

    local_base_path <- user_local_data_base_path()

    if(is.na(local_base_path)) {
        storage_base_path <- pkg_data_tmp_base_path()
    } else {
        storage_base_path <- local_base_path
    }

    result <- file.path(storage_base_path, common_acc_base_path)

    return(result)

}

accelerated_observed_phyloP_file_name <- function(chr, split_length) {
    result <- glue::glue("chr{chr}_obsPhyloP_{split_length}.csv")
    return(result)
}

accelerated_non_parametric_phyloP_distribution_file_name <- function(
        chr, split_length) {
    result <- glue::glue("chr{chr}_non_parametric_phyloP_distribution_{split_length}.csv")
    return(result)
}

accelerated_non_parametric_phyloP_file_name <- function(
        chr, split_length) {
    result <- glue::glue("chr{chr}_non_parametric_phyloP_{split_length}.csv")
    return(result)
}








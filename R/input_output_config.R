# ----------------------------------------------------------
# tmp dir

pkg_data_tmp_base_path <- function() {
    result <- tools::R_user_dir("cladeAcc", which = "data")
    return(result)
}

user_data_base_path <- function() {

    if(exists(user_base_path)) {
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

# TODO: terminar!!
multiz_alignment_paths <- function(alignment_id, chr) {

    file_name <- paste0('chr', chr, '.maf')

    tmp_storage_base_path <- pkg_data_tmp_base_path()
    tmp_storage_relative_path <- multiz_alignment_relative_path(alignment_id)
    tmp_multiz_alignment_path <- file.path(tmp_storage_base_path,
                                           tmp_storage_relative_path,
                                           file_name)
    tmp_multiz_alignment_gz_path <- paste0(tmp_multiz_alignment_path, '.gz')

    result <- list('tmp' = tmp_multiz_alignment_path,
                   'tmp_gz' = tmp_multiz_alignment_gz_path,
                   'local' = '',
                   'local_gz' = '',
                   'aws_gz' = '')

    return(result)

}

# TODO: hacer que esto devuelva tmp, local, aws
fasta_alignment_base_paths <- function(alignment_id, clade, feat_length, chr) {

    tmp_storage_base_path <- pkg_data_tmp_base_path()
    tmp_storage_relative_path <- fasta_alignment_relative_path(alignment_id,
                                                             clade, feat_length)
    base_path_tmp <- file.path(tmp_storage_base_path,
                                          tmp_storage_relative_path,
                                          paste0('chr', chr))

    result <- list('tmp' = base_path_tmp,
                   'tmp_gz' = '',
                   'local' = '',
                   'local_gz' = '',
                   'aws_gz' = '')

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

    tmp_fasta_alignment_path <- file.path(base_path$tmp, file_name)
    tmp_fasta_alignment_gz_path <- paste0(tmp_fasta_alignment_path, '.gz')

    result <- list('tmp' = tmp_fasta_alignment_path,
                   'tmp_gz' = tmp_fasta_alignment_gz_path,
                   'local' = '',
                   'local_gz' = '',
                   'aws_gz' = '')

    return(result)

}


fasta_alignment_relative_path <- function(alignment_id, clade, len) {
    result <- file.path(alignment_id, 'output', 'acceleration', clade,
                        'fasta_alignments', len)

    return(result)
}


# ----------------------------------------------------------
# consensus alignment / sequences


fasta_alignment_consensus_sequence_paths <- function(alignment_id, clade,
                                                     feat_length, chr) {

    file_name <- paste0('chr', chr, '_consensus_', feat_length,  '.txt')

    # ejemplo
    # outgroup_consensus_file_name_head: "chr"
    # outgroup_consensus_file_name_tail: "_25_consensus.txt"

    tmp_storage_base_path <- pkg_data_tmp_base_path()
    tmp_storage_relative_path <-
        fasta_alignment_consensus_sequence_relative_path(alignment_id, clade,
                                                         feat_length)
    tmp_file_path <- file.path(tmp_storage_base_path,
                               tmp_storage_relative_path,
                               file_name)
    tmp_file_gz_path <- paste0(tmp_file_path, '.gz')

    result <- list('tmp' = tmp_file_path,
                   'tmp_gz' = tmp_file_gz_path,
                   'local' = '',
                   'local_gz' = '',
                   'aws_gz' = '')

    return(result)


}

fasta_alignment_consensus_sequence_relative_path <- function(alignment_id,
                                                             clade, len) {
    result <- file.path(alignment_id, 'output', 'acceleration', clade,
                        'fasta_alignments', len, 'consensus_sequence')

    return(result)
}

# ----------------------------------------------------------
# scoring
# acc raw scoring, acc filtered scoring, cons raw scoring, raw filtered scoring

acc_raw_scoring_relative_path <- function(alignment_id, clade, len) {

    result <- file.path(alignment_id, 'output', 'acceleration', clade,
                        'score_raw', len)
    return(result)
}


# score_elements_filtered_base_path <- function(alignment_id, clade,
#                                              feat_length) {
custom_filtering_base_path <- function(alignment_id, clade) {

    base_path <- pkg_data_tmp_base_path()
    result <- file.path(base_path, alignment_id,
                        'output', 'custom_filtering', clade)
    return(result)
}

# custom_filtering_acc_base_path <- function(alignment_id, clade,
#                                            feat_length) {
#
#     base_path <- custom_filtering_base_path(alignment_id, clade,
#                                             feat_length)
#
#
# }
#
# custom_filtering_cons_base_path <- function(){
#
#
# }


acc_raw_scoring_paths <- function(alignment_id, clade, feat_length, chr){

    file_name <- paste0('chr', chr, '_score_', feat_length, '.csv')

    tmp_storage_base_path <- pkg_data_tmp_base_path()
    tmp_storage_relative_path <- acc_raw_scoring_relative_path(alignment_id,
                                                                clade, feat_length)
    tmp_file_path <- file.path(tmp_storage_base_path,
                               tmp_storage_relative_path,
                               file_name)
    tmp_file_gz_path <- paste0(tmp_file_path, '.gz')

    result <- list('tmp' = tmp_file_path,
                   'tmp_gz' = tmp_file_gz_path,
                   'local' = '',
                   'local_gz' = '',
                   'aws_gz' = '')

    return(result)
}

acc_filtered_scoring_path <- function(alignment_id, clade, feat_length, chr){

    base_path <- custom_filtering_base_path(alignment_id, clade)
    out_base_path <- file.path(base_path, 'acceleration', 'score_filtered',
                               feat_length)
    if(!dir.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }
    out_file_name <- paste0('chr', chr, '_score_', feat_length,
                            '_filtered_norm.csv')
    out_file_path <- file.path(out_base_path, out_file_name)

    result <- list('tmp' = out_file_path,
                   'tmp_gz' = paste0(out_file_path, '.gz'),
                   'local' = '',
                   'local_gz' = '',
                   'aws_gz' = '')

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

    tmp_storage_base_path <- pkg_data_tmp_base_path()
    tmp_storage_relative_path <- neutral_model_relative_path(alignment_id)
    tmp_neutral_model_path <- file.path(tmp_storage_base_path,
                                        tmp_storage_relative_path,
                                        file_name)

    # TODO: completar!!

    result <- list('tmp' = tmp_neutral_model_path,
                   'local' = '',
                   'aws_gz' = '')

    return(result)

}


# ----------------------------------------------------------
# conservation


# TODO: ver que esto puede estar en una serie de directorios, no solo tmp
conserved_elements_path <- function(alignment_id, clade) {

    base_path <- pkg_data_tmp_base_path()
    result <- file.path(base_path, alignment_id,
                        'output', 'conservation', clade)
    return(result)
}

# TODO: ver que esto puede estar en una serie de directorios, no solo tmp
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
    if(!dir.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    file_name <- paste0('chr', chr, '_score_', feat_length, '.csv')
    result <- file.path(out_base_path, file_name)
    return(result)
}



# ----------------------------------------------------------
# acceleration

# TODO: ver que esto puede estar en una serie de directorios, no solo tmp
accelerated_elements_path <- function(alignment_id, clade) {

    base_path <- pkg_data_tmp_base_path()
    acc_base_path <- file.path(base_path, alignment_id,
                               'output', 'acceleration', clade)

    result <- list(
        observed_phyloP = file.path(acc_base_path, 'observed_phyloP'),
        non_parametric_phyloP_distribution =
            file.path(acc_base_path, 'non_parametric_distribution_phyloP'),
        non_parametric_phyloP =
            file.path(acc_base_path, 'non_parametric_phyloP'))

    return(result)

}

acc_candidate_elements_path <- function(alignment_id, clade, len) {

    base_path <- pkg_data_tmp_base_path()
    result <- file.path(base_path, alignment_id,
                        'output', 'acceleration', clade, 'candidates', len)
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









pkg_data_tmp_base_path <- function() {
    result <- tools::R_user_dir("cladeAcc", which = "data")
    return(result)
}

multiz_alignment_relative_path <- function(alignment_id) {
    result <- file.path(alignment_id, 'raw')
    return(result)
}

neutral_model_relative_path <- function(alignment_id) {
    result <- file.path(alignment_id, 'neutral_model')
    return(result)
}


conserved_elements_path <- function(alignment_id, clade) {

    base_path <- pkg_data_tmp_base_path()
    result <- file.path(base_path, alignment_id,
                        'output', 'conservation', clade)
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

neutral_model_paths <- function(alignment_id) {

    config <- load_config()
    file_name <- config$conservation$neutral_model[[alignment_id]]$file_name

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

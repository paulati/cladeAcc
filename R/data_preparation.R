
# private functions

# alignment_id <- '100_way'
# chr <- 22
# alignment_id in {100_way, 77_way}
download_multiz_file <- function(alignment_id, chr) {

    config <- load_config()
    alignment_source_config <- config$data_preparation$source[[alignment_id]]
    url <- file.path(alignment_source_config$base_url,
                  glue::glue(alignment_source_config$file_name_pattern))
    # tmp_storage_relative_path <- file.path(alignment_id, 'raw')
    tmp_storage_relative_path <- multiz_alignment_relative_path(alignment_id)
    file_path <- download_data(url, tmp_storage_relative_path,
                               force_download = FALSE)
    return(file_path)
}

download_neutral_model_file <- function(alignment_id) {

    config <- load_config()
    neutral_model_config <- config$conservation$neutral_model[[alignment_id]]
    url <- file.path(neutral_model_config$base_url, neutral_model_config$file_name)
    tmp_storage_relative_path <- neutral_model_relative_path(alignment_id)
    file_path <- download_data(url, tmp_storage_relative_path,
                               force_download = FALSE)
    return(file_path)
}



get_md5_multiz_data <- function(alignment_id, force_download = FALSE) {

    config <- load_config()

    alignments_config <- config$data_preparation$source[[alignment_id]]

    url <- file.path(alignments_config$base_url, alignments_config$md5_file_name)

    tmp_storage_relative_path <- multiz_alignment_relative_path(alignment_id)
    file_path <- download_data(url, tmp_storage_relative_path,
                               force_download)

    md5_data_tmp <- read.delim(file_path, sep = " ", header = FALSE)
    data_col_names <- c('value', 'file')
    md5_data <- md5_data_tmp |> dplyr::select(V1, V3) |>
        setNames(dplyr::all_of(data_col_names))

    return(md5_data)
}


get_md5_neutral_model_data <- function(alignment_id, force_download = FALSE) {

    config <- load_config()
    neutral_model_config <- config$conservation$neutral_model[[alignment_id]]
    url <- file.path(neutral_model_config$base_url,
                     neutral_model_config$md5_file_name)
    tmp_storage_relative_path <- neutral_model_relative_path(alignment_id)
    file_path <- download_data(url, tmp_storage_relative_path,
                               force_download)

    md5_data_tmp <- read.delim(file_path, sep = " ", header = FALSE)
    data_col_names <- c('value', 'file')
    md5_data <- md5_data_tmp |> dplyr::select(V1, V3) |>
        setNames(dplyr::all_of(data_col_names))

    return(md5_data)

}

check_md5 <- function(file_path, md5_data) {


    colnames(md5_data) <- c('md5', 'file_id')

    md5_value <- tools::md5sum(file_path)
    file_name <- basename(file_path)
    target_value <- md5_data |>
        dplyr::filter(file_id == file_name) |>
        dplyr::select(md5) |>
        dplyr::pull()
    ok <- md5_value == target_value

    return(ok)

}

save_file_to_local_storage <- function(tmp_path) {


}

save_file_to_aws_storage <- function(tmp_path) {


}

# --------------------------------------------------------------------
# public functions
# --------------------------------------------------------------------

#' clean_data_preparation_tmp_folder
#' @param clean_tmp boolean Files and folders in local tmp will be deleted
#' @param clean_local_storage boolean Files and folders in local storage
#' will be deleted
#' @param clean_aws_storage boolean Files and folders in aws storage
#' will be deleted
#' @export
clean_data_preparation_folders <- function(clean_tmp = TRUE,
                                           clean_local_storage = FALSE,
                                           clean_aws_storage = FALSE) {

    if(clean_tmp) {
        # tmp_base_dir <- tools::R_user_dir("cladeAcc", which = "data")
        tmp_base_dir <- pkg_data_tmp_base_path()
        # delete all files and folder in this directory
        unlink(tmp_base_dir, recursive = TRUE)
    }

    if(clean_local_storage) {

        # TODO
    }

    if(clean_aws_storage) {

        # TODO
    }

}

#' @title prepare_alignment_files
#' @description Prepare alignment files for specified chromosomes
#' @param alignment_id char One of these possible values (100_way, 77_way)
#' @param chrs chr Vector containing alignment chromosome ids
#' @param local_storage boolean Save alignments to local storage
#' (default value FALSE)
#' @param aws_storage boolean Save alignments to aws storage
#' (default value FALSE)
#' @export
prepare_alignment_files <- function(alignment_id, chrs, local_storage = FALSE,
                          aws_storage = FALSE) {

    # reuse the same md5_data object for all chrs
    md5_data <- get_md5_multiz_data(alignment_id)

    lapply(chrs, function(x) prepare_alignment_file(alignment_id, x, md5_data,
                                                    local_storage, aws_storage))

}

#' Prepare alignment file for specified chromosome
#' @param alignment_id char One of these possible values (100_way, 77_way)
#' @param chr chr alignment chromosome id
#' @param md5_data data.frame Data frame containing target md5 values
#' (default value NA)
#' @param local_storage boolean Save alignments to local storage
#' (default value FALSE)
#' @param aws_storage boolean Save alignments to aws storage
#' (default value FALSE)
#' @export
prepare_alignment_file <- function(alignment_id, chr, md5_data = NA,
                         local_storage = FALSE, aws_storage = FALSE) {

    # debug:
    # alignment_id <- '100_way'
    # chr <- 17
    print(alignment_id)


    #download data
    file_path <- download_multiz_file(alignment_id, chr)

    # check md5
    if (anyNA(md5_data)) {
        md5_data <- get_md5_multiz_data(alignment_id)
    }

    download_ok <- check_md5(file_path, md5_data)

    print(download_ok)

    if(download_ok) {

        # save file local / aws


    } else {
        msg <- paste0('error downloading chr' , chr, ' ',  alignment_id,
                      ' alignment.')
        print(msg)
    }

}


#' Prepare neutral model for specified alignment
#' @param alignment_id char One of these possible values (100_way, 77_way)
#' @param md5_data data.frame Data frame containing target md5 values
#' (default value NA)
#' @param local_storage boolean Save alignments to local storage
#' (default value FALSE)
#' @param aws_storage boolean Save alignments to aws storage
#' (default value FALSE)
#' @export
prepare_neutral_model_file <- function(alignment_id, md5_data = NA,
                                local_storage = FALSE, aws_storage = FALSE){

    #download data
    file_path <- download_neutral_model_file(alignment_id)

    # check md5
    if (anyNA(md5_data)) {
        md5_data <- get_md5_neutral_model_data(alignment_id)
    }

    download_ok <- check_md5(file_path, md5_data)

    if(download_ok) {

        # save file local / aws


    } else {
        msg <- paste0('error downloading chr' , chr, ' ',  alignment_id,
                      ' alignment.')
        print(msg)
    }



}



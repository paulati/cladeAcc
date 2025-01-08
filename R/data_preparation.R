library(tools)
library(dplyr)


# where to save data?
# https://stackoverflow.com/questions/47101382/where-to-put-r-files-that-generate-package-data
# https://r-pkgs.org/data.html

#----------------

# https://www.ascend.io/blog/how-to-build-a-data-pipeline-in-six-steps/
# data pipeline
# https://www.simplilearn.com/what-is-data-processing-article
# data processing
# data preparation  https://www.talend.com/resources/what-is-data-preparation/
# https://aws.amazon.com/what-is/data-preparation/
# Collect data. Collecting data is the process of assembling all the data you need for ML. ...
# Clean data. Cleaning data corrects errors and fills in missing data as a step to ensure data quality. ...
# Label data. ...
# Validate and visualize.

# https://rstudio.github.io/config/articles/introduction.html#usage

# https://cran.r-project.org/web/packages/config/vignettes/config.html

# config package

# https://www.appsilon.com/post/r-config

# private functions

options(timeout = max(300, getOption("timeout")))

# generic function to download data from specified url
download_data <- function(url,
                          tmp_storage_relative_path = '',
                          force_download = FALSE) {

    tmp_base_dir <- tools::R_user_dir("cladeAcc", which = "data")
    dir.create(tmp_base_dir, recursive = TRUE, showWarnings = TRUE)

    if(nchar(tmp_storage_relative_path) == 0) {
        data_base_dir <- tmp_base_dir
    } else {
        data_base_dir <- file.path(tmp_base_dir, tmp_storage_relative_path)
        dir.create(data_base_dir, recursive = TRUE, showWarnings = TRUE)
    }

    setwd(data_base_dir)

    file_name <- basename(url)

    data_file_path <- file.path(data_base_dir, file_name)
    if(file.exists(data_file_path)) {
        if(force_download) {
            download.file(url, file_name)
        }
    } else {

        download.file(url, file_name)
    }

    local_data_file_path <- file.path(data_base_dir, file_name)

    return(local_data_file_path)
}

# alignment_id <- '100_way'
# chr <- 22
# alignment_id in {100_way, 77_way}
download_file <- function(alignment_id, chr) {

    config <- load_config()
    alignment_source_config <- config$data_preparation$source[[alignment_id]]
    url <- file.path(alignment_source_config$base_url,
                  glue::glue(alignment_source_config$file_name_pattern))
    tmp_storage_relative_path <- file.path(alignment_id, 'raw')
    file_path <- download_data(url, tmp_storage_relative_path,
                               force_download = FALSE)
    return(file_path)
}

get_md5_data <- function(alignment_id) {

    config <- load_config()
    alignment_source_config <- config$data_preparation$source[[alignment_id]]
    url <- file.path(alignment_source_config$base_url,
                     alignment_source_config$md5_file_name)
    tmp_storage_relative_path <- file.path(alignment_id, 'raw')
    file_path <- download_data(url, tmp_storage_relative_path,
                               force_download = FALSE)

    md5_data_tmp <- read.delim(file_path, sep = " ", header = FALSE)
    data_col_names <- c('value', 'file')
    md5_data <- md5_data_tmp |> select(V1, V3) |> setNames(all_of(data_col_names))

    return(md5_data)
}


check_md5 <- function(file_path, md5_data) {

    # debug:
    # file_path <- "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/chr22.maf.gz"
    # alignment_id <- '100_way'
    # md5_data <- get_md5_data(alignment_id)

    md5_value <- md5sum(file_path)
    file_name <- basename(file_path)
    target_value <- md5_data |>
        filter(file == file_name) |>
        select(value) |>
        pull()
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
        tmp_base_dir <- tools::R_user_dir("cladeAcc", which = "data")
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
#' (default value TRUE)
#' @export
prepare_alignment_files <- function(alignment_id, chrs, local_storage = FALSE,
                          aws_storage = TRUE) {

    # reuse the same md5_data object for all chrs
    md5_data <- get_md5_data(alignment_id)

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
#' (default value TRUE)
#' @export
prepare_alignment_file <- function(alignment_id, chr, md5_data = NA,
                         local_storage = FALSE, aws_storage = TRUE) {

    # debug:
    # alignment_id <- '100_way'
    # chr <- 22


    #download data
    file_path <- download_file(alignment_id, chr)

    # check md5
    if (is.na(md5_data)) {
        md5_data <- get_md5_data(alignment_id)
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

# --------------------------------------------------------------------




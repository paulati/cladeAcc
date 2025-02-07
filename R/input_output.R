
options(timeout = max(300, getOption("timeout")))

# generic function to download data from specified url
download_data <- function(url,
                          tmp_storage_relative_path = '',
                          force_download = FALSE) {

    #tmp_base_dir <- tools::R_user_dir("cladeAcc", which = "data")
    tmp_base_dir <- pkg_data_tmp_base_path()
    dir.create(tmp_base_dir, recursive = TRUE, showWarnings = FALSE)

    if(nchar(tmp_storage_relative_path) == 0) {
        data_base_dir <- tmp_base_dir
    } else {
        data_base_dir <- file.path(tmp_base_dir, tmp_storage_relative_path)
        dir.create(data_base_dir, recursive = TRUE, showWarnings = FALSE)
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


#--------------------------------------------------
# specific functions:

load_multiz_alignment <- function(alignment_id, chr, sequence_names = NULL) {

    # check order: tmp dir, local dir, aws dir

    paths <- multiz_alignment_paths(alignment_id, chr)
    tmp_file_path <- paths$tmp
    tmp_file_path_gz <- paths$tmp_gz
    local_file_path <- paths$local
    local_file_path_gz <- paths$local_gz
    aws_file_path <- paths$aws_gz

    use_this_path <- ''

    if(file.exists(tmp_file_path)) {
        use_this_path <- tmp_file_path
    } else if(file.exists(tmp_file_path_gz)) {

        # unzip:
        gunzip(filename = tmp_file_path_gz,
               destname = tmp_file_path,
               remove = FALSE)

        use_this_path <- tmp_file_path
    } else if(file.exists(local_file_path)) {
        use_this_path <- local_file_path

    } else if(file.exists(local_file_path_gz)) {
        # TODO
    } else if(file.exists(aws_file_path)) {

        # TODO
    }


    #setwd(local.base.folder.path)

    # if(! file.exists(align_local_file_name)){
    #     if(file.exists(file_name_gz))
    #     {
    #         gunzip(local.file_name_gz, align_local_file_name, remove = FALSE)
    #
    #     }
    #     else
    #     {
    #         download.from.s3(account.key, account.secret,
    #                          remote.base.folder.path, file_name_gz,
    #                          local.base.folder.path, file_name_gz,
    #                          align_local_file_name, bucket.name)
    #     }
    # }

    align <- rphast::read.msa(use_this_path, seqnames = sequence_names)

    return(align)


}

load_neutral_model <- function(alignment_id) {

    # check order: tmp dir, local dir, aws dir

    paths <- neutral_model_paths(alignment_id)
    tmp_file_path <- paths$tmp
    local_file_path <- paths$local
    aws_file_path <- paths$aws_gz

    use_this_path <- ''

    if(file.exists(tmp_file_path)) {
        use_this_path <- tmp_file_path
    } else if(file.exists(local_file_path)) {
        use_this_path <- local_file_path
    } else if(file.exists(aws_file_path)) {

        # TODO
    }

    result <- rphast::read.tm(use_this_path)

    return(result)

}


# elements <- conserved_elements
# alignment_id
# clade
# chr
save_conserved_elements <- function(elements, alignment_id, clade,
                                    output_file_name) {

    output_base_path <- conserved_elements_path(alignment_id, clade)
    if(!dir.exists(output_base_path)) {
        dir.create(output_base_path, recursive = TRUE)
    }
    output_file_path <- file.path(output_base_path, output_file_name)


    names <- str_replace(elements$attribute, pattern = "id \"([0-9]+)\"",
                replacement = "id_\\1")

    q <- nrow(elements)
    bed_data <- data.frame(chr = elements$seqname,
                           start = elements$start,
                           end = elements$end,
                           name = names,
                           score = elements$score,
                           strand = elements$strand)

    write.table(bed_data, output_file_path, sep="\t",
                col.names = FALSE, row.names = FALSE, quote=FALSE)


    return(output_file_path)
}

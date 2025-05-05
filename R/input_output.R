
# TODO: reescribir todos los save para que llamen a un save base
# los mismo para todos los loads que solo leen un dataframe

options(timeout = max(300, getOption("timeout")))


prepare_base_dir <- function() {

    custom_user_base_path <- user_local_data_base_path()

    if(is.na(custom_user_base_path)) {
        #tmp_base_dir <- tools::R_user_dir("cladeAcc", which = "data")
        result <- pkg_data_tmp_base_path()
    } else {
        result <- custom_user_base_path
    }

    print(result)

    if(! file.exists(result)) {
        dir.create(result, recursive = TRUE, showWarnings = FALSE)
    }

    return(result)

}



# generic function to download data from specified url
download_data <- function(url,
                          storage_relative_path = '',
                          force_download = FALSE) {


    base_dir <- prepare_base_dir()

    print(base_dir)

    if(nchar(storage_relative_path) == 0) {
        data_base_dir <- base_dir
    } else {
        data_base_dir <- file.path(base_dir, storage_relative_path)
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


bed_to_feat <- function(data_bed) {


    if('V1' %in% colnames(data_bed)) {

        result <- rphast::feat(seqname = data_bed$V1,
                           start = data_bed$V2,
                           end = data_bed$V3,
                           score = data_bed$V5,
                           strand = data_bed$V6,
                           attribute = data_bed$V4)

    } else {

        result <- rphast::feat(seqname = data_bed$seqname,
                               start = data_bed$start,
                               end = data_bed$end,
                               score = data_bed$score,
                               strand = data_bed$strand,
                               attribute = data_bed$attribute)

    }

    return(result)
}

#--------------------------------------------------
# specific functions:

load_multiz_alignment <- function(alignment_id = NA, chr = NA,
                                  sequence_names = NULL,
                                  file_path = NA) {

    if(is.na(file_path)) {

        paths <- multiz_alignment_paths(alignment_id, chr)
        data_file_path <- paths$path
        data_file_path_gz <- paths$path_gz

        use_this_path <- ''

        if(file.exists(data_file_path)) {
            use_this_path <- data_file_path

        } else if(file.exists(data_file_path_gz)) {
            # unzip:
            R.utils::gunzip(filename = data_file_path_gz,
                            destname = data_file_path,
                            remove = FALSE)

            use_this_path <- data_file_path

        }

    } else {

        use_this_path <- file_path

    }

    align <- rphast::read.msa(use_this_path, seqnames = sequence_names,
                              pointer.only = FALSE)

    return(align)


}

# step in "cons" "acc"
load_neutral_model <- function(alignment_id, step = "cons") {

    fle_path <- neutral_model_paths(alignment_id, step)

    result <- rphast::read.tm(fle_path)

    return(result)

}


save_labeled_neutral_model <- function(alignment_id, neutral_model) {

    path <- neutral_model_paths(alignment_id, step = "acc")

    file_name <- basename(path)
    base_path <- dirname(path)

    if(! file.exists(base_path)) {
        dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
    }

    neutral_model_labeled_file_name <- stringr::str_replace(file_name,
                                                   pattern = "\\.mod",
                                                   replacement = ".labeled.mod")
    neutral_model_labeled_file_path <- file.path(base_path,
                                                neutral_model_labeled_file_name)

    labeled_neutral_model_file_path <-
        rphast::write.tm(neutral_model, neutral_model_labeled_file_path,
                         append = FALSE)

    return(neutral_model_labeled_file_path)
}


# elements <- conserved_elements
# alignment_id
# clade
# chr

save_conserved_elements <- function(elements, alignment_id, clade,
                                    output_file_name, filtered = FALSE,
                                    format = 'feat') {

    if(filtered) {
        output_base_path <- conserved_elements_filtered_path(alignment_id,
                                                             clade)
    } else {
        output_base_path <- conserved_elements_path(alignment_id, clade)
    }

    print(output_base_path)

    if(!file.exists(output_base_path)) {
        dir.create(output_base_path, recursive = TRUE)
    }
    output_file_path <- file.path(output_base_path, output_file_name)

    q <- nrow(elements)
    elements_data_cols_in_bed <- c('seqnames', 'start', 'end')

    if(is.null(elements$seqnames)) {
        seqnames <- elements$seqname
    } else {
        seqnames <- elements$seqnames
    }

    bed_data <- data.frame(chr = seqnames,
                           start = elements$start,
                           end = elements$end)

    if(! is.null(elements$attribute)) {
        names <- stringr::str_replace(elements$attribute,
                                      pattern = "id \"([0-9]+)\"",
                                      replacement = "id_\\1")
        bed_data$names <- names
        elements_data_cols_in_bed <- c(elements_data_cols_in_bed, 'names')
    } else {

        bed_data$names <- ''
    }

    if(! is.null(elements$score)) {
        bed_data$score <- elements$score
        elements_data_cols_in_bed <- c(elements_data_cols_in_bed, 'score')
    } else {

        bed_data$score <- ''
    }

    if(! is.null(elements$strand)) {
        bed_data$strand <- elements$strand
        elements_data_cols_in_bed <- c(elements_data_cols_in_bed, 'strand')
    } else {
        bed_data$strand <- '*'
    }

    extra_columns <- setdiff(colnames(elements), elements_data_cols_in_bed)
    for(x in extra_columns) {
        bed_data[, x] <- elements[, x]
    }

    write_col_names <- format != 'feat'
    write.table(bed_data, output_file_path, sep="\t",
                col.names = write_col_names, row.names = FALSE, quote=FALSE)


    return(output_file_path)
}

# file_path <- cons_elem_path
load_conserved_elements <- function(alignment_id = NA, clade = NA, chr = NA,
                                    in_common = FALSE, filtered = FALSE,
                                    file_path = NA, format = "feat") {

    if(is.na(file_path)) {

        if(filtered) {
            base_path <- conserved_elements_filtered_path(alignment_id,
                                                          clade)
        } else {
            base_path <- conserved_elements_path(alignment_id, clade)
        }

        if(in_common) {
            file_name <- conserved_mostConserved_in_common_file_name(chr)
        } else {
            file_name <- conserved_mostConserved_file_name(chr)
        }

        file_path <- file.path(base_path, file_name)

    }

    read_col_names <- format != 'feat'
    data <- read.delim(file_path, sep = '\t', header = read_col_names)

    if(format == "feat") {
        result <- bed_to_feat(data)
    } else {
        result <- data
    }


    # doesn't work: score columns should contain integer
    # result <- rphast::read.feat(file_path, pointer.only = FALSE)

    return(result)

}

# TODO: cambiar a obs_phylop en lugar de accelerated
# save_accelerated_elements
save_obs_phyloP_elements <- function(elements, alignment_id, clade,
                                      feat_length, output_file_name) {

    acc_base_path <- accelerated_elements_path(alignment_id, clade)

    output_base_path <- file.path(acc_base_path$observed_phyloP, feat_length)
    if(!file.exists(output_base_path)) {
        dir.create(output_base_path, recursive = TRUE, showWarnings = FALSE)
    }
    output_file_path <- file.path(output_base_path, output_file_name)

    write.table(elements, output_file_path,
                col.names = TRUE, row.names = FALSE,
                sep = "\t", quote = FALSE)

    return(output_file_path)
}

load_obs_phyloP_elements <- function(file_path = NA, alignment_id = '',
                                     clade = '', feat_length = 25,
                                     file_name = '') {

    if(is.na(file_path)) {
        acc_base_path <- accelerated_elements_path(alignment_id, clade)
        file_path <- file.path(acc_base_path$observed_phyloP, feat_length,
                               file_name)
    }

    data <- read.delim(file_path, sep = '\t', header = TRUE)

    return(data)

}

save_non_parametric_phyloP_elements <- function(elements, alignment_id, clade,
                                                feat_length, chr) {


    file_name <- accelerated_non_parametric_phyloP_file_name(chr, feat_length)
    acc_base_path <- accelerated_elements_path(alignment_id, clade)
    output_base_path <- file.path(acc_base_path$non_parametric_phyloP,
                                  feat_length)
    if(!file.exists(output_base_path)) {
        dir.create(output_base_path, recursive = TRUE, showWarnings = FALSE)
    }
    file_path <- file.path(output_base_path, file_name)

    write.table(elements, file_path, col.names = TRUE, row.names = FALSE,
                sep = "\t", quote = FALSE)

    return(file_path)

}

save_non_parametric_phyloP_distribution <- function(data, alignment_id, clade,
                                                    feat_length, chr) {

    output_file_name <- accelerated_non_parametric_phyloP_distribution_file_name(
        chr, feat_length)

    base_path <- accelerated_elements_path(alignment_id, clade)

    output_base_path <- file.path(base_path$non_parametric_phyloP_distribution,
                                  feat_length)

    if(!file.exists(output_base_path)) {
        dir.create(output_base_path, recursive = TRUE)
    }
    output_file_path <- file.path(output_base_path, output_file_name)

    write.table(data,
                file = output_file_path,
                sep = '\t', col.names = TRUE, row.names = FALSE,
                quote = FALSE)

    return(output_file_path)

}


load_non_parametric_phyloP_distribution <- function(file_path = NA,
                                                    alignment_id = '',
                                                    clade = '', chr = '',
                                                    feat_length = 25) {


    if(is.na(file_path)) {

        base_path <- accelerated_elements_path(alignment_id, clade)

        file_name <- accelerated_non_parametric_phyloP_distribution_file_name(
            chr, feat_length)

        file_path <- file.path( base_path$non_parametric_phyloP_distribution,
                                feat_length, file_name)
    }

    if (file.exists(file_path)) {

        data <- read.delim(file_path, sep = '\t', header = TRUE)

    } else {

        data <- data.frame()
    }

    result <- list('data' = data, 'path' = file_path)

    return(result)
}

load_non_parametric_phyloP_elements <- function(file_path = NA,
                                                alignment_id = '',
                                                clade = '', chr = '',
                                                feat_length = 25) {

    if(is.na(file_path)) {

        base_path <- accelerated_elements_path(alignment_id, clade)

        file_name <- accelerated_non_parametric_phyloP_file_name(
            chr, feat_length)

        file_path <- file.path( base_path$non_parametric_phyloP,
                                feat_length, file_name)
    }

    if (file.exists(file_path)) {

        data <- read.delim(file_path, sep = '\t', header = TRUE)

    } else {

        data <- data.frame()
    }

    result <- list('data' = data, 'path' = file_path)

    return(result)

}

save_elements_consensus_info <- function(data, alignment_id,
                                         clade, feat_length, chr) {

    out_file_path <- fasta_alignment_consensus_sequence_paths(alignment_id,
                                                    clade, feat_length, chr)
    out_base_path <- dirname(out_file_path$path)

    if(! file.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    write.table(data, file = out_file_path$path, sep = '\t',
                col.names = TRUE, row.names = FALSE, quote = FALSE)


    return(out_file_path$path)

}

load_elements_consensus_info <- function(alignment_id, clade, feat_length,
                                         chr) {

    file_path <- fasta_alignment_consensus_sequence_paths(alignment_id,
                                                    clade, feat_length, chr)

    data <- read.delim(file_path$path, sep = '\t', header = TRUE)

    return(data)
}

# ------------------------------------------------------
# scroring / custom filtering

save_acc_raw_scoring <- function(data, alignment_id, clade, feat_length, chr) {

    file_path <- acc_raw_scoring_paths(alignment_id, clade, feat_length,
                                           chr)

    out_base_path <- dirname(file_path$path)
    if(! file.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    write.table(data, file = file_path$path, sep = '\t', col.names = TRUE,
                row.names = FALSE, quote = FALSE)

    return(file_path)

}



load_acc_raw_scoring <- function(file_path) {

    data <- read.delim(file_path, sep = '\t', header = TRUE)
    return(data)

}

save_acc_filtered_scoring <- function(data, alignment_id,
                                      clade, feat_length, chr) {

    out_file_path <- acc_filtered_scoring_path(alignment_id,
                                                    clade, feat_length, chr)
    out_base_path <- dirname(out_file_path$path)
    if(! file.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    #     /u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/scoring/25
    # chr22_score_25_filtered_elements_norm.csv
    # write.table(data_ingroup_sort, "join_filtered_elements_norm.csv", sep="\t",
    #             quote = FALSE, col.names = TRUE, row.names = FALSE)

    write.table(data, out_file_path$path, sep="\t",
                quote = FALSE, col.names = TRUE, row.names = FALSE)

    return(out_file_path$path)

}

load_acc_filtered_scoring <- function(file_path) {
    data <- read.delim(file_path, sep = '\t', header = TRUE)
    return(data)
}

# result from merge cons acc
save_candidate_elements_raw_scoring <- function(data, alignment_id, clade,
                                                feat_length, chr) {

    file_path <- candidate_elements_raw_scoring_path(alignment_id, clade,
                                        feat_length, chr)

    base_path <- dirname(file_path)
    if(! file.exists(base_path)) {
        dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
    }

    write.table(data, file = file_path, sep = '\t', col.names = TRUE,
                row.names = FALSE, quote = FALSE)

    return(file_path)

}

# result from merge cons acc
load_candidate_elements_raw_scoring <- function(file_path) {

    data <- read.delim(file_path, header = TRUE, sep = '\t')
    return(data)
}

# result from merge cons acc
save_candidate_elements_filtered_scoring <- function(data, alignment_id, clade,
                                                     feat_length, chr) {

    base_path <- custom_filtering_base_path(alignment_id, clade)

    out_base_path <- file.path(base_path, 'candidate_elements', 'score_filtered',
                               feat_length)
    if(!file.exists(out_base_path)) {
        dir.create(out_base_path, recursive = TRUE, showWarnings = FALSE)
    }

    file_name <- paste0('chr', chr, '_score_', feat_length, '.csv')
    file_path <- file.path(out_base_path, file_name)
    write.table(data, file = file_path, sep = '\t', col.names = TRUE,
                row.names = FALSE, quote = FALSE)

    return(file_path)

}

# result from merge cons acc
load_candidate_elements_filtered_scoring <- function() {

    data <- read.delim(file_path, header = TRUE, sep = '\t')
    return(data)

}


# revisar!!!
# save_filtered_score_elements <- function(data, alignment_id, clade,
#                                          feat_length, chr, step_name_suffix = '') {
#
#     base_path <- custom_filtering_base_path(alignment_id, clade,
#                                                    feat_length)
#     if(! dir.exists(base_path)) {
#         dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
#     }
#
#     file_name <- glue::glue('chr{chr}_score{step_name_suffix}.csv')
#     file_path <- file.path(base_path, file_name)
#
#     write.table(data, file = file_path, sep = '\t', col.names = TRUE,
#                 row.names = FALSE, quote = FALSE)
#
#     return(file_path)
#
#
# }


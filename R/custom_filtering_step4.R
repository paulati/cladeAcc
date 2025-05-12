

# alignment_id <- '100_way'
# ingroup_clade <- 'mammals'
# outgroup_clade <- 'sarcopterygii'
# feat_length <- 25
# chr <- 22


# step 4
filter_acc_scoring_len <- function(acc_scoring_path, alignment_id,
                                    clade, feat_length, chr) {


    acc_filtered_scoring_paths <- acc_filtered_scoring_path(alignment_id,
                                    clade, feat_length, chr)
    if(nchar(acc_scoring_path) == 0) {

        out_file_path <- ''

    } else {

        out_file_path <- acc_filtered_scoring_paths$path

        if(! file.exists(out_file_path)) {

            # read result from calculate_acc_scoring_data for specific chr
            data <- load_acc_raw_scoring(acc_scoring_path)

            #filter as join
            data_ingroup <- data |> dplyr::filter(hamming_in_out_excl_gaps > 0 &
                                                      acc_element_shift_count != 0)
            # combino dos pasos
            # el que se guardaria en join_filtered_elements.csv
            # y el normalize como en https://github.com/paulati/acc_regions_mammals_aves/blob/223eb53cf3f202baaf1323cc8718c2210e99dc32/aves/py_source/acc_scoring/analysis.R#L75

            data_ingroup$shift_count_rel <- data_ingroup$acc_element_shift_count /
                data_ingroup$acc_element_size

            # scale 0 1
            max_distance <- max(data_ingroup$hamming_in_out_excl_gaps)
            min_distance <- min(data_ingroup$hamming_in_out_excl_gaps)
            data_ingroup$hamming_in_out_excl_gaps_scl <-
                (data_ingroup$hamming_in_out_excl_gaps - min_distance) /
                (max_distance - min_distance)

            max_shift_count <- max(data_ingroup$shift_count_rel)
            min_shift_count <- min(data_ingroup$shift_count_rel)
            data_ingroup$shift_count_rel_scl <-
                (data_ingroup$shift_count_rel - min_shift_count) /
                (max_shift_count - min_shift_count)

            # plot(data_ingroup$shift_count_rel_scl, data_ingroup$hamming_in_out_scl)
            # hist(data_ingroup$shift_count_rel_scl)
            # plot(data_ingroup$shift_count_rel_scl)
            # hist(data_ingroup$hamming_in_out_scl)

            #max_distance <- data_ingroup[(data_ingroup$distance.rel > 0.8), ]
            #max_distance.shift <- max_distance[max_distance$shift.rel > 0.6, ]

            # calculo el modulo del vector formado por distance.rel y shift.rel
            data_ingroup$norm2_scl <- apply(
                data_ingroup[,
                    c("hamming_in_out_excl_gaps_scl", "shift_count_rel_scl")],
                MARGIN = 1,
                FUN = function(x) norm(x, type="2"))

            data_ingroup_sort <- data_ingroup[order(-data_ingroup$norm2_scl), ]


            out_file_path <- save_acc_filtered_scoring(data_ingroup_sort,
                                        alignment_id, clade, feat_length, chr)

            # result join_filtered_elements_norm.csv
            # normalizo los datos de los join y da como resultado join_filtered_elements_norm.csv segun este codigo
            # https://github.com/paulati/acc_regions_mammals_aves/blob/223eb53cf3f202baaf1323cc8718c2210e99dc32/aves/py_source/acc_scoring/analysis.R#L4

        }

    }

    return(out_file_path)

}

filter_acc_scoring <- function(acc_scoring_paths, alignment_id, ingroup_clade,
                               feat_lengths, chr) {

    result <- lapply(seq(1, length(feat_lengths)), function(x) {
        feat_length <- feat_lengths[x]
        acc_scoring_path <- acc_scoring_paths[[x]]
        filter_acc_scoring_len(acc_scoring_path$path, alignment_id,
                               ingroup_clade, feat_length, chr)
    })

    names(result) <- paste0('len_', feat_lengths)

    return(result)


}

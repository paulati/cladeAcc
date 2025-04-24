# step 3
calculate_acc_scoring_len <- function(alignment_id, ingroup_clade,
                                      outgroup_clade, feat_length, chr) {

    ingroup_consensus_data <- load_elements_consensus_info(alignment_id,
                                                           ingroup_clade,
                                                           feat_length, chr)

    ingroup_consensus_data$id <- basename(ingroup_consensus_data$path)

    outgroup_consensus_data <- load_elements_consensus_info(alignment_id,
                                                            outgroup_clade,
                                                            feat_length, chr)

    outgroup_consensus_data$id <- basename(outgroup_consensus_data$path)

    consensus_data <- dplyr::inner_join(ingroup_consensus_data,
                                        outgroup_consensus_data,
                                        by = 'id',
                                        suffix = c("_ingroup", "_outgroup"))

    shift_sequences <- get_shift_seqs(consensus_data$path_ingroup, alignment_id)

    shift_counts <- calculate_shifts_counts(consensus_data, shift_sequences)

    gap_counts <- calculate_gaps_counts(
        consensus_data$consensus_sequence_ingroup)

    hamming_ingroup_incl_gaps <- calculate_hamming_pair_distances(
        consensus_data$consensus_sequence_ingroup, shift_sequences)

    hamming_outgroup_incl_gaps <- calculate_hamming_pair_distances(
        consensus_data$consensus_sequence_outgroup, shift_sequences)

    # excluding gaps that are simultaneously shifts:
    data_clean <- clean_shifts_gaps_positions(consensus_data,
                                              shift_sequences)
    # hamming based on clean data
    hamming_ingroup_excl_gaps <- calculate_hamming_pair_distances(
        data_clean$clean_ingroup_consensus,
        data_clean$clean_ingroup_shift_sequence)

    hamming_outgroup_excl_gaps <-  calculate_hamming_pair_distances(
        data_clean$clean_outgroup_consensus,
        data_clean$clean_outgroup_shift_sequence)

    # length(hamming_ingroup_excl_gaps) == length(hamming_outgroup_excl_gaps)
    hamming_in_out_incl_gaps <- unlist(lapply(
        seq(1, length(hamming_ingroup_incl_gaps)), function(x){
            coord_x <- hamming_ingroup_incl_gaps[x]
            coord_y <- hamming_outgroup_incl_gaps[x]
            result <- calculate_distance_to_line_xeqy(coord_x, coord_y)
            return(result)
        }))

    hamming_in_out_excl_gaps <- unlist(lapply(
        seq(1, length(hamming_ingroup_excl_gaps)), function(x){
            coord_x <- hamming_ingroup_excl_gaps[x]
            coord_y <- hamming_outgroup_excl_gaps[x]
            result <- calculate_distance_to_line_xeqy(coord_x, coord_y)
            return(result)
        }))

    #----
    # format data and return

    acc_elements_info <- parse_element_ids(consensus_data$id)

    result <- data.frame('acc_element_bed_chr' = acc_elements_info$chr,
                 'acc_element_bed_start' = acc_elements_info$start,
                 'acc_element_bed_end' = acc_elements_info$end,
                 'id' = consensus_data$path_ingroup,
                 'acc_element_len' = acc_elements_info$len,
                 'acc_element_size' = unlist(lapply(shift_sequences, nchar)),
                 'acc_element_gap_count' = gap_counts,
                 'acc_element_shift_count' = shift_counts,
                 'hamming_ingroup_incl_gaps' = hamming_ingroup_incl_gaps,
                 'hamming_outgroup_incl_gaps' = hamming_outgroup_incl_gaps,
                 'hamming_in_out_incl_gaps' = hamming_in_out_incl_gaps,
                 'hamming_ingroup_excl_gaps' = hamming_ingroup_excl_gaps,
                 'hamming_outgroup_excl_gaps' = hamming_outgroup_excl_gaps,
                 'hamming_in_out_excl_gaps' = hamming_in_out_excl_gaps

    )

    scoring_path <- save_acc_raw_scoring(result, alignment_id, ingroup_clade,
                                         feat_length, chr)

    return(scoring_path)

}



calculate_acc_scoring <- function(alignment_id, ingroup_clade,
                           outgroup_clade, feat_lengths = c(25, 50), chr) {


    result <- lapply(feat_lengths, function(feat_length) {

        calculate_acc_scoring_len(alignment_id, ingroup_clade, outgroup_clade,
                                              feat_length, chr)

    })

    names(result) <- paste0('len_', feat_lengths)

    return(result)

}

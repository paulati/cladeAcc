# alignment_id <- '100_way'
# clade <- 'mammals'
# chr <- 22
# step 3 aux:
get_shift_seqs <- function(paths, alignment_id) {


    config <- load_config()
    filtering_config <- config$custom_filtering$filtering[[alignment_id]]
    shift_specie <- filtering_config[['shift_specie']]

    shift_sequences <- unlist(lapply(paths, function(path) {

        alignment <- load_multiz_alignment(file_path = path,
                                           sequence_names = shift_specie)

        shift_sequence <- alignment$seqs

        return(shift_sequence)
    }))

    return(shift_sequences)

}


calculate_gaps_counts <- function(sequences) {
    result <- unlist(lapply(sequences,
                            function(sequence){
        gap_count <- stringr::str_count(sequence, pattern = "N|\\-")
        return(gap_count)
    }))

    return(result)
}

calculate_gaps_mask <- function(sequences) {
    result <- lapply(sequences,
                function(sequence){
                    sequence_parts <- stringr::str_split(sequence, pattern = "")
                    gap_mask <- as.integer(unlist(lapply(
                        sequence_parts,
                        function(x) x %in% c('N', '-'))))
                    return(gap_mask)
                })

    return(result)
}


# sequence1 <-   "TTXGAGGAGGTXTAXTTCCGXATCC"
# sequence2 <- "TTCGAGGAGGTCTACTTCCGCATCC"
calculate_shifts <- function(sequence1, sequence2) {

    sequence1_parts <- unlist(stringr::str_split(sequence1, pattern = ""))
    sequence2_parts <- unlist(stringr::str_split(sequence2, pattern = ""))

    shifts <- c()

    if(length(sequence1_parts) == length(sequence2_parts)) {

        shifts <- unlist(lapply(

            seq(1, length(sequence1_parts)), function(x) {

            seq1_base <- sequence1_parts[x]
            seq2_base <- sequence2_parts[x]

            if(
                (seq1_base %in% c('X', 'N', '-')) |
                (seq2_base %in% c('X', 'N', '-'))) {

                result <- 0

            } else {

                result <- as.integer(seq1_base != seq2_base)

            }

            return(result)
        }))

    } else {

        print("error")

    }

    return(shifts)

}

calculate_shifts_mask <- function(sequences1, sequences2) {

    if(length(sequences1) == length(sequences2)) {

        shift_mask <- lapply(seq(1, length(sequences1)), function(x){
            mask <- calculate_shifts(sequences1[x], sequences2[x])
            return(mask)
        })

    } else {

        print('error')
    }

    return(shift_mask)
}


calculate_shifts_counts <- function(consensus_data,
                                   shift_sequences) {

    # shift occurs when outgroup consensus is different from shift_sequence
    # and
    # ingroup consensus is equals to shift_sequence

    # any base in c( X, N, or -) is NOT a shift

    result <- unlist(lapply(seq(1, nrow(consensus_data)), function(x) {

        outgroup_consensus <- consensus_data$consensus_sequence_outgroup[x]
        ingroup_consensus <- consensus_data$consensus_sequence_ingroup[x]
        shift_sequence <- shift_sequences[x]

        ingroup_shifts <- calculate_shifts(ingroup_consensus, shift_sequence)
        outgroup_shifts <- calculate_shifts(outgroup_consensus, shift_sequence)

        position_is_shift <- unlist(lapply(seq(1, length(ingroup_shifts)),
                                           function(x) {

            result <- as.integer(
                (ingroup_shifts[x] == 0) & (outgroup_shifts[x] == 1))

            return(result)

        }))

        shift_count <- sum(position_is_shift)

        return(shift_count)

    }))

    return(result)


}


calculate_hamming_pair_distances <- function(sequences1, sequences2) {

    if(length(sequences1) == length(sequences2)) {

        result <- unlist(lapply(
            seq(1, length(sequences1)),
            function(x){
                elem_distance <- e1071::hamming.distance(
                    sequences1[x],
                    sequences2[x])
                return(elem_distance)
            }))

    } else {

        print('error')

        result <- c()
    }

    return(result)

}

get_substring_by_mask <- function(original_str, mask) {

    original_parts <- unlist(stringr::str_split(original_str, pattern = ""))
    clean_parts <- original_parts[mask]
    clean_str <- paste0(clean_parts, collapse = "")
    return(clean_str)

}


# exclude gaps that are simultaneously shifts:
clean_shifts_gaps_positions <- function(consensus_data,
                                        shift_sequences) {

    ingroup_gaps_mask <- calculate_gaps_mask(
        consensus_data$consensus_sequence_ingroup)
    ingroup_shifts_mask <- calculate_shifts_mask(
        consensus_data$consensus_sequence_ingroup, shift_sequences)

    outgroup_gaps_mask <- calculate_gaps_mask(
        consensus_data$consensus_sequence_outgroup)
    outgroup_shifts_mask <- calculate_shifts_mask(
        consensus_data$consensus_sequence_outgroup, shift_sequences)

    clean_ingroup_data <- lapply(seq(1, nrow(consensus_data)), function(x) {

        # ingroup
        gaps_ingroup_mask <- ingroup_gaps_mask[[x]]
        shift_ingroup_mask <- ingroup_shifts_mask[[x]]
        gap_and_shift_ingroup_mask <- gaps_ingroup_mask & shift_ingroup_mask

        original_ingroup_consensus <- consensus_data$consensus_sequence_ingroup[x]
        clean_ingroup_consensus <- get_substring_by_mask(
            original_ingroup_consensus, ! gap_and_shift_ingroup_mask)

        shift_sequence <- shift_sequences[[x]]
        clean_ingroup_shift_sequence <- get_substring_by_mask(
            shift_sequence, ! gap_and_shift_ingroup_mask)

        # outgroup
        gaps_outgroup_mask <- outgroup_gaps_mask[[x]]
        shift_outgroup_mask <- outgroup_shifts_mask[[x]]
        gap_and_shift_outgroup_mask <- gaps_outgroup_mask & shift_outgroup_mask

        original_outgroup_consensus <- consensus_data$consensus_sequence_outgroup[x]
        clean_outgroup_consensus <- get_substring_by_mask(
            original_outgroup_consensus, ! gap_and_shift_outgroup_mask)

        clean_outgroup_shift_sequence <- get_substring_by_mask(
            shift_sequence, ! gap_and_shift_outgroup_mask)



        result <- c(clean_ingroup_consensus, clean_ingroup_shift_sequence,
                    clean_outgroup_consensus, clean_outgroup_shift_sequence)

    })

    clean_ingroup_data_mtx <- do.call(rbind, clean_ingroup_data)
    clean_ingroup_data_df <- data.frame(clean_ingroup_data_mtx)
    colnames(clean_ingroup_data_df) <- c('clean_ingroup_consensus',
                                         'clean_ingroup_shift_sequence',
                                         'clean_outgroup_consensus',
                                         'clean_outgroup_shift_sequence')

    return(clean_ingroup_data_df)


}

parse_element_ids <- function(ids) {

    # id <- 'acc_25_chr22_33381950_33381974.fasta'

    pattern = 'acc_(?<len>[0-9]+)_(?<chr>chr[0-9]+)_(?<start>[0-9]+)_(?<end>[0-9]+)\\.fasta'

    len <- stringr::str_extract(ids, pattern, group = 1) # len
    chr <- stringr::str_extract(ids, pattern, group = 2) # chr
    start <- stringr::str_extract(ids, pattern, group = 3) # start
    end <- stringr::str_extract(ids, pattern, group = 4) # end

    result <- list(
        'len' = as.integer(len),
        'chr' = chr,
        'start' = as.integer(start),
        'end' = as.integer(end)
    )

    return(result)

}

# cross product (vecA x vecB) equals norm(vecA) * norm(vecB) * sin(alfa)
# coord_x <- 3
# coord_y <- 5
calculate_distance_to_line_xeqy <- function(coord_x, coord_y){

    if(coord_x > coord_y) {
        sign <-  -1
    } else {
        sign <- 1
    }

    point <- c(coord_x, coord_y)
    crossprod_result <- crossprod(point, c(1, 1))
    sin_alfa <- crossprod_result /
        (norm(c(1,1), type="2") * norm(point, type="2"))
    distance_result <- sin_alfa * norm(point, type="2")

    result <- distance_result * sign

    return(result)

}

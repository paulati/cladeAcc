# step 5 ?
# alignment_id <- '100_way'
# clade <- 'mammals_sarcopterygii'
# chr <- 22
calculate_conserved_elements <- function(alignment_id, clade, chr) {

    output_file_name <- conserved_filtered_file_name(chr)
    output_base_path <- conserved_elements_filtered_path(alignment_id,
                                                         clade)
    out_file_path <- file.path(output_base_path, output_file_name)

    if(! file.exists(out_file_path)) {

        data <- load_conserved_elements(alignment_id, clade, chr,
                                        in_common = TRUE)

        data_gr <- GenomicRanges::makeGRangesFromDataFrame(data)

        distance_threshold <- 20
        length_threshold <- 100
        data_gr_reduced <- GenomicRanges::reduce(data_gr,
                                            drop.empty.ranges=FALSE,
                                            min.gapwidth = distance_threshold,
                                            with.revmap=FALSE,
                                            with.inframe.attrib=FALSE)
        # min.gapwidth Ranges separated by a gap of at least min.gapwidth positions are not merged
        data_gr_reduced_length_filter <- data_gr_reduced[
            IRanges::width(data_gr_reduced) >= length_threshold]

        data_gr_reduced_length_filter$id <- paste0('chr', chr, '_',
            c(1:length(data_gr_reduced_length_filter)))

        data_gr_reduced_length_filter <- data.frame(data_gr_reduced_length_filter)

        out_file_path <- save_conserved_elements(data_gr_reduced_length_filter,
                                                 alignment_id, clade,
                                                 output_file_name,
                                                 filtered = TRUE,
                                                 format = 'data.frame')
    }
    # out.file.path <- file.path(base_path, "mammals/data/results",
    #                            "intersection_modneu_phastCons100way_sarcopterygii_mammals_functional_elements.bed")
    # # "/paula/2018/acc_regions/resultados/intersection_modneu_phastCons100way_sarcopterygii_mammals_functional_elements.bed"
    # write.table(gr.util.df, out.file.path, sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    #return(gr.util)

    return(out_file_path)



}



# alignment_id <- '100_way'
# cons_clade <- 'mammals_sarcopterygii'
# feat_length <- 25
# chr <- 22
# acc_scoring_path <- "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/custom_filtering/mammals/acceleration/score_filtered/25/chr114_score_25_filtered_norm.csv"
#
# step 5 (needs step 4 results join_filtered_elements_norm.csv)
calculate_cons_scoring_data <- function(alignment_id, cons_clade, feat_length,
                                        chr, acc_scoring_path)  {

    out_file_path <- candidate_elements_raw_scoring_path(alignment_id,
                                            cons_clade, feat_length, chr)

    if(! file.exists(out_file_path)) {


        cons_elem_path <- calculate_conserved_elements(alignment_id,
                                                       cons_clade,
                                                       chr)

        conservation_data <- load_conserved_elements(
            file_path = cons_elem_path, format = 'data.frame')
        conservation_gr <- GenomicRanges::makeGRangesFromDataFrame(
            conservation_data, keep.extra.columns = TRUE)

        acc_scoring_data <- load_acc_filtered_scoring(
            file_path = acc_scoring_path)
        acc_scoring_gr <- GenomicRanges::makeGRangesFromDataFrame(
            acc_scoring_data,
            seqnames.field = "acc_element_bed_chr",
            start.field = "acc_element_bed_start",
            end.field = "acc_element_bed_end",
            keep.extra.columns = TRUE)

        combined_data <- combine_data(acc_scoring_gr, conservation_gr)

        # how many acc elements in each conserved element?
        acc_elements_count_in_cons_elem <- combined_data |>
            dplyr::group_by(cons_id) |>
            dplyr::count(name = 'acc_elements_count_in_cons_elem')

        # how many shifts elements in each conserved element?
        shift_count_in_cons_elem <- combined_data |>
            dplyr::group_by(cons_id) |>
            dplyr::summarise(
                shift_count_in_cons_elem = sum(acc_acc_element_shift_count))

        tmp <- dplyr::inner_join(acc_elements_count_in_cons_elem,
                                 shift_count_in_cons_elem,
                                 by = "cons_id")

        gap_count_in_cons_elem <- combined_data |>
            dplyr::group_by(cons_id) |>
            dplyr::summarise(
                gap_count_in_cons_elem = sum(acc_acc_element_gap_count))

        tmp <- dplyr::inner_join(tmp, gap_count_in_cons_elem, by = "cons_id")

        # given a cons element, list of shift counts for each acc element in the cons element
        acc_elements_shift_counts <- combined_data |>
            dplyr::group_by(cons_id) |>
            dplyr::summarise(acc_elements_shift_counts =
                                paste0(acc_acc_element_shift_count, collapse = " "))

        tmp <- dplyr::inner_join(tmp, acc_elements_shift_counts, by = "cons_id")

        # given a cons element, list of gap counts for each acc element in the cons element
        acc_elements_gap_counts <- combined_data |>
            dplyr::group_by(cons_id) |>
            dplyr::summarise(acc_elements_gap_counts =
                                 paste0(acc_acc_element_gap_count, collapse = " "))

        tmp <- dplyr::inner_join(tmp, acc_elements_gap_counts, by = "cons_id")
        # build a table with data associated with each functional element (phastCons)
        #functional.element.count <- nrow(elements.count.functional.element)
        # functional.element.data <- data.frame( id = numeric(functional.element.count),
        #                                        acc_elements_count = numeric(functional.element.count),
        #                                        functional_element_shift_count = numeric(functional.element.count),
        #                                        functional_element_gap_count = numeric(functional.element.count),
        #                                        acc_elements_shift_counts = character(functional.element.count), # lista de shifts por elementos
        #                                        acc_elements_gap_counts = character(functional.element.count), # lista de gaps por elementos
        #                                        stringsAsFactors = FALSE )


        # colnames(functional.element.data) <- c("phastCons_id", "phastCons_seqnames",
        #                                        "phastCons_start", "phastCons_end",
        #                                        "phastCons_width", "phastCons_strand",
        #                                        "acc_elements_count",
        #                                        "functional_element_shift_count",
        #                                        "functional_element_gap_count",
        #                                        "acc_elements_shift_counts",
        #                                        "acc_elements_gap_counts" )

        all_data <- dplyr::inner_join(tmp, combined_data, by = 'cons_id')

        conservation_result <- data.frame("phastCons_id" = all_data$cons_id,
              "phastCons_seqnames" = all_data$cons_seqnames,
              "phastCons_start" = all_data$cons_start,
              "phastCons_end" = all_data$cons_end,
              "phastCons_width" = all_data$cons_width,
              "phastCons_strand" = all_data$cons_strand,
              "acc_elements_count" = all_data$acc_elements_count_in_cons_elem,
              "conserved_element_shift_count" = all_data$shift_count_in_cons_elem,
              "conserved_element_gap_count" = all_data$gap_count_in_cons_elem,
              "acc_elements_shift_counts" = all_data$acc_elements_shift_counts,
              "acc_elements_gap_counts" = all_data$acc_elements_gap_counts,
              "acc_id" = all_data$acc_id,
              stringsAsFactors = FALSE)

        combined_data_to_save_cols <- c("acc_seqnames", "acc_start", "acc_end",
            "acc_width", "acc_strand", "acc_id", "acc_acc_element_len",
            "acc_acc_element_size", "acc_acc_element_gap_count",
            "acc_acc_element_shift_count", "acc_hamming_ingroup_incl_gaps",
            "acc_hamming_outgroup_incl_gaps", "acc_hamming_ingroup_excl_gaps",
            "acc_hamming_outgroup_excl_gaps",
            "cons_id")

        acceleration_conservation_result <- dplyr::inner_join(
            conservation_result,
            combined_data[, combined_data_to_save_cols],
            by = c('phastCons_id' = 'cons_id', 'acc_id' = 'acc_id'))


        out_file_path <- save_candidate_elements_raw_scoring(
            acceleration_conservation_result, alignment_id, cons_clade,
            feat_length, chr)
    }
    return(out_file_path)

}


calculate_cons_scoring <- function(alignment_id, cons_clade, feat_lengths,
                                   chr, acc_scoring_paths) {

    result <- lapply(seq(1, length(feat_lengths)), function(x) {
        feat_length <- feat_lengths[x]
        acc_scoring_path <- acc_scoring_paths[[x]]
        calculate_cons_scoring_data(alignment_id, cons_clade, feat_length,
                                    chr, acc_scoring_path)
    })

    names(result) <- paste0('len_', feat_lengths)

    return(result)

}

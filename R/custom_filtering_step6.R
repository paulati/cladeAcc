
# step 6 filter cons soring data
filter_cons_scoring_len <- function(score_cons_raw_file_path, alignment_id,
                                     clade, feat_length, chr) {
    # result join_filtered_elements_norm_cod_nonCod_oneacczerogap_functionalregions.csv
    # https://github.com/paulati/acc_regions_mammals_aves/blob/223eb53cf3f202baaf1323cc8718c2210e99dc32/aves/r_source/functional_regions/_phastCons.analyze.R#L249

    data <- load_candidate_elements_raw_scoring(score_cons_raw_file_path)

    # calculo cant de shifts sin gaps por elementos conservados:
    new_columns <- lapply(seq(1, nrow(data)), function(i) {

        acc_elements_gap_count <- data$acc_elements_gap_counts[i]
        acc_elements_shift_count <- data$acc_elements_shift_counts[i]

        gaps_count <- as.numeric(unlist(
            stringr::str_split(acc_elements_gap_count, " ")))
        shifts_count <- as.numeric(unlist(
            stringr::str_split(acc_elements_shift_count, " ")))

        index_no_gaps <- which(gaps_count == 0)
        # print(index.no.gaps)

        if(length(index_no_gaps) > 0) {
            conserved_element_shift_count_exclgt0gap  <-
                sum(shifts_count[index_no_gaps])
            acc_elements_count_exclgt0gap <- length(index_no_gaps)
        } else {
            conserved_element_shift_count_exclgt0gap <- 0
            acc_elements_count_exclgt0gap <- 0
        }

        result <- c(conserved_element_shift_count_exclgt0gap,
                    acc_elements_count_exclgt0gap)

        return(result)

    })

    new_columns_df <- do.call(rbind, new_columns)
    colnames(new_columns_df) <- c('conserved_element_shift_count_exclgt0gap',
                                  'acc_elements_count_exclgt0gap')

    data <- cbind(data, new_columns_df)

    #conservo los registros que tienen al menos un acelerado sin gaps
    conserved_element_data_oneacczerogap <- data |>
        dplyr::filter(acc_elements_count_exclgt0gap > 0)

    conserved_element_data_oneacczerogap$conserved_element_shift_count_exclgt0gap_rel <-
        conserved_element_data_oneacczerogap$conserved_element_shift_count_exclgt0gap /
        conserved_element_data_oneacczerogap$phastCons_width

    # fitered.data$functional.element.shift.count.exclgt0gap_rel <- fitered.data$functional_element_shift_count_exclgt0gap / fitered.data$width
    conserved_element_data_oneacczerogap$acc_elements_count_exclgt0gap_rel <-
        conserved_element_data_oneacczerogap$acc_elements_count_exclgt0gap /
        conserved_element_data_oneacczerogap$phastCons_width

    # plot(fitered.data$functional.element.shift.count.exclgt0gap_rel, fitered.data$acc.elements.count.exclgt0gap_rel)

    # hist(fitered.data$acc.elements.count.exclgt0gap_rel)

    # guardo los datos de regiones funcionales (como resumen de los nonCod)
    #data.to.save <- merge(nonCod.data.functional, fitered.data,
    #                       by.x=c("phastCons_id"), by.y = c("id"))

    # colnames(data.to.save) <- c( "phastCons_id", "acc_element_bed_chr", "acc_element_bed_start",
    #                              "acc_element_bed_end", "id", "acc_element_len",
    #                              "acc_element_size", "acc_element_gap_count", "acc_element_shift_count",
    #                              "hamming_ingroup_incl_gaps", "hamming_outgroup_incl_gaps", "hamming_in_out_incl_gaps",
    #                              "hamming_ingroup_excl_gaps", "hamming_outgroup_excl_gaps", "hamming_in_out_excl_gaps",
    #                              "acc_element_table_browser", "acc_element_shift_count_rel", "hamming_in_out_excl_gaps_scl",
    #                              "shift_count_rel_scl", "norm2_scl", "coding",
    #                              "acc_elements_in_phastCons_count", "functional_element_shift_count", "functional_element_gap_count",
    #                              "acc_elements_in_phastCons_shift_counts", "acc_elements_in_phastCons_gap_counts",
    #                              "phastCons_seqnames", "phastCons_start", "phastCons_end", "phastCons_width",
    #                              "phastCons_strand", "functional_elements_shift_rel",  "functional_elements_acc_elements_count_rel",
    #                              "functional_element_shift_count_exclgt0gap", "acc_elements_count_exclgt0gap")


    #out.file.name <- "/paula/2018/acc_regions/scoring/data/results/201901/filter/join/join_filtered_elements_norm_cod_nonCod_zerogap_functionalregions.csv"
    # out.file.name <- "/paula/2018/acc_regions/scoring/data/results/201901/filter/join/join_filtered_elements_norm_cod_nonCod_oneacczerogap_functionalregions.csv"
    # write.table(data.to.save, out.file.name, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


    candidate_elements_file_path <- save_candidate_elements_filtered_scoring(
        conserved_element_data_oneacczerogap, alignment_id,
        clade, feat_length, chr)

    return(candidate_elements_file_path)

}


filter_cons_scoring <- function(score_cons_raw_file_paths, alignment_id,
                                    clade, feat_lengths, chr) {

    result <- lapply(seq(1, length(feat_lengths)), function(x) {
        feat_length <- feat_lengths[x]
        score_cons_raw_file_path <- score_cons_raw_file_paths[[x]]
        if(nchar(score_cons_raw_file_path) == 0) {
            cons_scoring_path <- ''
        } else {
            cons_scoring_path <- filter_cons_scoring_len(score_cons_raw_file_path, alignment_id,
                                            clade, feat_length, chr)
        }
        return(cons_scoring_path)
    })

    names(result) <- paste0('len_', feat_lengths)

    return(result)

}

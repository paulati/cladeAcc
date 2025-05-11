

# acc_scoring <- acc_scoring_gr
# cons_filtered <- conservation_gr
combine_data <- function(acc_scoring, cons_filtered) {

    # overlaps scores conserved elements:
    overlaps <- GenomicRanges::findOverlaps(query = acc_scoring,
                             subject = cons_filtered, type = "any")

    acc_scoring_indexes <- S4Vectors::queryHits(overlaps)
    cons_filtered_indexes <- S4Vectors::subjectHits(overlaps)

    acc_data <- data.frame(acc_scoring[acc_scoring_indexes])
    colnames(acc_data) <- paste0('acc_', colnames(acc_data))
    cons_data <- data.frame(cons_filtered[cons_filtered_indexes])
    colnames(cons_data) <- paste0('cons_', colnames(cons_data))

    result <- cbind(acc_data, cons_data)

    return(result)


    # # using overlaps queryHits to index acc.elements.data
    # acc.elements.data.functional <- acc.elements.data[queryHits(overlaps), ]
    # phastcons_ids <- gr.functional.elements[subjectHits(overlaps)]@
    #     elementMetadata$id
    # # add associated phastcons id to each acc element
    # acc.elements.data.functional$phastCons_id <- phastcons_ids
    #
    # functional.element.data <- get.acc.functional.element.data(
    #     acc.elements.data.functional, gr.functional.elements.df)



}

process_filtering_chr <- function(alignment_id, chr, ingroup_clade,
                                  outgroup_clade, debug = FALSE) {

    feat_lengths <- c(25, 50)

    print('processing step 1 / 6')

    elements_alignments_paths <- extract_elements_alignments_fasta(
        alignment_id, chr, ingroup_clade, outgroup_clade,
        split_lengths = c(25, 50))

    print('processing step 2 / 6')

    consensus_sequences_paths <-
        calculate_elements_consensus_sequences(alignment_id, chr, ingroup_clade,
                                    outgroup_clade, feat_lengths)


    print('processing step 3 / 6')

    acc_scoring_paths <- calculate_acc_scoring(alignment_id, ingroup_clade,
                            outgroup_clade, feat_lengths, chr)

    print('processing step 4 / 6')

    filter_acc_scoring_path <- filter_acc_scoring(acc_scoring_paths,
                                    alignment_id, ingroup_clade, feat_lengths,
                                    chr)

    print('processing step 5 / 6')

    cons_clade <- paste0(ingroup_clade, '_', outgroup_clade)
    score_cons_raw_file_paths <- calculate_cons_scoring(alignment_id,
                        cons_clade, feat_lengths, chr, acc_scoring_paths)

    print('processing step 6 / 6')

    score_cons_file_paths <- filter_cons_scoring(score_cons_raw_file_paths,
                                alignment_id, cons_clade, feat_lengths, chr)


    print(score_cons_file_paths)

}

# debug:
# alignment_id <- '100_way'
# clade <- 'mammals'
# chrs <- c(22)
# chr <- 22

#' complete description
#' @param alignment_id char One of these possible values (100_way, 77_way)
#' @param clade to complete
#' @param chrs description
#' @export
process_filtering <- function(alignment_id, clade, chrs) {

    if(clade == 'mammals' || clade == 'mammals_sarcopterygii') {

        ingroup_clade <- 'mammals'
        outgroup_clade <- 'sarcopterygii'


    } else if(clade == 'aves' || clade == 'aves_sarcopterygii') {

        ingroup_clade <- 'aves'
        outgroup_clade <- 'sarcopterygii'

    }

    lapply(chrs, function(chr) {
        process_filtering_chr(alignment_id, chr, ingroup_clade, outgroup_clade)
    })



}





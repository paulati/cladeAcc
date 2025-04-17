

# step 3 or 4 ?
# alignment_id <- '100_way'
# clade <- 'mammals_sarcopterygii'
# chr <- 22
calculate_conserved_elements <- function(alignment_id, clade, chr) {

    data <- load_conserved_elements(alignment_id, clade, chr, in_common = TRUE)

    data_gr <- GenomicRanges::makeGRangesFromDataFrame(data)

    distance_threshold <- 20
    length_threshold <- 100
    data_gr_reduced <- GenomicRanges::reduce(data_gr, drop.empty.ranges=FALSE,
                         min.gapwidth = distance_threshold,
                         with.revmap=FALSE, with.inframe.attrib=FALSE)
    # min.gapwidth Ranges separated by a gap of at least min.gapwidth positions are not merged
    data_gr_reduced_length_filter <- data_gr_reduced[
        IRanges::width(data_gr_reduced) >= length_threshold]

    data_gr_reduced_length_filter$id <- c(1:length(data_gr_reduced_length_filter))

    data_gr_reduced_length_filter <- data.frame(data_gr_reduced_length_filter)

    output_file_name <- conserved_filtered_file_name(chr)
    out_file_path <- save_conserved_elements(data_gr_reduced_length_filter,
                                        alignment_id, clade,
                                        output_file_name, filtered = TRUE,
                                        format = 'data.frame')

    # out.file.path <- file.path(base_path, "mammals/data/results",
    #                            "intersection_modneu_phastCons100way_sarcopterygii_mammals_functional_elements.bed")
    # # "/paula/2018/acc_regions/resultados/intersection_modneu_phastCons100way_sarcopterygii_mammals_functional_elements.bed"
    # write.table(gr.util.df, out.file.path, sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

    #return(gr.util)

    return(out_file_path)



}


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

    #TODO: hacer una funcion que me devuelva los paths para debug
    # y para usar de params en cada uno de los pasos siguentes
    # hacer algunos check para no recalcular lo que ya existe


    print(elements_alignments_paths)

    print('processing step 2 / 6')

    consensus_sequences_paths <-
        calculate_elements_consensus_sequences(alignment_id, chr, ingroup_clade,
                                    outgroup_clade, feat_lengths)

    if(debug) {

        consensus_sequences_paths <- list(
            'len_25' = list(
                'consensus_ingroup_file_path' = '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/fasta_alignments/25/consensus_sequence/chr17_consensus_25.txt',
                'consensus_outgroup_file_path' = '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/sarcopterygii/fasta_alignments/25/consensus_sequence/chr17_consensus_25.txt'
            ),
            'len_50' = list(
                'consensus_ingroup_file_path' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/fasta_alignments/50/consensus_sequence/chr17_consensus_50.txt",
                'consensus_outgroup_file_path' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/sarcopterygii/fasta_alignments/50/consensus_sequence/chr17_consensus_50.txt"
            )
        )
    }

    print(consensus_sequences_paths)

    print('processing step 3 / 6')

    acc_scoring_paths <- calculate_acc_scoring(alignment_id, ingroup_clade,
                            outgroup_clade, feat_lengths, chr)

    if(debug) {
        acc_scoring_paths <- list(
            'len_25' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/score_raw/25/chr17_score_25.csv",
            'len_50' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/score_raw/50/chr17_score_50.csv"
        )
    }

    print(acc_scoring_paths)

    print('processing step 4 / 6')

    filter_acc_scoring_path <- filter_acc_scoring(acc_scoring_paths,
                                    alignment_id, ingroup_clade, feat_lengths,
                                    chr)

    if(debug) {

        filter_acc_scoring_path <- list(
            'len_25' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/custom_filtering/mammals/acceleration/score_filtered/25/chr17_score_25_filtered_norm.csv",
            'len_50' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/custom_filtering/mammals/acceleration/score_filtered/50/chr17_score_50_filtered_norm.csv"
            )
    }

    print(filter_acc_scoring_path)

    print('processing step 5 / 6')

    cons_clade <- paste0(ingroup_clade, '_', outgroup_clade)
    score_cons_raw_file_paths <- calculate_cons_scoring(alignment_id,
                        cons_clade, feat_lengths, chr, acc_scoring_paths)

    if(debug) {

        score_cons_raw_file_paths <- list(
            'len_25' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/custom_filtering/mammals_sarcopterygii/candidate_elements/score_raw/25/chr17_score_25.csv",
            'len_50' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/custom_filtering/mammals_sarcopterygii/candidate_elements/score_raw/50/chr17_score_50.csv"
        )
    }

    print(score_cons_raw_file_paths)

    print('processing step 6 / 6')

    score_cons_file_paths <- filter_cons_scoring(score_cons_raw_file_paths,
                                alignment_id, cons_clade, feat_lengths, chr)

    if(debug) {
        score_cons_file_paths <- list(
            'len_25' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/custom_filtering/mammals_sarcopterygii/candidate_elements/score_filtered/25/chr17_score_25.csv",
            'len_50' = "/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/custom_filtering/mammals_sarcopterygii/candidate_elements/score_filtered/50/chr17_score_50.csv"
        )
    }

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

    if(clade == 'mammals') {

        ingroup_clade <- 'mammals'
        outgroup_clade <- 'sarcopterygii'


    } else if(clade == 'aves') {

        ingroup_clade <- 'aves'
        outgroup_clade <- 'sarcopterygii'

    }

    lapply(chrs, function(chr) {
        process_filtering_chr(alignment_id, chr, ingroup_clade, outgroup_clade)
    })



}

#-----------------------------------------
#TODO: analizar todo por cromosoma, al final hacer los joins de todos los archivos, arreglar los id de phastcons, que quede un unico archivo final
#-----------------------------------------



# main <- function() {
#
#     # preparation_aves()
#     # preparation_mammals()
#
#
#
#     # consenso:
#
#     # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/77_way/output/acceleration/aves/fasta_alignments/25/chr24'
#     #
#     # calculate_consensus_sequences(alignment_id = '77_way', clade = 'aves',
#     #                               fasta_base_path = fasta_base_path,
#     #                               chr = 24, feat_length = 25)
#     #
#     #
#     # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/77_way/output/acceleration/sarcopterygii/fasta_alignments/25/chr24'
#     #
#     # calculate_consensus_sequences(alignment_id = '77_way', clade = 'sarcopterygii',
#     #                               fasta_base_path = fasta_base_path,
#     #                               chr = 24, feat_length = 25)
#     #
#     #
#     # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/77_way/output/acceleration/aves/fasta_alignments/50/chr24'
#     #
#     # calculate_consensus_sequences(alignment_id = '77_way', clade = 'aves',
#     #                               fasta_base_path = fasta_base_path,
#     #                               chr = 24, feat_length = 50)
#     #
#     # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/77_way/output/acceleration/sarcopterygii/fasta_alignments/50/chr24'
#     #
#     # calculate_consensus_sequences(alignment_id = '77_way', clade = 'sarcopterygii',
#     #                               fasta_base_path = fasta_base_path,
#     #                               chr = 24, feat_length = 50)
#
#
# #---------------------------
#
#     # mammals
#
#     # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/fasta_alignments/25/chr22'
#     #
#     # calculate_consensus_sequences(alignment_id = '100_way', clade = 'mammals',
#     #                               fasta_base_path = fasta_base_path,
#     #                               chr = 22, feat_length = 25)
#     #
#     #
#     # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/sarcopterygii/fasta_alignments/25/chr22'
#     #
#     # calculate_consensus_sequences(alignment_id = '100_way', clade = 'sarcopterygii',
#     #                               fasta_base_path = fasta_base_path,
#     #                               chr = 22, feat_length = 25)
#     #
#     #
#     # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/fasta_alignments/50/chr22'
#     #
#     # calculate_consensus_sequences(alignment_id = '100_way', clade = 'mammals',
#     #                               fasta_base_path = fasta_base_path,
#     #                               chr = 22, feat_length = 50)
#     #
#     # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/sarcopterygii/fasta_alignments/50/chr22'
#     #
#     # calculate_consensus_sequences(alignment_id = '100_way', clade = 'sarcopterygii',
#     #                               fasta_base_path = fasta_base_path,
#     #                               chr = 22, feat_length = 50)
#
#     # alignment_id <- '100_way'
#     # ingroup_clade <- 'mammals'
#     # outgroup_clade <- 'sarcopterygii'
#     # feat_length <- 25
#     # chr <- 22
#     # acc_scoring_path <- calculate_acc_scoring_data(alignment_id, ingroup_clade,
#     #                     outgroup_clade, feat_length, chr)
#
#     # filter_acc_scoring_path <- filter_acc_scoring_data(acc_scoring_path,
#     #                                 alignment_id, ingroup_clade, feat_length,
#     #                                 chr)
#
#     # # cons scoring data:
#     # score_cons_raw_file_path <- calculate_cons_scoring_data(
#     #                             alignment_id = alignment_id,
#     #                             cons_clade = 'mammals_sarcopterygii',
#     #                             feat_length = feat_length,
#     #                             chr = chr,
#     #                             acc_scoring_path = filter_acc_scoring_path)
#
#
#     # score_cons_file_path <- filter_cons_scoring_data(
#     #     score_cons_raw_file_path, alignment_id, 'mammals_sarcopterygii', feat_length,
#     #     chr)
#
#     #-----------------------------------------
#     #TODO: analizar todo por cromosoma, al final hacer los joins de todos los archivos, arreglar los id de phastcons, que quede un unico archivo final
#     #-----------------------------------------
#
#     # feat <- rphast::feat(seqname = 'galGal6', start = 1077205, end = 1077229)
#     #
#     # tmp1 <- alignment[[1]]
#     # for (i in seq(1, length(tmp1))) {
#     #     tmp_i <- tmp1[i]
#     #     if(str_detect(tmp_i, pattern = "\\-")) {
#     #         print(i)
#     #         # print(tmp_i)
#     #     }
#     # }
#
#     # 77_way: https://genome.ucsc.edu/cgi-bin/hgTables?db=galGal6&hgta_group=compGeno&hgta_track=cons77way&hgta_table=phastConsElements77way&hgta_doSchema=describe+table+schema
#     # If the space is insufficient and the gap size is a multiple of 3, a "*" is displayed; other gap sizes are indicated by "+"
#
#     # Aviso:
#     #     In .Call2("fasta_index", filexp_list, nrec, skip, seek.first.rec,  :
#     #                   reading FASTA file /u01/home/pbeati/.local/share/R/cladeAcc/77_way/output/acceleration/aves/fasta_alignments/25/chr24/acc_25_chr24_1077205_1077229.fasta: ignored 104 invalid one-letter sequence codes
#
#
#
#     # scoring
#     # alignment_id <- '100_way'
#     # ingroup_clade <- 'mammals'
#     # outgroup_clade <- 'sarcopterygii'
#     # feat_length <- 25
#     # chr <- 22
#     #
#     # calculate_acc_scoring_data(alignment_id, ingroup_clade, outgroup_clade,
#     #                                    feat_length, chr)
#
# }



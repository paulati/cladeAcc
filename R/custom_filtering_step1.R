
# step 1
extract_elements_alignments_fasta <- function(alignment_id, chr, ingroup_clade,
                                              outgroup_clade,
                                              split_lengths = c(25, 50)) {

    base_path <- accelerated_elements_path(alignment_id, ingroup_clade)

    elements_alignments_paths <- lapply(split_lengths, function(split_length) {

        file_name <- accelerated_non_parametric_phyloP_file_name(chr,
                                                                 split_length)

        elements_file_path <- file.path(base_path$non_parametric_phyloP,
                                        split_length,
                                        file_name)

        # ingroup
        ingroup_alignments_paths <- extract_fasta_alignments(elements_file_path,
                                                             alignment_id, chr,
                                                             ingroup_clade)

        # outgroup
        outgroup_alignemnts_paths <- extract_fasta_alignments(elements_file_path,
                                                              alignment_id, chr,
                                                              outgroup_clade)

    })

    names(elements_alignments_paths) <- paste0('len_', split_lengths)


    return(elements_alignments_paths)

}


extract_fasta_alignments <- function(elements_file_path, alignment_id, chr,
                                     clade) {

    data <- load_non_parametric_phyloP_elements(file_path = elements_file_path)
    data_acc <- data$data
    mask <- data_acc$non_parametric_FDR < 0.05

    # sum(mask)
    # nrow(data_acc)

    data_acc_significative <- data_acc[mask, ]

    config <- load_config()
    filtering_config <- config$custom_filtering$filtering[[alignment_id]]
    filtering_clade_config <- filtering_config[[clade]]
    sequence_names <- filtering_clade_config$species

    #alignment <- load_multiz_alignment(alignment_id, chr, sequence_names)
    alignment <- NA #delay load of alignment

    feat_seq_name <- sequence_names[1] # first sequence is reference sequence

    feats <- rphast::feat(seqname = feat_seq_name, # sin esto no funciona extract.feature.msa
                          start = data_acc_significative$start,
                          end = data_acc_significative$end,
                          feature = data_acc_significative$feature,
                          score =  data_acc_significative$score)

    # checking existence of output directory:
    dummy_start <- feats[1, 'start']
    dummy_end <- feats[1, 'end']
    dummy_chr <- chr
    feat_length <- dummy_end - dummy_start + 1
    dummy_file_path <- fasta_alignment_paths(alignment_id, clade, feat_length,
                                             dummy_chr, dummy_start, dummy_end)
    base_path <- dirname(dummy_file_path$tmp)
    if(! dir.exists(base_path)) {
        dir.create(base_path, recursive =  TRUE, showWarnings = FALSE)
    }

    #alignment_paths <- apply(feats, MARGIN = 1, FUN = function(feat) {
    alignment_paths <- unlist(lapply(seq(1, nrow(feats)), function(x) {

        feat <- feats[x, ]
        start <- as.integer(feat$start)
        end <- as.integer(feat$end)
        feat_length <- end - start + 1

        out_file_path <- fasta_alignment_paths(alignment_id, clade, feat_length,
                                               chr, start, end)

        if(! file.exists(out_file_path$tmp)) {

            if(is.na(alignment)) {
                alignment <- load_multiz_alignment(alignment_id, chr, sequence_names)
            }

            feat_alignment <- rphast::extract.feature.msa(rphast::copy.msa(alignment),
                                                          feat, do4d = FALSE,
                                                          pointer.only = FALSE)

            # replace all * by N (DNA_ALPHABET
            # "A" "C" "G" "T" "M" "R" "W" "S" "Y" "K" "V" "H" "D" "B" "N" "-" "+" ".")
            # in order to rad fasta format later

            replace_asterix_by_N <- stringr::str_replace_all(feat_alignment$seqs,
                                                    pattern = "\\*",
                                                    replacement = "N")

            feat_alignment$seqs <- replace_asterix_by_N

            rphast::write.msa(feat_alignment,
                              file = out_file_path$tmp,
                              format = "FASTA",
                              pretty.print = FALSE)

            rm()
            gc()

        }

        result <- out_file_path$tmp

        return(result)
    }))


    return(alignment_paths)

}


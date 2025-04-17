


#step 2
calculate_consensus_sequences <- function(alignment_id, clade,
                                          fasta_base_path, chr, feat_length) {

    out_file_path <- fasta_alignment_consensus_sequence_paths(alignment_id,
                                                    clade, feat_length, chr)

    if(! file.exists(out_file_path$tmp)) {

        config <- load_config()
        filtering_config <- config$custom_filtering$filtering[[alignment_id]]
        filtering_clade_config <- filtering_config[[clade]]
        threshold <- as.numeric(filtering_clade_config$consensus_threshold)

        fasta_files_paths <- list.files(fasta_base_path, full.names = TRUE)

        consensus_sequences <- unlist(lapply(fasta_files_paths,
            function(fasta_files_path) {

                x <- Biostrings::readDNAStringSet(file = fasta_files_path,
                                                   format = "fasta")
                consensus_seq <- Biostrings::consensusString(x,
                                                          threshold = threshold,
                                                          ambiguityMap = 'X')
                return(consensus_seq)
            }))

        result <- data.frame('path' = fasta_files_paths,
                             'consensus_sequence' = consensus_sequences)

        out_file_path <- save_elements_consensus_info(result, alignment_id,
                                                      clade, feat_length, chr)

    }

    return(out_file_path)
}



calculate_elements_consensus_sequences <- function(alignment_id, chr,
        ingroup_clade, outgroup_clade, feat_lengths = c(25, 50)) {


    result <- lapply(feat_lengths, function(feat_length) {

        fasta_base_path_in <- fasta_alignment_base_paths(alignment_id,
                                            ingroup_clade, feat_length, chr)

        # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/fasta_alignments/25/chr22'

        consensus_in_file_path <- calculate_consensus_sequences(alignment_id,
                                    ingroup_clade, fasta_base_path_in$tmp,
                                    chr, feat_length)


        fasta_base_path_out <- fasta_alignment_base_paths(alignment_id,
                                                         outgroup_clade,
                                                         feat_length, chr)

        # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/fasta_alignments/25/chr22'

        consensus_out_file_path <- calculate_consensus_sequences(alignment_id, outgroup_clade,
                                      fasta_base_path_out$tmp, chr, feat_length)

        result <- list('consensus_ingroup_file_path' = consensus_in_file_path,
                       'consensus_outgroup_file_path' = consensus_out_file_path )

        return(result)
    })

    names(result) <- paste0('len_', feat_lengths)

    return(result)


    # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/sarcopterygii/fasta_alignments/25/chr22'
    #
    # calculate_consensus_sequences(alignment_id = '100_way', clade = 'sarcopterygii',
    #                               fasta_base_path = fasta_base_path,
    #                               chr = 22, feat_length = 25)
    #
    #
    # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/mammals/fasta_alignments/50/chr22'
    #
    # calculate_consensus_sequences(alignment_id = '100_way', clade = 'mammals',
    #                               fasta_base_path = fasta_base_path,
    #                               chr = 22, feat_length = 50)
    #
    # fasta_base_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/output/acceleration/sarcopterygii/fasta_alignments/50/chr22'
    #
    # calculate_consensus_sequences(alignment_id = '100_way', clade = 'sarcopterygii',
    #                               fasta_base_path = fasta_base_path,
    #                               chr = 22, feat_length = 50)



}

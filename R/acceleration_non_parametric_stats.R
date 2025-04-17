# all these methods are internal

# clade mammals_sarcopterygii or aves_sarcopterygii
get_informative_elements <- function(alignment_id, clade, chr) {

    candidate_regions <- load_conserved_elements(alignment_id,
                                                 clade, chr,
                                                 in_common = TRUE)

    return(candidate_regions)

}

simulate_elements <- function(alignment, informative_elements, feat_length, nrep) {

    informative_elements_seqname <- unique(informative_elements$seqname)

    align_copy <- rphast::copy.msa(alignment)

    # lo necesito para que coincidan los nombres de secuencias en extract.feature.msa
    ref_specie_name <- rphast::names.msa(align_copy)[1] # reference is first specie
    informative_elements$seqname <- ref_specie_name

    informative_elements_alignment <- rphast::extract.feature.msa(align_copy,
                                         informative_elements,
                                         pointer.only = TRUE)

    simulated_msa <- rphast::sample.msa(x = informative_elements_alignment,
                                  size = nrep * feat_length,
                                  replace = TRUE)

    # produce features allowing long alignment to be interpreted as
    # concatenation of shorter alignments
    start_indexes <- seq(from = 1, by = feat_length, length.out = nrep)

    #length(startIdx)
    #max(startIdx)

    features <- rphast::feat(
                    seqname = paste0("sim_", informative_elements_seqname),
                     src = "sim",
                     feat=".",
                     start = start_indexes,
                     end = start_indexes + feat_length - 1)

    result <- list(msa = simulated_msa, features = features)

    return(result)


    # --------------------------

    # informative_elements <- get_informative_elements(alignment_id, clade, chr)
    #
    # alignment_copy <- copy.msa(alignment)
    #
    # element_align <- extract.feature.msa(alignment_copy, informative_elements,
    #                                      pointer.only = TRUE)


}



calculate_non_parametric_distribution_phyloP <- function(alignment_id, alignment,
                                               neutral_model_labeled,
                                               acc_clade, branch_label,
                                               informative_elements_clade,
                                               chr, feat_length, nrep) {

    # distribution_data_file_name <- paste0('chr', chr)

    non_parametric_distribution <-
        load_non_parametric_phyloP_distribution(alignment_id = alignment_id,
                                                clade = acc_clade, chr = chr,
                                                feat_length = feat_length)

    if(nrow(non_parametric_distribution$data) == 0) {

        informative_elements <- get_informative_elements(alignment_id,
                                            informative_elements_clade, chr)

        simulated_elements <- simulate_elements(alignment, informative_elements,
                                            feat_length, nrep)


        non_parametric_distribution_data <- rphast::phyloP(
                                mod = neutral_model_labeled,
                                msa = simulated_elements$msa,
                                mode="ACC",
                                features = simulated_elements$features,
                                branches = branch_label)

        non_parametric_distribution_file_path <-
            save_non_parametric_phyloP_distribution(
                non_parametric_distribution_data, alignment_id,
                acc_clade, feat_length, chr)


    } else {

        non_parametric_distribution_file_path <-
            non_parametric_distribution$path

    }

    return(non_parametric_distribution_file_path)

}

# informative_elements_clade <- candidate_regions_clade
# called from acceleration.R
calculate_non_parametric_stats <- function(alignment_id, alignment,
                                           neutral_model_labeled,
                                           acc_clade, branch_label,
                                           informative_elements_clade,
                                           chr, feat_lengths, nrep) {


    paths <- lapply(feat_lengths, function(x) {
        distribution_data_path <- calculate_non_parametric_distribution_phyloP(
            alignment_id, alignment, neutral_model_labeled, acc_clade, branch_label,
            informative_elements_clade, chr, x, nrep)
        return(distribution_data_path)
    })

    names(paths) <- paste0('len_', feat_lengths)

    return(paths)


}

# clade <- acc_clade
# called from acceleration.R

calculate_non_parametric_significance <- function(obs_phyloP_result_paths,
                                                  non_parametric_stats_paths,
                                                  alignment_id, clade, chr) {

    obs_phyloP_elements <- lapply(obs_phyloP_result_paths, function(x){
        obs_phyloP_elements_len <- load_obs_phyloP_elements(
            file_path = x)
        return(obs_phyloP_elements_len)
    })
    names(obs_phyloP_elements) <- names(obs_phyloP_result_paths)

    non_parametric_phyloP_distribution <- lapply(non_parametric_stats_paths,
        function(x){
            non_parametric_phyloP_distribution_len <-
                load_non_parametric_phyloP_distribution(file_path = x)
            return(non_parametric_phyloP_distribution_len)
        })
    names(non_parametric_phyloP_distribution) <-
        names(non_parametric_stats_paths)


    non_parametric_phyloP_elements_paths <- lapply(names(obs_phyloP_elements),
                                                   function(len) {

        non_parametric_pval <- unlist(sapply(obs_phyloP_elements[[len]]$lnlratio,
            function(x) {
                distribution_lnlratio <-
                    non_parametric_phyloP_distribution[[len]]$data$lnlratio
                pval <-
                    sum(x <= distribution_lnlratio) / length(distribution_lnlratio)
                return(pval)
            }))

        non_parametric_pval_FDR <- p.adjust(non_parametric_pval, method = "BH")

        obs_phyloP_elements[[len]]$non_parametric_FDR <- non_parametric_pval_FDR

        # output_file_name <- paste0('chr', chr)

        feat_length <- as.integer(stringr::str_replace(len,
                                              pattern = 'len_',
                                              replacement = ''))

        non_parametric_phyloP_elements_path <- save_non_parametric_phyloP_elements(
            obs_phyloP_elements[[len]],
            alignment_id, clade,
            feat_length, chr)


        return(non_parametric_phyloP_elements_path)
    })

    names(non_parametric_phyloP_elements_paths) <- names(obs_phyloP_elements)


    return(non_parametric_phyloP_elements_paths)

}



test <- function() {

    alignment_id <- '100_way'
    acc_clade <- 'mammals'
    acc_branch <- 'MAMMALS'
    informative_elements_clade <- 'mammals_sarcopterygii'
    chr <- 22
    feat_length <- 25
    nrep <- 100000


    config <- load_config()
    filtering_config <- config$acceleration$filtering[[alignment_id]]
    sequence_names <- filtering_config[[acc_clade]]
    alignment <- load_multiz_alignment(alignment_id, chr, sequence_names)

    label_neutral_model_func <- label_neutral_model_tree_mammals
    neutral_model <- load_neutral_model(alignment_id, step = "acc")
    neutral_model_labeled <- label_neutral_model_func(neutral_model)


    # ejemplo: https://github.com/paulati/acc_regions_mammals_aves/blob/223eb53cf3f202baaf1323cc8718c2210e99dc32/mammals/acceleration/nonparasim_acceleration_mammals.R

    non_parametric_distribution_data <-
        non_parametric_distribution_phyloP(alignment_id, alignment,
                                           neutral_model_labeled,
                                           acc_clade, acc_branch,
                                           informative_elements_clade,
                                           chr, feat_length, nrep)


    # #esto esta ordenado por score descendente desde acceleration_mammals.R
    # setwd(obsphyloP.output.base.path)
    # obsPhyloP.file.name <- paste("chr", chr.id, "_obsPhyloP_", split.length, ".csv", sep="")
    # obsPhyloP <- read.table(obsPhyloP.file.name, sep="\t", header=TRUE)
    #
    # #nrow(obsPhyloP)
    # #max(obsPhyloP$end)
    #
    # nonParaPval <- sapply(obsPhyloP$lnlratio, empirical.pval, nonParaPhyloP$lnlratio)
    # nonParaFDR <- p.adjust(nonParaPval, method="BH")
    #
    # #length(nonParaFDR)
    # #nonParaFDR[5000:10000]
    #
    # #esto esta ordenado por posicion (start end), no coincide con el orden de obsPhyloP
    # setwd(split.feats.base.path)
    # file.name <- paste("chr", chr.id, "_", split.feats.file.name, "_", split.length, ".csv", sep="")
    # split.feats <- read.table(file.name, sep="\t", header = TRUE)
    # #ordeno split.feats en el orden de obsPhyloP:
    # order.df <- obsPhyloP[, c("chr", "start", "end")]
    # order.df$sort_id <- c(1:nrow(order.df))
    # tmp <- merge(split.feats, order.df, by.x = c("seqname", "start", "end"), by.y = c("chr", "start", "end"))
    # split.feats.sort <- tmp[order(tmp$sort_id), ]
    #
    # #tiene el orden de obsPhyloP
    # indexes <- (nonParaFDR < 0.05)
    # nonParaSigFeats <- split.feats.sort[indexes,]
    #
    # nonParaSigFeats$feature <- "mammalsAR"
    # nonParaSigFeats$score <- obsPhyloP$lnlratio[indexes]
    # nonParaSigFeats$seqname <- paste("hg38.chr", chr.id, sep="")
    #
    # #max(nonParaSigFeats$end)
    #
    # #ordeno por score:
    # nonParaSigFeats.sort <- nonParaSigFeats[order(- nonParaSigFeats$score), ]
    # file.name <- paste("chr", chr.id, "_mar_", split.length, ".gff", sep="")
    # setwd(nonparasimphyloP.gff.output.base.path)
    # #no usar write feat porque redondea el score a 0
    # #write.feat(nonParaSigFeats.sort, file.name)
    # write.table(nonParaSigFeats.sort, file.name, sep="\t",
    #             col.names = FALSE, row.names = FALSE, quote = FALSE)
    #
    # q <- nrow(nonParaSigFeats.sort)
    # result.bed <- data.frame(chr=character(q),
    #                          start=integer(q),
    #                          end=character(q))
    #
    # chr.name <- paste("chr", chr.id, sep="")
    # result.bed$chr <- rep(chr.name, q)
    # result.bed$start <- nonParaSigFeats.sort$start
    # result.bed$end <- nonParaSigFeats.sort$end
    #
    # setwd(nonparasimphyloP.bed.output.base.path)
    # file.name <- paste("chr", chr.id, "_acc_", split.length, ".bed", sep="")
    # write.table(result.bed, file.name, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



}





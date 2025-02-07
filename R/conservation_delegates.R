required_species_features_mammals_100way <- function(align)
{

    # Primate subset: al menos hg38 y otoGar3. (No sé si incorporar a otro primate, como Rhesus)
    # Euarchontoglires subset: al menos mm10, oryCun2, ochPri3
    # Laurasiatheria subset: al menos susScr3, turTru2, bosTau8, felCat8, myoLuc2
    # Afrotheria subset: al menos loxAfr3 o echTel2
    # Mammal subset: al menos dasNov3 o monDom5 o macEug2, ornAna1.

    cons_req_primates <- c("hg38", "otoGar3" )
    cons_req_euarchontoglires <- c("mm10", "oryCun2", "ochPri3")
    cons_req_laurasiatheria <- c("susScr3", "turTru2", "bosTau8", "felCat8", "myoLuc2")
    cons_req_mammal <- c("ornAna1")

    cons_req <- c(cons_req_primates, cons_req_euarchontoglires, cons_req_laurasiatheria, cons_req_mammal)

    #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
    req_regions_1 <- rphast::informative.regions.msa(align, min.numspec=11, spec=cons_req,
                                             refseq="hg38", gaps.inf=FALSE)

    cons_opt_afrotheria <- c("loxAfr3", "echTel2")
    req_regions_2 <- rphast::informative.regions.msa(align, min.numspec=1, spec=cons_opt_afrotheria,
                                             refseq="hg38", gaps.inf=FALSE)

    cons_opt_mammal <- c("dasNov3", "monDom5", "macEug2")
    req_regions_3 <- rphast::informative.regions.msa(align, min.numspec=1, spec=cons_opt_mammal,
                                             refseq="hg38", gaps.inf=FALSE)

    #coverage.feat: Any features object passed into this function which is stored as a pointer
    #to an object stored in C may be reordered (sorted) by this function.
    intersection_regions <- rphast::coverage.feat(req_regions_1,
                                          req_regions_2,
                                          req_regions_3,
                                          or=FALSE, get.feats=TRUE)

    intersection_regions_order <- sort(intersection_regions, decreasing = FALSE)

    return(intersection_regions_order)

}

required_species_features_sarcopterygii_100way <- function(align)
{
    cons_req_species <- c("allMis1", "anoCar2")
    cons_req_turtles <- c("cheMyd1", "chrPic2", "pelSin1", "apaSpi1")
    #esta no esta en el alineamiento

    #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
    req_regions <- rphast::informative.regions.msa(align, min.numspec=2,
                                           spec = cons_req_species,
                                           refseq = "hg38", gaps.inf = FALSE)

    #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
    turtle_regions <- rphast::informative.regions.msa(align, min.numspec=1,
                        spec=cons_req_turtles, refseq="hg38", gaps.inf = FALSE)

    #calculo la interseccion entre req_regions y turtle_regions
    #coverage.feat: Any features object passed into this function which is stored as a pointer
    #to an object stored in C may be reordered (sorted) by this function.
    intersection_regions <- rphast::coverage.feat(req_regions, turtle_regions,
                                                  or = FALSE, get.feats = TRUE)

    intersection_regions_order <- sort(intersection_regions, decreasing = FALSE)

    return(intersection_regions_order)

}

required_species_features_aves_77way <- function(align) {

    alignment_species <- rphast::names.msa(align)

    falcon_eagle <- c()
    if ("falChe1" %in% alignment_species) {
        falcon_eagle <- c(falcon_eagle, "falChe1")
    }
    if("falPer1" %in% alignment_species) {
        falcon_eagle <- c(falcon_eagle, "falPer1")
    }
    if("halLeu1" %in% alignment_species) {
        falcon_eagle <- c(falcon_eagle, "halLeu1")
    }

    zebrafinch_medium_ground_finch <- c()
    if("taeGut2" %in% alignment_species) {
        zebrafinch_medium_ground_finch <- c(zebrafinch_medium_ground_finch, "taeGut2")
    }
    if("geoFor1" %in% alignment_species) {
        zebrafinch_medium_ground_finch <- c(zebrafinch_medium_ground_finch, "geoFor1")
    }

    ostrich_tinamou <- c()
    if("strCam1" %in% alignment_species) {
        ostrich_tinamou <- c(ostrich_tinamou, "strCam1")
    }
    if("tinGut2" %in% alignment_species) {
        ostrich_tinamou <- c(ostrich_tinamou, "tinGut2")
    }

    #required 11 birds

    # #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
    req_regions_1 <- rphast::informative.regions.msa(align, min.numspec = 1,
                                            spec = falcon_eagle,
                                            refseq = "galGal6",
                                            gaps.inf = FALSE)

    req_regions_2 <- rphast::informative.regions.msa(align, min.numspec = 1,
                                        spec = zebrafinch_medium_ground_finch,
                                        refseq = "galGal6",
                                        gaps.inf = FALSE)

    req_regions_3 <- rphast::informative.regions.msa(align, min.numspec = 1,
                                            spec = ostrich_tinamou,
                                            refseq = "galGal6",
                                            gaps.inf = FALSE)

    req_regions_4 <- rphast::informative.regions.msa(align, min.numspec = 11,
                                             spec = NULL,
                                            refseq = "galGal6",
                                            gaps.inf = FALSE)
    # spec: the default value of NULL implies use all species in the alignment.
    #req_regions_2 <- informative.regions.msa(align.2, min.numspec=10, spec=NULL,
    #                                       refseq="oryzias_latipes", gaps.inf=FALSE)

    #coverage.feat: Any features object passed into this function which is stored as a pointer
    #to an object stored in C may be reordered (sorted) by this function.
    # result <- coverage.feat(req_regions, get.feats=TRUE)

    intersection_regions <- rphast::coverage.feat(req_regions_1,
                                          req_regions_2,
                                          req_regions_3,
                                          req_regions_4,
                                          or=FALSE, get.feats=TRUE)


    intersection_regions_order <- sort(intersection_regions, decreasing = FALSE)

    return(intersection_regions_order)

}


required_species_features_sarcopterygii_77way <- function(align) {

    alignment_species <- rphast::names.msa(align)

    cons_req_species <- c()
    if ("allMis1" %in% alignment_species) {
        cons_req_species <- c(cons_req_species, "allMis1")
    }
    if ("anoCar2" %in% alignment_species) {
        cons_req_species <- c(cons_req_species, "anoCar2")
    }

    cons_req_turtles <- c()
    if ("cheMyd1" %in% alignment_species) {
        cons_req_turtles <- c(cons_req_turtles, "cheMyd1")
    }
    if ("chrPic2" %in% alignment_species) {
        cons_req_turtles <- c(cons_req_turtles, "chrPic2")
    }
    if ("pelSin1" %in% alignment_species) {
        cons_req_turtles <- c(cons_req_turtles, "pelSin1")
    }
    if ("apaSpi1" %in% alignment_species) {
        cons_req_turtles <- c(cons_req_turtles, "apaSpi1")
    }

    #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
    req_regions <- rphast::informative.regions.msa(align, min.numspec = 2,
                                            spec = cons_req_species,
                                            refseq = "galGal6",
                                            gaps.inf = FALSE)

    #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
    turtle_regions <- rphast::informative.regions.msa(align, min.numspec = 1,
                                                spec = cons_req_turtles,
                                                refseq = "galGal6",
                                                gaps.inf = FALSE)

    #calculo la interseccion entre req_regions y turtle_regions
    #coverage.feat: Any features object passed into this function which is stored as a pointer
    #to an object stored in C may be reordered (sorted) by this function.
    intersection_regions <- rphast::coverage.feat(req_regions, turtle_regions,
                                                  or = FALSE, get.feats = TRUE)

    intersection_regions.order <- sort(intersection_regions, decreasing = FALSE)

    return(intersection_regions.order)

}

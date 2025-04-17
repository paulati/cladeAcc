# private function

label_neutral_model_tree <- function(neutral_model, clade, label) {

    original_tree <- neutral_model$tree

    decimal_group_pattern <- "(\\d+\\.?\\d*)"

    if(clade == 'mammals') {

        # original_tree <- ":0.21065):0.0768943,ornAna1:0.501668):0.124913,((((("

        pattern_head <- paste0("(?<head>.*ornAna1\\:",
                               decimal_group_pattern,
                               "\\))")

    } else if(clade == 'aves') {

        # original_tree <- "tinGut2:0.162739):0.0403754):0.186966,("
        # "tinGut2:0.162739):0.0403754)AVES:0.186966,("

        pattern_head <- paste0("(?<head>.*tinGut2\\:",
                               decimal_group_pattern,
                               "\\)\\:",
                               decimal_group_pattern,
                               "\\))")
    } else {

        print("error: clade should be mammals or aves")

    }

    head <- stringr::str_extract(original_tree, pattern_head)
    tail <- stringr::str_sub(original_tree, start = nchar(head)+1,
                    end = nchar(original_tree))

    labeled_tree <- paste0(head, label, tail)

    neutral_model$tree <- labeled_tree

    return(neutral_model)

}

# ------------------------------------------------------------
# public functions

#' Delegate function for ....
#' @param neutral_model description
#' @export
label_neutral_model_tree_mammals <- function(neutral_model) {
    clade <- 'mammals'
    label <- 'MAMMALS'
    result <- label_neutral_model_tree(neutral_model, clade, label)
    return(result)
}

#' Delegate function for ....
#' @param neutral_model description
#' @export
label_neutral_model_tree_aves <- function(neutral_model) {
    clade <- 'aves'
    label <- 'AVES'
    result <- label_neutral_model_tree(neutral_model, clade, label)
    return(result)
}


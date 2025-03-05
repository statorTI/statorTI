#' @title
#' Pairwise Distance; Jaccard
#'
#' @description
#' Retruns a distance for the jaccard distance of two sets
#'
#' @param cat1 vector binary list of cells in/not in cat1
#' @param cat2 vector binary list of cells in/not in cat2
#'
#' @importFrom magrittr "%>%"
#'
#' @export
#' @importFrom magrittr "%>%"
#'
#' @export

jaccard_distance <- function(cat1, cat2) {
    df <- data.frame(cat1 = cat1, cat2 = cat2) %>%
        dplyr::mutate(both = cat1 * cat2)

    dist <- 1 - sum(df$both) / (sum(df$cat1) + sum(df$cat2))

    return(round(dist, 3))
}


#' @title
#' Pairwise Distance; Hypergeometric
#'
#' @description
#' Retruns a p-value for a one-sided hypergeometric test on two subsets
#'
#' @param cat1 vector binary list of cells in/not in cat1
#' @param cat2 vector binary list of cells in/not in cat2
#'
#' @importFrom magrittr "%>%"
#'
#' @export

hypergeo_test <- function(cat1, cat2) {
    if (identical(cat1, cat2)) {
        return(1)
    }

    df <- data.frame(cat1 = cat1, cat2 = cat2) %>%
        dplyr::mutate(both = cat1 * cat2)

    p_val <- stats::phyper(
        sum(df$both) - 1,
        sum(df$cat1),
        nrow(df) - sum(df$cat1),
        sum(df$cat2),
        lower.tail = FALSE
    )

    return(p_val)
}


#' Pairwise Distance; Chi Square
#'
#' @description
#' Retruns a p-value for a chi square test of independence on two subsets
#'
#' @param cat1 vector binary list of cells in/not in cat1
#' @param cat2 vector binary list of cells in/not in cat2
#'
#' @details
#' Not recommended for use, as it tests both postive and negative
#' correlation. For TI we usually only want to test positive correlation.
#'
#' @importFrom magrittr "%>%"
#'
#' @export

chisqr_test <- function(cat1, cat2) {
    warning("This function is not reccomended for use, see details")
    return(stats::chisq.test(cat1, cat2, correct = TRUE)$p.value)
}


#' @importFrom magrittr "%>%"

pairwise_distance <- function(df, fn) {
    CellID <- NULL

    df <- df %>%
        dplyr::transmute(dplyover::across2x(!CellID, !CellID,
            .fns = fn,
            .comb = "all"
        )) %>%
        unique()

    return(df)
}

#' @importFrom magrittr "%>%"

create_graph <- function(clusters, gene_cell_matrix, fn) {
    values <- pairwise_distance(fn = fn, df = clusters)
    dim <- sqrt(length(values))

    m <- matrix(values, dim, dim)
    names <- colnames(clusters[, -1])
    colnames(m) <- names
    rownames(m) <- names

    graph <- igraph::graph_from_adjacency_matrix(m,
        weighted = TRUE,
        mode = "max"
    ) %>%
        igraph::simplify(remove.multiple = FALSE)

    return(graph)
}


holm_bonferroni <- function(x, alpha) {
    sorted_pvals <- sort(x)
    correction <- alpha / rev(seq_along(x))
    index <- match(FALSE, (sorted_pvals <= correction))
    if (is.integer(index) && length(index) == 0L) {
        return(1)
    }

    max_pval <- sorted_pvals[[index]]

    return(max_pval)
}


simplify_graph <- function(graph, method, alpha, bound) {

    if (method == "mst") {
        graph <- igraph::mst(graph)
    } else if (method == "cyclic") {
        if (is.double(alpha) && is.double(bound)) {
            stop("Only one of `alpha` or `bound` should be set")
        } else if (is.double(alpha)) {
            pval <- holm_bonferroni(x = igraph::E(graph)$weight, alpha = alpha)

            graph <- igraph::delete_edges(
                graph, igraph::E(graph)[igraph::E(graph)$weight >= pval]
            )
        } else if (is.double(bound)) {
            graph <- igraph::delete_edges(
                graph, igraph::E(graph)[igraph::E(graph)$weight >= bound]
            )
        }
    } else if (method == "full") {
        return(graph)
    }

    return(graph)
}

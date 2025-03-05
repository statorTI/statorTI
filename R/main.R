#' statorTI
#'
#' @description
#' Perform Stator based trajectory inference
#'
#' @inheritParams cluster_func
#' @param gene_cell_matrix matrix a gene cell matrix, not currently used
#' @param fn string the function to use to calculate the distance between
#' clusters on the graph
#' @inheritParams igraph_plotter
#' @param alpha float the significance level for p-testing, should not be set
#' if bound is set
#' @param bound float an upper bound for pruning, should not be set if alpha
#' is set
#'
#' @details
#'
#' By default fn = jaccard_distance, a built in function. Other build in
#' options include hypergeo_test and chisqr_test. The later is not recommended.
#' Can also be a user defined function which takes numerical vectors as
#' input and returns a scalar.
#'
#' Use bound to set the upper bound for pruning if using a metric based
#' function. use alpha if using a statistical test based function. Both of
#' these values are only used if graph_type is "cyclic".
#'
#' @export

statorTI <- function(state_list,
                     cell_list,
                     gene_cell_matrix,
                     fn = jaccard_distance,
                     graph_type = "full",
                     edge_labels = FALSE,
                     alpha = FALSE,
                     bound = FALSE,
                     root = "") {
    clusters <- cluster_func(cell_list = cell_list, state_list = state_list)

    graph <- create_graph(
        clusters = clusters, gene_cell_matrix = gene_cell_matrix,
        fn = fn
    )

    graph <- simplify_graph(
        graph = graph, method = graph_type, alpha = alpha,
        bound = bound
    )

    igraph_plotter(graph,
        graph_type = graph_type, edge_labels = edge_labels,
        root = root
    )

    return(graph)
}

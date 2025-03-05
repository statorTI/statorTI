#' @title
#' Plot statorTI Graphs
#'
#' @description
#' Plot statorTI graphs using the reccomended settings.
#'
#' @param graph igraph object a statorTI graph to be plotted
#' @param graph_type string one of "mst", "cyclic", "full".
#' @param edge_labels bool include edge weights on the graph? Default FALSE
#' @param root string Only used is graph_type = "mst", set the root vertex
#'
#' @importFrom magrittr "%>%"
#'
#' @export

igraph_plotter <- function(graph, graph_type, edge_labels, root) {
    if (grepl("mst", graph_type, fixed = TRUE)) {
        if (root != "") {
            graph_layout <- igraph::layout_as_tree(graph, root = root)
        } else {
            graph_layout <- igraph::layout_as_tree
        }
    } else if (grepl("cyclic", graph_type, fixed = TRUE)) {
        graph_layout <- igraph::layout_with_graphopt(
            graph,
            charge = 100,
            mass = 20,
            spring.constant = 0.5,
            spring.length = 0.5,
            max.sa.movement = 20
        )
    } else if (grepl("full", graph_type, fixed = TRUE)) {
        graph_layout <- igraph::layout_in_circle
    }

    if (edge_labels) {
        edge_labels <- round(igraph::E(graph)$weight, 5)
    } else {
        edge_labels <- NULL
    }

    igraph::plot.igraph(graph,
        layout = graph_layout,
        vertex.size = 40, vertex.size2 = 10, vertex.shape = "rectangle",
        vertex.color = "#00BFC4", vertex.label.color = "black",
        vertex.frame.color = "#00BFB4", vertex.label.family = "sans",
        edge.color = "#F8766D", edge.label = edge_labels
    )
}


#' @title
#' Clusters upset plot
#'
#' @description
#' create an upset plot showing number of overlapping cluster combinations.
#'
#' @inheritParams cluster_func
#' @param n integer the number of interactions to plot
#'
#' @returns A ggplot object.
#'
#' @importFrom magrittr "%>%"
#'
#' @export

clusters_upset <- function(cell_list, state_list, n = Inf) {
    CellID <- State <- count <- NULL

    data <- data.frame(CellID = cell_list, State = state_list)

    clusters <- data %>%
        dplyr::group_by(CellID) %>%
        dplyr::mutate(dplyr::across(State, ~ paste0(., collapse = " "))) %>%
        dplyr::mutate(dplyr::across(State, ~ strsplit(., split = " "))) %>%
        unique() %>%
        ggplot2::ggplot(ggplot2::aes(x = as.list(State))) +
        ggplot2::geom_bar() +
        ggplot2::geom_text(
            stat = "count",
            ggplot2::aes(label = ggplot2::after_stat(count)),
            vjust = -1
        ) +
        ggupset::scale_x_upset(n_intersections = n) +
        ggplot2::xlab("State") +
        ggplot2::theme_classic()

    return(clusters)
}


#' @title
#' Bubble plot for state overlap
#'
#' @description
#' Bubble proportionality plot for pairwise state variables.
#'
#' @inheritParams cluster_func
#' @param state1 string
#' @param state2 string
#'
#' @returns a shiny object of the bubble plot
#'
#' @importFrom magrittr "%>%"
#'
#' @export

bubble_func <- function(cell_list, state_list, state1, state2) {
    State <- CellID <- NULL

    data <- data.frame(CellID = cell_list, State = state_list)

    clusters <- data %>%
        stats::filter(State == state1 | State == state2) %>%
        dplyr::group_by(CellID) %>%
        dplyr::mutate(counts = dplyr::n()) %>%
        dplyr::mutate(State = dplyr::case_when(
            counts == 2 ~ "both",
            .default = State
        )) %>%
        dplyr::mutate(State = as.factor(State)) %>%
        unique()

    state_count <- summary(clusters$State)
    state_percents <- as.character(round(state_count / sum(state_count) * 100))
    state_labels <- paste(
        as.character(sort(unique(clusters$State))), ":",
        state_percents, "%"
    )

    plot <- bubbles::bubbles(
        value = state_count,
        label = state_labels
    )

    return(plot)
}


#' @title
#' Skree Plot for PCA
#'
#' @description
#' Create a skree plot for pca.
#'
#' @param n integer number of pca components to plot
#' @param sdev numerical vector of standard deviations
#'
#' @returns a ggplot object displaying a skree plot
#'
#' @importFrom magrittr "%>%"
#'
#' @export

scree_plot <- function(n, sdev) {
    variance <- sdev^2 / sum(sdev^2)

    plot <- ggplot2::ggplot(data = NULL, ggplot2::aes(
        x = c(1:n),
        y = variance[1:n]
    )) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 4) +
        ggplot2::xlab("Principal Component") +
        ggplot2::ylab("Variance Explained") +
        ggplot2::ggtitle("Scree Plot") +
        ggplot2::ylim(0, max(variance[1:n]) + 0.25)

    return(plot)
}

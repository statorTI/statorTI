#' @title
#' One-Hot-Encode Cluster Annotations
#'
#' @description
#' One-hot-encode a metadata file for use with statorTI. In order for statorTI
#' to work effectively a large number of cells should belong to multiple
#' clusters, as opposed to in traditional one-hot-encoding.
#'
#' @param cell_list vector list of cells, same length as state_list
#' @param state_list vector list of state, same length as state_list
#'
#' @returns dataframe A C x S dataframe of cells with one-hot-encoded states
#'
#' @importFrom magrittr "%>%"
#'
#' @export

cluster_func <- function(cell_list, state_list) {

    State <- value <- State <- NULL

    data <- tibble::tibble(State = state_list, CellID = cell_list) %>%
        dplyr::mutate(value = 1)  %>%
        tidyr::spread(State, value, fill = 0, sep = ":")

    return(data)
}

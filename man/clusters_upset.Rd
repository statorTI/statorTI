% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphing.R
\name{clusters_upset}
\alias{clusters_upset}
\title{Clusters upset plot}
\usage{
clusters_upset(cell_list, state_list, n = Inf)
}
\arguments{
\item{cell_list}{vector list of cells, same length as state_list}

\item{state_list}{vector list of state, same length as state_list}

\item{n}{integer the number of interactions to plot}
}
\value{
A ggplot object.
}
\description{
create an upset plot showing number of overlapping cluster combinations.
}

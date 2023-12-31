% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MDSplot.R
\name{MDSplot}
\alias{MDSplot}
\title{Multidimensional scaling plot}
\usage{
MDSplot(dist, titlename = NULL, type = NULL, nc = NULL)
}
\arguments{
\item{dist}{An object of class "dist"}

\item{titlename}{An title name of multidimensional scaling (MDS) plot. If it is \code{NULL},
the title of ggplot will not be shown in the output.}

\item{type}{MDS methods implemented. The default values is metric MDS \code{type = "mMDS"}.
Users can also choose classical MDS \code{type = "cMDS"} and non-metric MDS \code{type = "nMDS"}.}

\item{nc}{Number of colours. If it is \code{NULL}, we will use 3 gradient colours for scatter plot.}
}
\value{
A ggplot output
}
\description{
Multidimensional scaling plot
}

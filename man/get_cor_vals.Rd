% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Correlation.R
\name{get_cor_vals}
\alias{get_cor_vals}
\title{Calculate correlation between all columns of matrices X and Y}
\usage{
get_cor_vals(
  X,
  Y,
  fdr_filter = 0.1,
  n_cores = 1,
  self_only = FALSE,
  method = "pearson",
  pbar = FALSE
)
}
\arguments{
\item{X}{matrix}

\item{Y}{matrix}

\item{fdr_filter}{FDR filter for correlation values}

\item{self_only}{only calculate variables with same name in X and Y. example use would be}

\item{method}{'method' argument to cor and cor.test
correlated DNA methylation for a gene and expression of same gene.}

\item{pbar}{print progress bar}
}
\description{
Calculate correlation between all columns of matrices X and Y
}

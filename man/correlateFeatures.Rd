% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Correlation.R
\name{correlateFeatures}
\alias{correlateFeatures}
\title{Correlate Features pulled from SummExpDR object}
\usage{
correlateFeatures(
  x,
  features_x,
  features_y,
  prune_x = NULL,
  prune_y = NULL,
  assayKey = NULL,
  ...
)
}
\arguments{
\item{x}{SummExpDR object}

\item{features_x}{features to pull that will be represented in X in get_cor_vals}

\item{features_y}{features to pull that will be represented in Y in get_cor_vals}

\item{prune_x}{regex to remove from values in features_x.
use when different data modalities have measurements for same gene and you want to set self_only to TRUE}

\item{prune_y}{regex to remove from values in features_y.
use when different data modalities have measurements for same gene}

\item{assayKey}{assay to pull if using assay data}

\item{...}{other args to get_cor_vals}
}
\description{
Correlate Features pulled from SummExpDR object
}

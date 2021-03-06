% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DimReduction.R
\name{mPCA}
\alias{mPCA}
\title{Run PCA on PCA (multistage PCA')}
\usage{
mPCA(x, assay_use = "imputed_mat", min_var = 0.05, max_dims = 5)
}
\arguments{
\item{x}{multiExp object}

\item{assay_use}{assay to pull from}

\item{min_var}{= minimum pct variance explained per PC in PCA on each dataset}

\item{max_dims}{maximum number of dimensions to pull out}
}
\description{
Intended as a solution to capturing major sources of signal in multiple datasets while preventing
signal from being diminished due to presence of more features in data.
}

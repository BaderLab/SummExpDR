% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiExpIntegration.R
\name{createMultiExp}
\alias{createMultiExp}
\title{create multiExp object}
\usage{
createMultiExp(summ_exp_list, assays_use = NULL, ...)
}
\arguments{
\item{summ_exp_list}{list of SummarizedExperiment objects}

\item{assays_use}{character vector or numeric vector with names corresponding to names in summ_exp_list}

\item{...}{args to pass to KNN_impute()}
}
\description{
Creates a new SummarizedExperiment object for summ_exp slot, containing combined data from multiple experiments
E.g. RNA, DNA methylation, CNV, etc. Data is standard normalized prior to performing imputation. Choice of
word 'Experiment' to describe multiple data modalities is to avoid confusion with concept of 'assays' in
SummarizedExperiment, which assumes identical set of biological feature names with different measurements
}

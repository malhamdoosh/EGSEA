#' @name arraydata
#' @title Example dataset for \code{egsea.ma}
#'
#' @description This is dataset is provided as an example only, for the \code{egsea.ma} function.
#'
#' @docType data
#' @format A Named List containing two elements, \code{arrays} and \code{targets}
#' 
#' \code{arrays} is a \code{limma::EList} with 4 values: source, E, genes and other.
#' 
#' \code{targets} is a data frame with 12 rows and 6 variables: Array, SampleID, Condition, Chip, Section and Experiment.
#'
#' @source URL
"arraydata"


loadArrayData <- function() {
    data(arraydata)
    if (is.null(arraydata)) message("arraydata failed to load")
}
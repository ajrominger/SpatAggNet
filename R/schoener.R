#' @title Schoener similarity
#'
#' @description Calculate Schoener similarity metric
#'
#' @details Faster implementation than what is available in
#'     \code{spaa::niche.overlap}.
#'
#' @param x site by species matrix (sites as rows, species as columns,
#'     abundances in cells)
#'
#' @return an \code{ncol(x)} by \code{ncol(x)} matrix of species-species
#'     similarities
#'
#' @export

schoener <- function(x) {
    y <- t(x)
    y <- y / rowSums(y)
    x <- t(y)

    o <- lapply(1:ncol(x), function(i) {
        colSums(abs(x[, i] - x)) / 2
    })

    o <- matrix(1 - unlist(o), nrow = ncol(x))
    diag(o) <- 0

    return(o)
}

#' @title Z-squared test statistic of the negative binomial distribution
#'
#' @description Exact test of the negative binomial
#'
#' @details Compares the observed log likelihood of the nagative binomial distribution
#' given the data to the sampling distribution of log likelihoods given the negative
#' binomial distribution is true.  This test statistic is Chi-square distributed with
#' d.f. = 1, meaning value larger than \code{qchisq(0.95, 1)} rejects the "null"
#' hypothesis that the data come from a negative binomial.
#'
#' @param x a numeric vector
#' @param size the size parameter
#' @param mu the mean parameter
#'
#' @return The z-squared statistic (a numeric of length 1)
#'
#' @export

nb_zsq <- function(x, size, mu) {
    ll_obs <- sum(dnbinom(x, size = size, mu = mu, log = TRUE))
    n <- length(x)

    # sampling distribution of log probabilities
    p0 <- p_naught(size, mu)

    # sampling mean and var
    m <- sum(p0 * exp(p0)) * n
    v <- sum((m/n - p0)^2 * exp(p0)) * n

    # z^2-value
    z <- ((ll_obs - m) / sqrt(v))^2

    return(as.numeric(z))
}

# helper function to compute sampling distribution of log probabilities
p_naught <- function(size, mu, mod = 'nbinom') {
    n0 <- 0:10^5

    if(mod == 'nbinom') {
        p0 <- dnbinom(n0, size = size, mu = mu, log = TRUE)
    } else {
        p0 <- dpois(n0, lambda = mu, log = TRUE)
    }

    p0 <- p0[is.finite(p0)]

    if(exp(p0[length(p0)]) > .Machine$double.eps^0.75) {
        n0add <- (n0[length(p0)] + 1):10^6
        p0add <- dnbinom(n0add, size = size, mu = mu, log = TRUE)

        if(mod == 'nbinom') {
            p0add <- dnbinom(n0add, size = size, mu = mu, log = TRUE)
        } else {
            p0add <- dpois(n0add, lambda = mu, log = TRUE)
        }

        p0add <- p0add[is.finite(p0add)]

        p0 <- c(p0, p0add)
    }

    return(p0)
}

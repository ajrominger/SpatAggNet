#' @title Fit the negative binomial
#'
#' @description Fits the negative binomial SSAD to a site x species matrix
#'
#' @details If likelihood optimization fails for a given species, \code{NA}
#'     will be returned for all row elements corresponding to that species.
#'
#' @param x a site by species matrix
#'
#' @return A \code{data.frame} with columns giving the species IDs, the two
#'     parameters of the negative binomial distribution, the log likelihood
#'     of the negative binomial distribution, and the log likelihood of the
#'     Poisson distribution.  Each row is one species and the species ID refers
#'     to the column index of that species in the site by species matrix
#'     \code{x}.
#'
#' @export

nb_fit <- function(x) {
    stats <- sapply(1:ncol(x), function(j) {
        # list of fitted params and max log lik
        fit <- nb_mle(x[, j])

        # max log lik for poisson
        pois_loglik <- sum(dpois(x[, j], mean(x[, j]), log = TRUE))

        # z-squared goodness of fit statistic for nb
        # if estimated k is not finite, replace it with 0.001
        # for z-sq calculation purposes
        k2use <- ifelse(is.finite(fit$estimate[1]),
                        fit$estimate[1],
                        0.001)

        z2 <- nb_zsq(x[, j], k2use, fit$estimate[2])

        # results to return
        f <- c(fit$estimate,
               nb_loglik = fit$loglik,
               pois_loglik = pois_loglik,
               delta_AIC = 2 * (2 - fit$loglik - (1 - pois_loglik)),
               zsq = z2,
               abund = sum(x[, j]),
               nocc = sum(x[, j] > 0))

        return(f)
    })

    stats <- rbind(sp_ID = 1:ncol(x), stats) |>
        t() |>
        as.data.frame()

    return(stats)
}




#' @title Maximum likelihood estimation of negative binomial distribution
#'
#' @description Finds the maximum likelihood estimates of the dispersion (k)
#'     and mean (mu) parameters of the negative binomial distribution.
#'
#' @details This function profiles the log likelihood assuming
#'     \code{mu = mean(x)} and uses \code{optimize} to find \code{k} that
#'     maximizes the log likelihood function built from \code{dnbinom}.
#'
#'     Maximizing the profile likelihood is faster and avoids computational
#'     issues when \code{var(x)} is close to \var{mean(x)}
#'
#' @param x an integer vector
#'
#' @return A \code{lsit} with components \code{estiamte} giving the MLE and
#'     \code{loglik} giving the value of the log likelihood function at the
#'     maximum. If \code{var(x) <= mean(x)} then the MLE does not exist and the
#'     function returns \code{size = Inf}.
#'
#' @export

nb_mle <- function(x) {
    x <- x[!is.na(x)]

    fit <- optimize(nb_ll, c(.Machine$double.eps, 1000), maximum = TRUE,
                    mu = mean(x), dat = x)

    names(fit) <- c('estimate', 'loglik')

    if(fit$estimate >= 990) fit$estimate <- Inf

    fit$estimate <- c(size = fit$estimate, mu = mean(x))

    return(fit)
}

# log likelihood function for nb distribution
nb_ll <- function(k, mu, dat) {
    sum(dnbinom(dat, size = k, mu = mu, log = TRUE))
}


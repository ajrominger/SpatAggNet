#' @title Simulate positive and negative interactions from purely random data
#'
#' @description Simulate spatial replicates of abundance data and apply the same
#' analytical pipeline to those data that would be applied to real data
#'
#' @details \code{simPlusMinus} draws random SAD and SSAD shapes from the raw data
#' and uses these to simulate more data and calculate network statistics on those
#' simulated data. \code{simpleSim} assumes one SAD and one SSAD and simulates data from
#' those, again calculated network statistics.
#'
#' Note: any value passed to \code{ssad_type} other than \code{'nbinom'} results
#' in a Poisson SSAD (i.e., there are only two options, negative binomial specified by
#' \code{'nbinom'} or Poisson specified by anything else)
#'
#' @param sad_stats a \code{data.frame} with columns \code{mod}, \code{par1}, \code{par2}
#' @param nsite number of sites to simulate
#' @param nspp number of species to simulate
#' @param mc_cores number of cores to use in \code{parallel::mclapply}
#' @param ssad_type string specifying SSAD shape (e.g. \code{'nbinom'})
#' @param sadfun function to generate random SAD sample
#' @param ssadfun function to generate random SSAD sample
#' @param k_fun function to relate k parameter of the SSAD to abundance
#' @param nsim number of simulations to run
#' @param dist_fun character name of the distance function to use (default \code{'schoener'})
#'
#' @return a \code{data.frame} with \code{<= nsim} rows (some simulations may be
#' thrown out if they do not meet data filtering requirements), and columns corresponding
#' to summary statistics about the positive and negative network characteristics
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export
#' @rdname simPlusMinus

simPlusMinus <- function(sad_stats, mc_cores,
                         ssad_type = 'nbinom', k_fun, nsim,
                         dist_fun = 'schoener') {
    # indexes for SAD data and nsite, nspp data
    idx_SAD <- sample(nrow(sad_stats), nsim, replace = TRUE)

    # loop over replicates and make communities, then calculate networks
    # o <- parallel::mclapply(1:nsim, mc.cores = mc_cores, FUN = function(i) {
    o <- lapply(1:nsim, FUN = function(i) {
        j <- idx_SAD[i]
        nsite <- sad_stats$nsite[j]
        nspp <- sad_stats$nspp[j]

        # make SAD
        rfun <- get(paste0('r', sad_stats$mod[j]))
        pars <- as.numeric(sad_stats[j, 2:3])
        pars <- pars[!is.na(pars)]
        abund <- do.call(rfun, c(list(nspp), as.list(pars)))

        # mean per site abund for each spp
        mu <- abund / nsite

        # make site by spp matrix by sampling ssad
        if(ssad_type == 'nbinom') {
            # calculate k (size param)
            k <- k_fun(abund)

            mat <- make_mat(nsite, nspp, rnbinom, size = k, mu = mu)
        } else {
            mat <- make_mat(nsite, nspp, rpois, lambda = mu)
        }

        return(.simCleanup(mat, dist_fun))
    })

    return(.outCleanup(o))
}



#' @export
#' @rdname simPlusMinus

simpleSim <- function(nsite, nspp, mc_cores, sadfun, ssadfun, nsim, dist_fun = 'schoener') {
    # make SAD
    ii <- rep(1:nsim, each = nspp)
    X <- sadfun(nspp * nsim)
    X <- split(X, ii)

    o <- parallel::mclapply(X, mc.cores = mc_cores, FUN = function(abund) {
        # calculate known quantities from abund
        J <- sum(abund)
        mu <- abund / nsite

        # loop over abundances and generate ssad
        mat <- make_mat(nsite, nspp, ssadfun, mu = mu)

        return(.simCleanup(mat, dist_fun))
    })

    return(.outCleanup(o))
}


make_mat <- function(nsite, nspp, rfun, ...) {
    # the way `rfun` recycles params, we need to fill matrix
    # by rows
    matrix(rfun(nsite * nspp, ...), nrow = nsite, byrow = TRUE)
}


.simCleanup <- function(mat, dist_fun) {
    defaultNames <- c('all.v',
                      'pos.n', 'pos.v', 'pos.rho.rho', 'pos.p', 'pos.m', 'pos.wm',
                      'neg.n', 'neg.v', 'neg.rho.rho', 'neg.p', 'neg.m', 'neg.wm')
    mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]

    if(any(dim(mat) < 10)) {
        o <- rep(NA, length(defaultNames))
    } else {
        o <- unlist(plusMinus(mat, dist_fun = dist_fun))
    }

    names(o) <- defaultNames

    return(o)
}

.outCleanup <- function(o) {
    o <- as.data.frame(do.call(rbind, o))
    o[is.na(o$pos.rho.rho) | is.na(o$neg.rho.rho), ] <- NA
    o <- o[!is.na(o$pos.n), ]

    return(o)
}

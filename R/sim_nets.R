#' @title Simulate site by species matrices and compute networks from them
#'
#' @description Networks are built from simulated data that mirrors real
#'     SADs and SSADs but has no real interactions/associations. Then summary
#'     statistics are calculated for those networks
#'
#' @details \code{sim_nets} draws random SAD and SSAD shapes from real data
#'     and uses these to simulate more data and calculate network statistics on
#'     those simulated data. \code{simple_sim_nets} assumes one SAD and one S
#'     SAD and simulates data from those, again calculated network statistics.
#'
#'     Note: any value passed to \code{ssad_type} other than \code{"nbinom"}
#'     results in a Poisson SSAD (i.e., there are only two options, negative
#'     binomial specified by \code{'nbinom'} or Poisson specified by anything
#'     else)
#'
#' @param sad_stats a \code{data.frame} with columns \code{mod}, \code{par1},
#'     \code{par2}
#' @param nsite number of sites to simulate
#' @param nspp number of species to simulate
#' @param mc_cores number of cores to use in \code{parallel::mclapply}
#' @param ssad_type string specifying SSAD shape (e.g. \code{"nbinom"})
#' @param sadfun function to generate random SAD sample
#' @param ssadfun function to generate random SSAD sample
#' @param k_fun function to relate k parameter of the SSAD to relative abundance
#' @param nsim number of simulations to run
#' @param alpha significance cut off for comparing observed species similarities
#'     to null similarities, defaults to \code{0.05}
#' @param dist_fun character name of the distance function to use for
#'     calculating species similarities (default \code{"schoener"})
#'
#' @return a \code{data.frame} with \code{<= nsim} rows (some simulations may be
#'     thrown out if they do not meet data filtering requirements), and columns
#'     corresponding to summary statistics about network characteristics
#'
#' @export
#' @rdname sim_nets

sim_nets <- function(sad_stats,
                     mc_cores,
                     ssad_type = "nbinom",
                     k_fun,
                     nsim,
                     alpha = 0.05,
                     B = 999,
                     dist_fun = "schoener") {

    # indexes for randomly sampled SADs
    idx_SAD <- sample(nrow(sad_stats), nsim, replace = TRUE)

    # loop over replicates and make site by spp matrix,
    # then calculate networks and network stats
    # o <- parallel::mclapply(1:nsim, mc.cores = mc_cores, FUN = function(i) {
    o <- lapply(1:nsim, FUN = function(i) {
        j <- idx_SAD[i]
        nsite <- sad_stats$nsite[j]
        nspp <- sad_stats$nspp[j]

        # get function to randomly sample SAD
        rfun <- paste0("r", sad_stats$mod[j]) |>
            getExportedValue("pika", name = _)

        # get and clean pars to pass to rfun
        pars <- as.numeric(sad_stats[j, c("par1", "par2")])
        pars <- pars[!is.na(pars)]

        # make SAD
        abund <- do.call(rfun, c(list(nspp), as.list(pars)))

        # mean per site abund for each spp will be used for
        # sampling from ssad
        mu <- abund / nsite

        # make site by spp matrix by sampling ssad
        if(ssad_type == 'nbinom') {
            # calculate k (size param)
            k <- k_fun(abund)

            mat <- make_mat(nsite, nspp, rnbinom, size = k, mu = mu)
        } else {
            mat <- make_mat(nsite, nspp, rpois, lambda = mu)
        }

        res <- run_net_stats(mat, alpha = alpha, B = B,
                             dist_fun = dist_fun)
        return(res)
    })

    # simplify output to a data.frame
    o <- do.call(rbind, o) |>
        as.data.frame()

    return(o)
}






make_mat <- function(nsite, nspp, rfun, ...) {
    # the way `rfun` recycles params, we need to fill matrix
    # by rows
    matrix(rfun(nsite * nspp, ...), nrow = nsite, byrow = TRUE)
}


# function to run `net_stats` with some safety checks
run_net_stats <- function(mat, ...) {

    # names of output
    pn_net_names <- c("num_v",
                      "num_e",
                      "rho",
                      "m_abund",
                      "wm_abund")

    pn_net_names <- paste(rep(c("pos", "neg"),
                              each = length(pn_net_names)),
                          pn_net_names,
                          sep = "_")

    out_names <- c("connectance",
                   "cluster_coef",
                   "apl",
                   "modularity",
                   "H_degree",
                   "mean_degree",
                   pn_net_names)

    # remove empty rows and cols
    mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]

    # if, after removing empty rows/cols, matrix is too small
    # (i.e. dim[i] < 10) just return NAs
    if(any(dim(mat) < 10)) {
        o <- rep(NA, length(out_names))
        names(o) <- out_names
    } else {
        o <- net_stats(mat, ...)
    }

    return(o)
}






# @export
# @rdname simPlusMinus

# simpleSim <- function(nsite, nspp, mc_cores, sadfun, ssadfun, nsim, dist_fun = 'schoener') {
#     # make SAD
#     ii <- rep(1:nsim, each = nspp)
#     X <- sadfun(nspp * nsim)
#     X <- split(X, ii)
#
#     o <- parallel::mclapply(X, mc.cores = mc_cores, FUN = function(abund) {
#         # calculate known quantities from abund
#         J <- sum(abund)
#         mu <- abund / nsite
#
#         # loop over abundances and generate ssad
#         mat <- make_mat(nsite, nspp, ssadfun, mu = mu)
#
#         return(.simCleanup(mat, dist_fun))
#     })
#
#     return(.outCleanup(o))
# }

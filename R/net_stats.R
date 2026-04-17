#' @title Make a network and calculate statistics
#'
#' @description Make a network based on a distance-like function and compute
#'     network summary statistics
#'
#' @param x site by species matrix (sites as rows, species as columns,
#'     abundances in cells)
#' @param alpha the significance cut off level
#' @param B number of replicates for the null model permutations
#' @param dist_fun character name of the distance function to use
#'     (default \code{"schoener"}). Note that the function must
#'     calculate distances between columns, not rows
#'
#' @return an \code{ncol(x)} by \code{ncol(x)} matrix of species-species similarities
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

net_stats(obsdat[[1]])

net_stats <- function(x, alpha = 0.05, B = 999, dist_fun = "schoener") {
    browser()
    dfun <- get(dist_fun)
    dfun <- function(x) vegan::vegdist(t(x)) |> as.matrix()

    # metacommunity abundance
    metaX <- colSums(x)

    # observed distance
    obs_d <- dfun(x)
    obs_d <- obs_d[lower.tri(obs_d)]

    # null distances
    null_mats <- stats::simulate(vegan::nullmodel(x, "r2dtable"), nsim = B)
    null_d <- lapply(1:dim(null_mats)[3], function(i) {
        d <- dfun(null_mats[, , i])
        return(d[lower.tri(d)])
    })

    # add observed to nulls and make the whole thing
    # a matrix where columns represent replicates
    null_d <- cbind(obs_d, do.call(cbind, null_d))

    # probabilities that the observed species-species distances
    # are >= or <= the null
    ppos <- rowMeans(null_d >= null_d[, 1])
    pneg <- rowMeans(null_d <= null_d[, 1])

    # significantly positive and negative edges
    epos <- as.integer(ppos <= alpha)
    eneg <- as.integer(pneg <= alpha)

    # number of significantly positive and negative edges
    npos <- sum(epos)
    nneg <- sum(eneg)

    # edge list
    elist <- cbind(t(combn(1:ncol(x), 2)), epos, eneg)

    # igraph stats
    all_net_stats <- g_stats(elist)

    # centrality v. abundance corr and means
    corMeanPos <- plus_minus_stats(elist[elist[, 3] == 1, 1:2, drop = FALSE],
                                   metaX)
    corMeanNeg <- plus_minus_stats(elist[elist[, 4] == 1, 1:2, drop = FALSE],
                                   metaX)

    # number of species in each network
    vpos <- length(unique(as.vector(elist[elist[, 3] == 1, 1:2])))
    vneg <- length(unique(as.vector(elist[elist[, 4] == 1, 1:2])))

    return(list(all = c(v = ncol(x)),
                pos = c(n = npos, corMeanPos),
                neg = c(n = nneg, corMeanNeg)))
}


# function to calculate statistics relevant to sub networks
# of positive or negative interactions:
#    - number of vertices
#    - corr between abundance and centrality
#    - mean abundances, unweighted and weighted by centrality

# and mean abundances, unweighted and weighted by centrality
plus_minus_stats <- function(elist, abund) {
    if(nrow(elist) < 3) {
        return(c(v = NA, rho = NA, p = NA, m = NA, wm = NA))
    } else {
        # centrality
        cen <- table(c(elist[, 1], elist[, 2]))

        # abundance ordered same as centrality vector
        x <- abund[as.integer(names(cen))]

        # correlation output
        o <- cor(x, cen, method = 'spearman')

        # relative abundance
        relx <- x / sum(abund)

        return(c(
            # number vertices
            num_v = length(cen),

            # centrality-abund corr
            rho = o,

            # abundance means
            m = mean(relx),
            wm = weighted.mean(relx, cen)
        ))
    }
}


g_stats <- function(elist) {
    # subset edge list to only connected vertices
    elist <- elist[elist[, 3] > 0 | elist[, 4] > 0, 1:2, drop = FALSE]

    # make sure vertices not interpreted as indexes
    elist[, 1] <- as.character(elist[, 1])
    elist[, 2] <- as.character(elist[, 2])

    if(nrow(elist) < 3) {
        return(c(
            connectance = NA,
            cluster_coef = NA,
            apl = NA,
            modularity = NA,
            H_degree = NA,
            mean_degree = NA
        ))
    } else {

        g <- igraph::graph_from_edgelist(
            elist,
            directed = FALSE
        )

        dd <- igraph::degree_distribution(g)

        return(c(
            connectance = igraph::edge_density(g),
            cluster_coef = igraph::transitivity(g, type = "global"),
            apl = igraph::mean_distance(g, directed = FALSE),
            modularity = igraph::cluster_louvain(g) |>
                igraph::membership() |>
                igraph::modularity(g, membership = _),
            H_degree = vegan::diversity(dd),
            mean_degree = sum((1:length(dd)) * dd)
        ))

    }
}

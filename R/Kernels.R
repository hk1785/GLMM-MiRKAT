Kernels <-
function (otu.tab, tree, distance = c("Jaccard", "BC", "U.UniFrac", "G.UniFrac.0.5", "W.UniFrac")) {
	pa.otu.tab <- otu.tab
	for (i in 1:nrow(otu.tab)) {
		ind <- which(otu.tab[i,] > 0)
		pa.otu.tab[i,ind] <- 1
	}
	jac <- as.matrix(bcdist(pa.otu.tab))
    bc <- as.matrix(bcdist(otu.tab))
    unifs <- GUniFrac(otu.tab, tree, alpha = c(0.5, 1))$unifracs
    u.unif <- unifs[, , "d_UW"]
    g.unif.0.5 <- unifs[, , paste("d_", 0.5, sep = "")]
    w.unif <- unifs[, , "d_1"]
	for (j in 1:nrow(otu.tab)) {
        ind <- is.na(u.unif[j, ])
        if (sum(ind) != 0) u.unif[, ind] <- 0
        ind <- is.na(g.unif.0.5[j, ])
        if (sum(ind) != 0) g.unif.0.5[, ind] <- 0
        ind <- is.na(w.unif[j, ])
        if (sum(ind) != 0) w.unif[, ind] <- 0
    }
    jac.k <- D2K(jac)
    bc.k <- D2K(bc)
    u.unif.k <- D2K(u.unif)
    g.unif.0.5.k <- D2K(g.unif.0.5)
    w.unif.k <- D2K(w.unif)
    Ks <- list(Jaccard = jac.k, BC = bc.k, U.UniFrac = u.unif.k, G.UniFrac.0.5 = g.unif.0.5.k, W.UniFrac = w.unif.k)
    return(Ks[!is.na(match(names(Ks), distance))])
}

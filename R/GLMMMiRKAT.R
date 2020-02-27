GLMMMiRKAT <-
function (y, cov = NULL, id, time.pt = NULL, Ks, model, slope = FALSE, n.perm = 5000) {
    if (is.null(time.pt) & slope) {
        stop("time.pt is required for the random slope model")
    }
    if (is.null(time.pt)) {
        id <- as.character(id)
        ind <- order(id)
        id <- id[ind]
        y <- y[ind]
		if (!is.null(cov)) {
			if (is.matrix(cov) | is.data.frame(cov)) {
				cov <- as.data.frame(cov)[ind, ]
			} else {
				cov <- as.data.frame(cov[ind])
			}
		}
        for (i in 1:length(Ks)) {
            Ks[[i]] <- Ks[[i]][ind, ind]
        }
    } else {
        ind <- order(as.character(time.pt))
        id <- as.character(id)
        id <- id[ind]
        y <- y[ind]
        if (!is.null(cov)) {
			if (is.matrix(cov) | is.data.frame(cov)) {
				cov <- as.data.frame(cov)[ind, ]
			} else {
				cov <- as.data.frame(cov[ind])
			}
		}
        time.pt <- time.pt[ind]
        for (i in 1:length(Ks)) {
            Ks[[i]] <- Ks[[i]][ind, ind]
        }
        ind <- order(id)
        id <- id[ind]
        y <- y[ind]
		if (!is.null(cov)) {
			if (is.matrix(cov) | is.data.frame(cov)) {
				cov <- as.data.frame(cov)[ind, ]
			} else {
				cov <- as.data.frame(cov[ind])
			}
		}
        time.pt <- time.pt[ind]
        for (i in 1:length(Ks)) {
            Ks[[i]] <- Ks[[i]][ind, ind]
        }
    }
    if (model == "gaussian") {
        y <- as.numeric(scale(y))
        if (is.null(cov)) {
            if (!is.null(time.pt) & slope) {
                fit <- lmer(y ~ (time.pt | id))
            } else {
                fit <- lmer(y ~ (1 | id))
            }
        } else {
			cov <- as.data.frame(cov)
			if (ncol(cov) == 1) {
				colnames(cov) <- "cov1"
				if (!is.null(time.pt) & slope) {
					formula.0 <- as.formula("y ~ (time.pt | id) + cov1")
					fit <- lmer(formula.0, data = cov)
				} else {
					formula.0 <- as.formula("y ~ (1 | id) + cov1")
					fit <- lmer(formula.0, data = cov)
				}
			} else {
				colnames(cov) <- paste("cov", 1:ncol(cov), sep = "")
				if (!is.null(time.pt) & slope) {
					formula.0 <- as.formula(paste(paste(c("y ~ (time.pt | id)", 
					paste("cov", 1:ncol(cov), sep = "")), collapse = " + "), 
					collapse = " + "))
					fit <- lmer(formula.0, data = cov)
				} else {
					formula.0 <- as.formula(paste(paste(c("y ~ (1 | id)", 
					paste("cov", 1:ncol(cov), sep = "")), collapse = " + "), 
					collapse = " + "))
					fit <- lmer(formula.0, data = cov)
				}
			}
        }
    } else {
        if (is.null(cov)) {
            if (!is.null(time.pt) & slope) {
                fit <- glmer(y ~ (time.pt | id), family = model)
            } else {
                fit <- glmer(y ~ (1 | id), family = model)
            }
        } else {
			cov <- as.data.frame(cov)
			if (ncol(cov) == 1) {
				colnames(cov) <- "cov1"
				if (!is.null(time.pt) & slope) {
					formula.0 <- as.formula("y ~ (time.pt | id) + cov1")
					fit <- glmer(formula.0, family = model, data = cov)
				} else {
					formula.0 <- as.formula("y ~ (1 | id) + cov1")
					fit <- glmer(formula.0, family = model, data = cov)
				}
			} else {
				colnames(cov) <- paste("cov", 1:ncol(cov), sep = "")
				if (!is.null(time.pt) & slope) {
					formula.0 <- as.formula(paste(paste(c("y ~ (time.pt | id)", 
					paste("cov", 1:ncol(cov), sep = "")), collapse = " + "), 
					collapse = " + "))
					fit <- glmer(formula.0, family = model, data = cov)
				} else {
					formula.0 <- as.formula(paste(paste(c("y ~ (1 | id)", 
					paste("cov", 1:ncol(cov), sep = "")), collapse = " + "), 
					collapse = " + "))
					fit <- glmer(formula.0, family = model, data = cov)
				}
			}
        }
    }
    fe <- as.numeric(getME(fit, "X") %*% getME(fit, "fixef"))
    re <- as.numeric(getME(fit, "Z") %*% getME(fit, "b"))
    if (model == "gaussian") {
        mu <- fe + re
        inv.de.mu <- diag(length(y))
    }
    if (model == "binomial") {
        mu <- exp(fe + re)/(1 + exp(fe + re))
        inv.de.mu <- diag(as.numeric(mu * (1 - mu)))
    }
    if (model == "poisson") {
        mu <- exp(fe + re)
        inv.de.mu <- diag(as.numeric(mu))
    }
    y.star <- as.numeric(fe + re + ginv(inv.de.mu) %*% (y - mu))
    r <- y.star - fe
    r.cov <- as.matrix((getME(fit, "Z") %*% Matrix::crossprod(getME(fit, 
        "Lambdat")) %*% getME(fit, "Zt")))
    if (model == "gaussian") {
        e.var <- diag(rep(var(mu), length(y)))
    }
    if (model == "binomial") {
        e.var <- diag(rep(1, length(y)))
    }
    if (model == "poisson") {
        e.var <- diag(rep(1, length(y)))
    }
    v.cov <- r.cov + e.var
    inv.vcov <- ginv(v.cov)
    r.s <- list()
    if (!slope) {
        id.time.pt <- id
        p.ind <- as.matrix(getPermuteMatrix(n.perm, length(r), 
            strata = id.time.pt))
        clust.sizes <- as.numeric(names(table(table(id.time.pt))))
        block.r.ids <- list()
        for (l in 1:n.perm) {
            block.r.id <- list()
            for (j in 1:length(clust.sizes)) {
                block.r.id.clust <- list()
                u.id.clust <- names(which(table(id.time.pt) == 
                  clust.sizes[j]))
                r.u.id.clust <- u.id.clust[shuffle(length(u.id.clust))]
                for (k in 1:length(u.id.clust)) {
                  block.r.id.clust[[k]] <- which(id.time.pt == 
                    r.u.id.clust[k])
                }
                block.r.id[[j]] <- unlist(block.r.id.clust)
            }
            unlist.block.r.id <- rep(NA, length(id.time.pt))
            for (j in 1:length(clust.sizes)) {
                unlist.block.r.id[id.time.pt %in% names(which(table(id.time.pt) == 
                  clust.sizes[j]))] <- block.r.id[[j]]
            }
            block.r.ids[[l]] <- unlist.block.r.id
        }
        for (l in 1:n.perm) {
            p.ind[l, ] <- p.ind[l, block.r.ids[[l]]]
        }
        for (j in 1:n.perm) {
            r.s[[j]] <- r[p.ind[j, ]]
        }
    }
    if (!is.null(time.pt) & slope) {
        id.names <- names(table(id))
        id.time.pt <- rep(NA, length(id))
        for (j in 1:length(id.names)) {
            ind <- which(id == id.names[j])
            id.time.pt[ind] <- paste(c(id.names[j], as.character(time.pt[ind])), 
                collapse = ".")
        }
        noid.time.pt <- rep(NA, length(id))
        for (j in 1:length(id.names)) {
            ind <- which(id == id.names[j])
            noid.time.pt[ind] <- paste(as.character(time.pt[ind]), 
                collapse = ".")
        }
        ex.clust <- names(table(noid.time.pt))
        ids <- rep(NA, length(id))
        for (j in 1:length(id.names)) {
            ind <- which(id == id.names[j])
            ids[ind] <- id.names[j]
        }
        if (length(names(which(table(noid.time.pt) == 1))) != 
            0) {
            warn.ids <- unique(ids[noid.time.pt %in% names(which(table(noid.time.pt) == 
                1))])
            if (length(warn.ids) == 1) {
                warning("The cluster id (", warn.ids, ") is a silent block which is not exchangeable with any other blocks. It did not contribute to the analysis.")
            }
            if (length(warn.ids) > 1) {
                warning("The cluster ids (", paste(warn.ids, 
                  collapse = ", "), ") are silent blocks which are not exchangeable with any other blocks. They did not contribute to the analysis.")
            }
        }
        block.r.ids <- list()
        for (l in 1:n.perm) {
            block.r.id <- list()
            for (j in 1:length(ex.clust)) {
                block.r.id.clust <- list()
                u.id.clust <- unique(id.time.pt[which(noid.time.pt == 
                  ex.clust[j])])
                r.u.id.clust <- u.id.clust[shuffle(length(u.id.clust))]
                for (k in 1:length(u.id.clust)) {
                  block.r.id.clust[[k]] <- which(id.time.pt == 
                    r.u.id.clust[k])
                }
                block.r.id[[j]] <- unlist(block.r.id.clust)
            }
            unlist.block.r.id <- rep(NA, length(id.time.pt))
            for (j in 1:length(ex.clust)) {
                unlist.block.r.id[id.time.pt %in% unique(id.time.pt[which(noid.time.pt == 
                  ex.clust[j])])] <- block.r.id[[j]]
            }
            block.r.ids[[l]] <- unlist.block.r.id
        }
        for (j in 1:n.perm) {
            r.s[[j]] <- r[block.r.ids[[j]]]
        }
    }
    Qs <- rep(NA, length(Ks))
    for (j in 1:length(Ks)) {
        Qs[j] <- (t(r) %*% inv.vcov %*% Ks[[j]] %*% inv.vcov %*% 
            r)/(t(r) %*% inv.vcov %*% r)
    }
    Q0s <- list()
    for (j in 1:length(Ks)) {
        Q0s.inv <- rep(NA, n.perm)
        for (k in 1:n.perm) {
            Q0s.inv[k] <- (t(r.s[[k]]) %*% inv.vcov %*% Ks[[j]] %*% 
                inv.vcov %*% r.s[[k]])/(t(r.s[[k]]) %*% inv.vcov %*% 
                r.s[[k]])
        }
        Q0s[[j]] <- Q0s.inv
    }
    pvs <- rep(NA, length(Ks))
    for (j in 1:length(Ks)) {
        pvs[j] <- (length(which(Q0s[[j]] > Qs[[j]])) + 1)/(n.perm + 
            1)
    }
    T <- min(pvs)
    T0 <- rep(NA, n.perm)
    for (l in 1:n.perm) {
        T0.s.n <- list()
        for (m in 1:length(Ks)) {
            T0.s.n[[m]] <- Q0s[[m]][-l]
        }
        a.Ts <- unlist(lapply(Ks, function(x) return((t(r.s[[l]]) %*% 
            inv.vcov %*% x %*% inv.vcov %*% r.s[[l]])/(t(r.s[[l]]) %*% 
            inv.vcov %*% r.s[[l]]))))
        a.pvs <- unlist(mapply(function(x, y) (length(which(x > 
            y)) + 1)/n.perm, T0.s.n, a.Ts))
        T0[l] <- min(a.pvs)
    }
    pv.opt <- (length(which(T0 < T)) + 1)/(n.perm + 1)
    names(pvs) <- names(Ks)
    return(list(ItembyItem = pvs, aGLMMMiRKAT = pv.opt))
}

getPermuteMatrix <-
function (perm, N, strata = NULL) {
    if (length(perm) == 1) {
        perm <- how(nperm = perm)
    }
    if (!missing(strata) && !is.null(strata)) {
        if (inherits(perm, "how") && is.null(getBlocks(perm))) 
            setBlocks(perm) <- strata
    }
    if (inherits(perm, "how")) 
        perm <- shuffleSet(N, control = perm)
    else {
        if (!is.integer(perm) && !all(perm == round(perm))) 
            stop("permutation matrix must be strictly integers: use round()")
    }
    if (is.null(attr(perm, "control"))) 
        attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"), 
            nperm = nrow(perm)), class = "how")
    perm
}

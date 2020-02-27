D2K <-
function(D) {
    n <- nrow(D)
    centerM <- diag(n) - 1/n
    K <- -0.5 * centerM %*% (D * D) %*% centerM
    eK <- eigen(K, symmetric = TRUE)
    K <- eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
    return(K)
}

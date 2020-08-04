# source("http://bioconductor.org/biocLite.R")
# biocLite("survcomp")

combine_test <- function (p, weight, method = c("fisher", "z.transform", "logit"), hetero = FALSE, na.rm = FALSE) 
{
    if (hetero) {
        k <- length(p)
        if (missing(weight)) {
            weight <- rep(1, k)
        }
        cc.ix <- !is.na(p)
        if (!all(cc.ix) && !na.rm) {
            stop("missing values are present!")
        }
        p <- p[cc.ix]
        weight <- weight[cc.ix]
        z <- qnorm(p, lower.tail = FALSE)
        Q <- sum(weight * (z - mean(z))^2)
        qpv <- pchisq(Q, df = k - 1, lower.tail = FALSE)
        return(list(Q = Q, p.value = qpv))
    }
    method <- match.arg(method)
    na.ix <- is.na(p)
    if (any(na.ix) && !na.rm) {
        stop("missing values are present!")
    }
    if (all(na.ix)) {
        return(NA)
    }
    p <- p[!na.ix]
    k <- length(p)
    if (k == 1) {
        return(p)
    }
    if (missing(weight)) {
        weight <- rep(1, k)
    }
    switch(method, fisher = {
        cp <- pchisq(-2 * sum(log(p)), df = 2 * k, lower.tail = FALSE)
    }, z.transform = {
        z <- qnorm(p, lower.tail = FALSE)
        cp <- pnorm(sum(weight * z)/sqrt(sum(weight^2)), lower.tail = FALSE)
    }, logit = {
        tt <- (-sum(log(p/(1 - p))))/sqrt(k * pi^2 * (5 * k + 
            2)/(3 * (5 * k + 4)))
        cp <- pt(tt, df = 5 * k + 4, lower.tail = FALSE)
    })
    return(cp)
}

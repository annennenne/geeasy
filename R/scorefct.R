#' Internal functions for the MESS package
#'
#' @param object input geepack object from a geeglm fit.
#' @param beta The estimated parameters. If set to \code{NULL} then
#'     the parameter estimates are extracted from the model fit object
#'     object.
#' @param testidx Indices of the beta parameters that should be tested
#'     equal to zero
#' @param sas Logical. Should the SAS version of the score test be
#'     computed. Defaults to \code{FALSE}.
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @importFrom MESS qdiag
#' @keywords manip
scorefct <- function(object, beta=NULL, testidx=NULL, sas=FALSE) {

    sinv <- function(obj) {
        return(chol2inv(chol(obj)))
    }
    
#
#   Should be necessary any more
#
    # Check that ids are correctly ordered
#    if (!ordered.clusters(object$id)) {
#        stop("clusters in the gee model are not ordered and contiguous. They really should be since otherwise geepack will consider non-contiguous clusters with same id as different. Reorder your dataset or redefine the cluster id variable and run the gee fit again.")
#    }

    if (any(object$weights != 1)) {
        stop("Haven't thought about if there is a problem with weights so will not do any computations")
    }

    clusters <- unique(object$id)
    nclusters <- length(clusters)
    if (is.null(beta)) {
        beta <- coef(object)
    }

    # Offsets handled correctly?
    y <- object$y

    x <- model.matrix(object)

    linear.predictors <- x%*%beta
    if (!is.null(object$offset))
        linear.predictors <- linear.predictors + object$offset
    mui <- object$family$linkinv(linear.predictors)

    invert <- if ("MASS" %in% loadedNamespaces()) {
#        MASS::ginv
        sinv
    } else { sinv }

    myres <- lapply(clusters, function(cluster) {
        # Individiuals in cluster
        idx <- (object$id == cluster)

        # Cluster size
        csize <- sum(idx)

        # Di is r*k
#        Di <- t(x[idx,,drop=FALSE]) %*% MESS:::qdiag(object$family$mu.eta(linear.predictors[idx]), nrow=csize)
        Di <- crossprod(x[idx,,drop=FALSE], diag(object$family$mu.eta(linear.predictors[idx]), nrow=csize))
        A  <- diag(sqrt(object$family$variance(mui[idx])), nrow=csize)
        Rmat <- diag(csize)
        Ralpha <- switch(object$corstr,
                         independence = Rmat,
                         exchangeable = matrix(rep(object$geese$alpha, csize^2), csize),
                         ar1          = object$geese$alpha^abs(row(Rmat) - col(Rmat)))
        if (object$corstr=="exchangeable")
            diag(Ralpha) <- 1

        V <- outer(MESS::qdiag(A), MESS::qdiag(A)) * Ralpha * object$geese$gamma  # Ok

        # V inverse
        Vinv <- invert(V)
        DiVinv <- tcrossprod(Di, Vinv)
        list(score = DiVinv %*% (y[idx] - mui[idx]),
             DUD   = tcrossprod(DiVinv, Di))
    })

    ## Speed improvement?
    # S <- apply(sapply(myres, function(oo) oo[[1]]), 1, sum)
    S <- rowSums(sapply(myres, function(oo) oo[[1]]))

    Vsand <- Reduce("+", lapply(myres, function(oo) { tcrossprod(oo[[1]])})) # I_1
    VDUD  <- Reduce("+", lapply(myres, function(oo) { oo[[2]] })) #  I_0
    iVDUD <- invert(VDUD)

    if(is.null(testidx)) {
        cat("Should not really be here")
        as.numeric(S %*% invert(invert(VDUD) %*% Vsand %*% invert(VDUD)) %*% S / nclusters)
    } else {
        if (sas) {
            as.numeric(t(S) %*% iVDUD[,testidx] %*% invert( (iVDUD %*%  Vsand  %*% iVDUD)[testidx,testidx] ) %*% iVDUD[testidx,] %*% S)
        }
        else {
            myvar <- Vsand[testidx,testidx] - Vsand[testidx, -testidx] %*% invert(Vsand[-testidx,-testidx]) %*% Vsand[-testidx,testidx]
            as.numeric(t(S[testidx]) %*% invert( myvar ) %*% S[testidx])
        }
    }
}




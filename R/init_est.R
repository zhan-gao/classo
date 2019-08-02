#' Initial estimator
#'
#' @export init_est

init_est <- function(X, y, TT){

    p <- ncol(X)
    N <- length(y) / TT

    beta0 <- matrix(0, N, p);
    for(i in 1:N){
        ind <- ( (i-1)*TT+1 ):(i*TT);
        yy <- y[ind]
        XX <- X[ind, ]

        # beta0[i, ] <- solve( t(XX) %*% XX ) %*% ( t(XX) %*% yy )
        beta0[i, ] <- lsfit(x = XX, y = yy, intercept = FALSE)$coefficients
    }

    return(beta0)
}

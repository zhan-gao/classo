#' Half panel Jackknife
#' develop on Dhaene and Jochmans (2015)
#'
#' @param t # of time periods
#' @param y (n*t vector) response variable
#' @param x (nt * p matrix)independent variable
#'
#' @export
#'
SPJ_PLS <- function(t, y, x){

    x <- as.matrix(x)
    n <- length(y) / t
    p <- ncol(x)

    period1_i <- c(rep(1, floor(t/2)), rep(0, ceiling(t/2)))
    period1 <- as.logical( rep(period1_i, n) )
    period2 <- !period1

    theta_bar <- V <- matrix(0, p, 2)



    for(tt in c(1,2)){

        if(tt == 1){
            x_half <- x[period1, ]
            y_half <- as.matrix(y[period1])
        } else {
            x_half <- x[period2, ]
            y_half <- as.matrix(y[period2])
        }

        t_half <- length(y_half) / n

        y_demean <- demean(y_half, n, t_half)
        x_demean <- demean(x_half, n, t_half)

        # b <- lsfit(x_demean, y_demean, intercept = FALSE)$coefficients
        b <- MASS::ginv( t(x_demean) %*% x_demean ) %*% ( t(x_demean) %*% y_demean )
        theta_bar[, tt] <- b

    }

    return(rowMeans(theta_bar))

}

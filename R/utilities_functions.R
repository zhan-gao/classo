#' @export data.normalization
#' @export demean

demean <- function(yy, N, TT){

    # Output is the demeaned data with the same dimension as input
    # NT * 1 or NT * p

    if(dim(yy)[1] != N*TT) print("Error! Dimension of
                                 inputs in demean is wrong!")

    p <- dim(yy)[2];

    if( p == 1){
        y.temp <- matrix(yy, nrow = TT);
        m <- apply(y.temp, 2, mean);
        y.temp <- y.temp - matrix( rep(m, times = TT), nrow = TT,
                                   ncol = N, byrow = TRUE);
        y <- matrix(y.temp, nrow = N*TT)
        return(y)
    }
    else{
        y <- matrix(0, N*TT, p);
        for(j in 1:p){
            y.temp <- matrix( yy[,j], nrow = TT);
            m <- apply(y.temp, 2, mean);
            y.temp <- y.temp - matrix( rep(m, times = TT), nrow = TT,
                                       ncol = N, byrow = TRUE);
            y[,j] <- matrix(y.temp, nrow = N*TT);
        }
        return(y)
    }
}

data.normalization <- function(yy, N, TT){

    # Output is the demeaned data with the same dimension as input
    # NT * 1 or NT * p

    if(dim(yy)[1] != N*TT) print("Error! Dimension of
                                 inputs in demean is wrong!")

    p <- dim(yy)[2];

    if( p == 1){
        y.temp <- robustHD::standardize( matrix(yy, nrow = TT),
                                         scaleFun = function(x) {
                                             sqrt(sum((x - mean(x)) ^ 2) / length(x))
                                         })
        y.mean <- matrix( rep( colMeans( matrix(yy, nrow = TT)), TT), nrow = TT, byrow = TRUE)
        y.raw.temp <- y.temp + y.mean

        y <- matrix(y.temp, nrow = N * TT)
        y.raw <- matrix(y.raw.temp, nrow = N*TT)

    }
    else{
        y <- matrix(0, N*TT, p);
        y.raw <- matrix(0, N*TT, p);
        for(j in 1:p){
            y.temp <- robustHD::standardize( matrix(yy[, j], nrow = TT),
                                             scaleFun = function(x) {
                                                 sqrt(sum((x - mean(x)) ^ 2) / length(x))
                                             })

            y.mean <- matrix( rep(colMeans( matrix(yy[, j], nrow = TT)), TT), nrow = TT, byrow = TRUE)
            y.raw.temp <- y.temp + y.mean

            y[, j] <- matrix(y.temp, nrow = N * TT)
            y.raw[, j] <- matrix(y.raw.temp, nrow = N*TT)

        }
    }

    return(list(y = y, y.raw = y.raw))
}

group.coerce <- function(group.est, a.out, group0, a0, N, N.frac, K, p){
    
    a.est <- matrix(0,K,p);
    group <- rep(0,N);
    
    for(k in 1:K){
        dif <- apply( (matrix( a.out[k,], K, p, byrow = TRUE) - a0)^2, 1, sum);
        gg <- which.min(dif);
        
        a.est[gg, ] <- a.out[k, ];
        group[group.est == k] <- gg;
    }
    
    if( sum(group == 0) ){
        noquote("Error: some inidividuals are not classfied!")
    }
    
    # correct ratio
    ratio <- sum( group == group0 ) / N;
    # rmse
    weighted.se <- sum( ( a.est[,1] - a0[,1] )^2 * N.frac );
    
    return( list(ratio = ratio, se = weighted.se,
                 a = a.est, group = group) );
}


data.normalization <- function(yy, N, TT){
    
    # Output is the demeaned data with the same dimension as input
    # NT * 1 or NT * p
    
    require("robustHD")
    
    if(dim(yy)[1] != N*TT) print("Error! Dimension of
                                 inputs in demean is wrong!")
    
    p <- dim(yy)[2];
    
    if( p == 1){
        y.temp <- standardize( matrix(yy, nrow = TT) );
        y <- matrix(y.temp, nrow = N*TT)
        return(y)
    }
    else{
        y <- matrix(0, N*TT, p);
        for(j in 1:p){
            y.temp <- standardize( matrix( yy[,j], nrow = TT) );
            y[,j] <- matrix(y.temp, nrow = N*TT);
        }
        return(y)
    }
}
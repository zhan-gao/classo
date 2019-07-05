#' PLS estimation by the iterative algorithm via Rmosek
#'
#' @param N individual dimension
#' @param TT time dimension
#' @param y y(TN * 1)
#' @param X X(TN * P)
#' @param K number of groups
#' @param lambda Tuning variable
#' @param R Maximum # of iteration
#' @param tol convergence criterion
#'
#' @return A list contains estimated coeffcients and group struncture
#' \item{b.est}{N * p matrix containing estimated slope for each individual}
#' \item{a.out}{K * p matrix containing estimated slope for each group}
#' \item{group.est}{group_id for each individual}
#' \item{converge}{A boolean indicating whether convergence criteria is met }
#'
#' @export
#'
#'
PLS.mosek <- function(N, TT, y, X, K, lambda, R, tol = 1e-4){

    p <- dim(X)[2];

    # Use individual regression result as the initial value
    beta0 = matrix(0, N, p);
    for(i in 1:N){
        ind <- ( (i-1)*TT+1 ):(i*TT);
        yy <- y[ind, ]
        XX <- X[ind, ]
        beta0[i, ] <- solve( t(XX) %*% XX ) %*% ( t(XX) %*% yy );
    }



    b.out <- array( beta0, c(N,p,K) );
    a.out <- matrix(0,K,p);

    b.old <- matrix(1,N,p);
    a.old <- matrix(1,1,p);

    for(r in 1:R){

        for(k in 1:K){

            # N * 1: consider it as gamma
            penalty.out <- pen.generate(b.out, a.out, N, p, K, k);
            # do optimization
            mosek.out <- opt.mosek(y, X, penalty.out, N, TT,
                                   K, p, lambda);
            a.out[k,] <- mosek.out$alpha;
            b.out[ , ,k] <- matrix(mosek.out$beta, N, p, byrow = TRUE);

        }

        # Check the convergence criterion
        a.new <- a.out[K,];
        b.new <- b.out[ , ,K];

        if(criterion(a.old,a.new,b.old,b.new,tol)){
            break;
        }
        # Update
        a.old <- a.out[K,];
        b.old <- b.out[ , ,K];
    }

    # put b.out to nearest a.out and get the group estimation

    a.out.exp <- aperm( array(a.out, c(K,p,N)), c(3,2,1) );
    d.temp <- (b.out - a.out.exp)^2;
    dist <- sqrt( apply(d.temp, c(1,3), sum) );
    group.est <- apply(dist,1,which.min);


    if( sum( as.numeric( table(group.est) ) > p ) == K ){
        # Post-Lasso estimation
        a.out <- post.lasso( group.est, y, X, K, p, N, TT );
    }


    b.est <- matrix(999, N, p);
    for(i in 1:N){
        group <- group.est[i];
        b.est[i,] <- a.out[group, ];
    }

    result <- list( b.est = b.est, a.out = a.out,
                    group.est = group.est, converge = (r < R));

    return(result)
}

#' PLS estimation by the iterative algorithm via CVXR
#'
#' @inheritParams PLS.mosek
#' @param solver "ECOS" or "MOSEK"
#'
#' @return A list contains estimated coeffcients and group struncture
#' \item{b.est}{N * p matrix containing estimated slope for each individual}
#' \item{a.out}{K * p matrix containing estimated slope for each group}
#' \item{group.est}{group_id for each individual}
#' \item{converge}{A boolean indicating whether convergence criteria is met }
#'
#' @export
#'


PLS.cvxr <- function(N, TT, y, X, K, lambda, R, tol = 1e-4, solver = "ECOS"){

    p <- dim(X)[2];

    # Use individual regression result as the initial value
    beta0 = matrix(0, N, p);
    for(i in 1:N){
        ind <- ( (i-1)*TT+1 ):(i*TT);
        yy <- y[ind, ]
        XX <- X[ind, ]
        beta0[i, ] <- solve( t(XX) %*% XX ) %*% ( t(XX) %*% yy );
    }



    b.out <- array( beta0, c(N,p,K) );
    a.out <- matrix(0,K,p);

    b.old <- matrix(1,N,p);
    a.old <- matrix(1,1,p);

    for(r in 1:R){

        for(k in 1:K){

            # print(c(r,k))

            # N * 1: consider it as gamma
            gamma <- pen.generate(b.out, a.out, N, p, K, k)

            # Commented out: Suggested by Dr. Narasimhan
            # b = Variable(N*p)
            # a = Variable(p)
            #
            # obj = 0
            # End Commented out

            X.list = list()
            for(i in 1:N){
                ind = ( (i-1)*TT + 1 ):(i*TT)
                id = ( (i-1)*p + 1 ):(i*p)
                X.list[[i]] =  X[ind, ]
                # Commented out: Suggested by Dr. Narasimhan
                # obj = obj + gamma[i] * norm2( b[id] - a )
                # End Commented out
            }

            ## Code added
            b = Variable(p, N)
            a = Variable(p)
            A <- matrix(1, nrow = 1, ncol = N)
            obj1 <- norm2(b - a %*% A, axis = 2) %*% gamma
            ## End Code added

            XX = bdiag(X.list)

            ## Original commented out
            ## obj = Minimize( sum_squares(y - XX %*% b)/(N * TT) + obj*(lambda/N) )
            ## End Original commented out

            ## Code added and modified
            obj = Minimize( sum_squares(y - XX %*% vec(b))/(N * TT) + obj1*(lambda/N) )
            ## End Code added and modified
            Prob = Problem(obj)
            cvxr.out = solve(Prob, solver = solver)

            a.out[k, ] = cvxr.out$getValue(a)
            b.out[ , , k] = matrix( cvxr.out$getValue(b), N, p, byrow = TRUE)

        }

        # Check the convergence criterion
        a.new <- a.out[K,];
        b.new <- b.out[ , ,K];

        if(criterion(a.old,a.new,b.old,b.new,tol)){
            break;
        }
        # Update
        a.old <- a.out[K,];
        b.old <- b.out[ , ,K];
    }

    # put b.out to nearest a.out and get the group estimation

    a.out.exp <- aperm( array(a.out, c(K,p,N)), c(3,2,1) );
    d.temp <- (b.out - a.out.exp)^2;
    dist <- sqrt( apply(d.temp, c(1,3), sum) );
    group.est <- apply(dist,1,which.min);


    if( sum( as.numeric( table(group.est) ) > p ) == K ){
        # Post-Lasso estimation
        a.out <- post.lasso( group.est, y, X, K, p, N, TT );
    }


    b.est <- matrix(999, N, p);
    for(i in 1:N){
        group <- group.est[i];
        b.est[i,] <- a.out[group, ];
    }

    result <- list( b.est = b.est, a.out = a.out,
                    group.est = group.est, converge = (r < R));

    return(result)
}


##########################################################
pen.generate <- function(b, a, N, p, K, kk){

    # generate the known part of the penalty term
    # output a N*1 vector

    a.out.exp <- aperm( array(a, c(K,p,N)), c(3,2,1) );
    p.temp <- ( b - a.out.exp )^2;
    p.norm <- sqrt( apply(p.temp, c(1,3), sum) );

    ind <- setdiff(1:K,kk);

    pen <- apply(p.norm[, ind], 1, prod);
    return(pen)

}

##########################################################
opt.mosek <- function(y, X, penalty, N, TT, K, p, lambda){

    # call Mosek solver to solve the optimization
    # conic programming


    # INPUT Arg:
    #   dimensions N, TT, K, p
    #   tuning parameter lambda
    #   data y(TN * 1), X(TN * P)
    #   parameter penalty (N*1) numeric


    # set sense of optim and tolerance
    prob <- list(sense = "min");
    prob$dparam$intpnt_nl_tol_rel_gap <- 1e-5;

    # objective: coeffiects c
    # order of variables:
    # beta_i (i = 1,2,..N), nu_i (i=1,2,...,N) , mu_i (i=1,2,...,N)
    # alpha_k, s_i (i=1,2,...,N), r_i (i=1,2,...,N), t_i (i=1,2,...,N),
    # w_i (i=1,2,...,N)

    prob$c <- c( rep( 0, N*(2*p+TT+2)+p ), rep( 1/(N*TT), N ),
                 penalty*lambda/N );

    # lieanr constraint: matrix A

    # There must be some smart methods to split the matrix
    # without invoke loops. Try to modify it later
    X.split <- list();
    for(i in 1:N){
        ind <- ( (i-1)*TT + 1 ):(i*TT);
        X.split[[i]] <-  X[ind, ];
    }

    A.y <- cbind( bdiag(X.split), Diagonal(TT*N),
                        Matrix(0, TT*N, N*(p+4)+p));

    A.0 <- cbind( Diagonal(N*p), Matrix(0, N*p, TT*N),
                  -Diagonal(N*p), -matrix(diag(p),N*p,p,byrow = TRUE),
                  Matrix(0,N*p,N*4));

    A.nhalf <- cbind( Matrix(0,N,N*(2*p+TT)+p), Diagonal(N),
                      Matrix(0,N,N), -Diagonal(N)/2, Matrix(0,N,N) );

    A.phalf <- cbind( Matrix(0,N,N*(2*p+TT)+p), Matrix(0,N,N),
                      Diagonal(N), -Diagonal(N)/2, Matrix(0,N,N) );

    A <- rbind(A.y, A.0, A.nhalf, A.phalf);
    prob$A <- as(A,"CsparseMatrix");

    # linear constraint: upperbound and lowerbound
    prob$bc <- rbind( blc = c( y, rep(0,N*p), rep(-1/2,N), rep(1/2,N) ),
                      buc = c( y, rep(0,N*p), rep(-1/2,N), rep(1/2,N) ) );
    # t_i \geq 0
    prob$bx <- rbind( blx = c( rep(-Inf,N*(2*p+TT+2)+p), rep(0,N)
                               , rep(-Inf,N) ),
                      bux = c( rep(Inf,N*(2*p+TT+4)+p) ) );

    # Conic constraints
    CC <- list();
    bench <- N*(2*p + TT) +p;

    for(i in 1:N){
        s.i <- bench + i;
        r.i <- bench + N + i;
        nu.i <- (N*p + (i-1)*TT + 1):(N*p + i*TT);
        w.i <- bench + 3*N + i;
        mu.i <- (N*(TT+p) + (i-1)*p + 1):(N*(TT+p) + i*p);
        CC <- cbind(CC, list("QUAD", c(r.i, nu.i, s.i)),
                     list("QUAD", c(w.i, mu.i)) );
    }
    prob$cones <- CC;
    rownames(prob$cones) <- c("type", "sub");

    # Invoke mosek solver

    mosek.out <- mosek(prob, opts = list(verbose = 0) );

    est <- mosek.out$sol$itr$xx;
    beta <- est[1:(N*p)];
    alpha <- est[( N*(2*p+TT)+1 ):( N*(2*p+TT)+p )];
    result <- list(beta = beta, alpha = alpha);
    return(result)
}

##########################################################
criterion <- function( a.old, a.new, b.old, b.new, tol ){

    d <- FALSE;

    a.cri <- sum( abs( a.old - a.new ) ) / (sum( abs(a.old) ) + 0.0001);
    b.cri <- mean( abs( b.old - b.new ) ) / (mean( abs(b.old) ) + 0.0001);

    if(a.cri < tol & b.cri < tol){
        d <- TRUE;
    }

    return(d)

}

##########################################################
post.lasso <- function( group.est, y, X,  K, p, N, TT ){

    a.out <- matrix(0,K,p);

    for( k in 1:K ){
        Ind <- 1:N; group.ind <- Ind[group.est == k];
        data.ind <- as.numeric( sapply(group.ind, extend.ind) );
        yy <- y[data.ind,];
        XX <- X[data.ind,];

        a.out[k, ] <- solve( t(XX) %*% XX ) %*% ( t(XX) %*% yy );
    }
    return(a.out)
}

extend.ind <- function(i){
    return( ((i-1)*TT+1):(i*TT) );
}

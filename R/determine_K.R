#' Determine number of groups using information criterion
#'
#' @param N individual dimension
#' @param TT time dimension
#' @param y y(TN * 1)
#' @param X X(TN * P)
#' @param y_raw y_raw(TN * 1) raw data without standardization
#' @param X_raw X_raw(TN * P) raw data without standardization
#' @param lambda_seq Candidate tuning variable
#' @param K_max Maximum number of groups
#' @param rho Tuning parameter in the IC
#' @param FUN choose which function to be used in PLS estimation
#' @param beta0 N*p matrix. initial estimator
#' @param MaxIter Maximum # of iteration
#' @param tol convergence criterion
#'
#' @return A list contains optimal K and lambda
#' \item{lambda}{Optimal lambda}
#' \item{K}{Optimal K}
#'
#' @export

determine_K <- function(N,
                        TT,
                        y,
                        X,
                        y_raw,
                        X_raw,
                        lambda_seq,
                        K_max,
                        rho = 2 / 3 * (N * TT)^(-0.5),
                        FUN = PLS.cvxr,
                        beta0 = NULL,
                        MaxIter = 500,
                        tol = 1e-4) {

    p <- ncol(X)

    num_lambda <- length(lambda_seq)
    IC_total <- matrix(0, K_max, num_lambda)

    for(ll in 1:num_lambda){
        for(K in 1:K_max){

            print(paste(as.character(ll),
                        "/",
                        as.character(num_lambda),
                        "th parameter; K = ",
                        as.character(K),
                        "/",
                        as.character(K_max)))

            lambda <- lambda_seq[ll]

            if(K == 1){
                a <- lsfit(X, y, intercept = FALSE)$coefficients
                bias <- SPJ_PLS(TT, y_raw, X_raw)
                a_corr <- 2*a - bias

                IC_total[K, ] <- mean( (y - X %*% a_corr)^2 )

                next
            }

            pls_out <- FUN(N,
                           TT,
                           y,
                           X,
                           K,
                           lambda,
                           beta0 = beta0,
                           R = MaxIter,
                           tol = tol,
                           post_est = FALSE,
                           bias_corr = FALSE)

            Q <- rep(1e10, K)

            # Post-estimation
            for (k in 1:K){

                group_k <- (pls_out$group.est == k)

                if(sum(group_k) > 2*p/TT){

                    Ind <- 1:N
                    group_ind <- Ind[group_k]
                    data_ind <- as.numeric( sapply(group_ind, function(i){((i-1)*TT+1):(i*TT)}) )
                    yy_k <- y[data_ind]
                    XX_k <- X[data_ind, ]
                    yy_raw_k <- y_raw[data_ind]
                    XX_raw_k <- X_raw[data_ind, ]

                    # bias correction
                    bias_k <- SPJ_PLS(TT, yy_raw_k, XX_raw_k)
                    a_k <- lsfit(XX_k, yy_k, intercept = FALSE)$coefficients
                    a_corr_k <- 2*a_k - bias_k

                } else {
                    a_corr_k <- pls_out$a.out[k, ]
                }

                Q[k] <- sum( (yy_k - XX_k %*% a_corr_k)^2 )
            }

            IC_total[K, ll] <- sum(Q) / (N*TT)

        }
    }

    IC <- log(IC_total) + rho * p * matrix(rep(1:K_max, num_lambda), nrow = K_max)


    min_ind <- which.min(IC)
    lambda_opt <- lambda_seq[ceiling(min_ind / K_max)]
    K_temp <- min_ind %% K_max
    if (K_temp == 0) K_temp <- K_max
    K_opt <- K_temp

    return(list(lambda = lambda_opt, K = K_opt))

}

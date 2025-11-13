# Functions for finding posterior modes

#' Sample n times from a multivariate normal where ||x[1:dim]|| < pi. Uses accept/reject strategy.
#' We call this a circular multivariate normal distribution
#'
#' @keywords internal
#' @noRd
rcnorm_accept_reject <- function(n, mu, Sigma, r, dim = length(mu)) {

  # sampled_mat will contain samples
  num_accepted <- 0
  p <- nrow(Sigma)
  sampled_mat <- matrix(rep(0, p * n), nrow = n)

  # Sample from multivariate normal and see if norm is less than r. If yes, add to
  # sample_mat
  while (num_accepted < n) {

    one_sample <- mvtnorm::rmvnorm(1, mu, Sigma)

    if (sum(one_sample[1:dim]^2) < r^2) {
      num_accepted <- num_accepted + 1
      sampled_mat[num_accepted, ] <- one_sample
    }

  }

  return(sampled_mat)
}



#' Uses BFGS to find local minimum of negative log likelihood or negative log likelihood
#' + w_0 *log prior
#'
#' @keywords internal
#' @noRd
find_local_minimum <- function(input_X, input_Y, num_betas, first_direction,
                               zero_is_zero, interpolation, b_vec,
                               prior_mean, prior_sd, w_0) {

  found_result <- FALSE

  # Keeps trying until local mode is reached
  while(!found_result){

    bfgs_result <- tryCatch({

      initial_point <- rcnorm_accept_reject(1, prior_mean, diag(prior_sd^2),
                                            pi, num_betas)

      # BFGS Optimization with initial starting point
      final_res <- stats::optim(as.vector(initial_point),
                         fn = neg_bayes_value,
                         method = 'BFGS', control = list(maxit = 10000),
                         input_X = input_X, input_Y = input_Y,
                         num_betas = num_betas, first_direction = first_direction,
                         zero_is_zero = zero_is_zero, interpolation = interpolation,
                         b_vec = b_vec, prior_mean = prior_mean, prior_sd = prior_sd,
                         weight_vector = NULL,
                         w_0 = w_0,
                         hessian = TRUE)


      if (sum(eigen(final_res$hessian)$values <= 0) == 0){
        found_result <- TRUE
      }

      # Return parameter, the optimization cost, and Inverse of Hessian
      # (Laplace Approximation)
      list(min_par = final_res$par,
           min_cost = final_res$value,
           laplace_approx = solve(final_res$hessian))

    },  error = function(e){


    })

  }

  return(bfgs_result)
}

#' Runs BFGS multiple times to find modes of likelihood or posterior distribution.
#' Returns the modes, likelihood value and inverse Hessian at each point
#' (Uses parallelization)
#'
#' @keywords internal
#' @noRd
find_modes <- function(num_searches, input_X, input_Y, num_betas,
                       first_direction, zero_is_zero, interpolation, b_vec,
                       prior_mean, prior_sd, w_0 = 1, num_cores){

  `%dopar%` <- foreach::`%dopar%`

  print("Finding Modes:")
  if(num_cores == "ALL"){
    num_cores <- parallel::detectCores()
  }

  cl <- parallel::makeCluster(num_cores)


  doSNOW::registerDoSNOW(cl)
  funcs <- ls(envir = .GlobalEnv)

  pb <- utils::txtProgressBar(max = num_searches, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Run BFGS multiple times using parallelization
  mode_samples <- foreach::foreach(i = 1:num_searches, .packages = c('pracma', 'stats', 'mvtnorm' ),
                              .export = funcs, .options.snow = opts) %dopar% {


                                mode_result <- find_local_minimum(input_X, input_Y, num_betas, first_direction,
                                                                  zero_is_zero, interpolation, b_vec,
                                                                  prior_mean, prior_sd, w_0)

                                mode_result


                              }

  close(pb)
  parallel::stopCluster(cl)

  return(mode_samples)

}


# Functions for MCMC (Simple Metrpolis-Hasitngs)


#' Obtain the normalizing constant (1/ P(||x[1:dim]|| < r)) for a circular normal distribution
#' Needed for log likelihood calculations for MCMC
#'
#' @keywords internal
#' @noRd
get_normalize_const_cnorm <- function(mu, Sigma, r, dim) {

  eigen_decomp <- eigen(Sigma[1:dim, 1:dim])

  B_matrix <- eigen_decomp$vectors %*% sqrt(diag(eigen_decomp$values))

  a_vector <- solve(B_matrix, mu[1:dim])

  prob <- 1 - CompQuadForm::imhof(r^2, eigen_decomp$values, h = rep(1, dim),
                    delta = a_vector^2,
                    epsabs = 10^(-6), epsrel = 10^(-6), limit = 10000)$Qq

  return(1/prob)
}

#' Obtains n chain of a single chain of Metropolis Hastings with proposal N(x_{t -1}, Sigma)
#' for the Diffeomorphism Model. proposal_cov is usually scaled, and proposal mean refers to
#' the first sample
#'
#' @keywords internal
#' @noRd
mcmc_one_chain <- function(n, proposal_mean, proposal_cov, input_X, input_Y, num_betas,
                           first_direction, zero_is_zero, interpolation, b_vec,
                           prior_mean, prior_sd, dim_correction = TRUE,
                           dim_parameter = 1) {

  # Apply correction to proposal covariance
  if(dim_correction) {
    proposal_cov <- dim_parameter * (2.4^2)/nrow(proposal_cov) * proposal_cov
  }

  # Initialize MCMC samples
  mcmc_samples <- matrix(rep(0, n *nrow(proposal_cov) ), nrow = n)

  # First sample, using posterioal mean
  old_samp <- rcnorm_accept_reject(1, proposal_mean, proposal_cov, pi, num_betas)

  # Log Posterior of old sample
  old_post_like <- bayes_value(old_samp[1, ], input_X, input_Y, num_betas,
                               first_direction, zero_is_zero,
                               interpolation, b_vec, prior_mean, prior_sd)

  num_accepted <- 0
  for (i in 1:n){

    # Get new sample
    new_samp <- rcnorm_accept_reject(1, old_samp, proposal_cov, pi, num_betas)

    # Log Posterior of new sample
    new_post_like <- bayes_value(new_samp[1, ], input_X, input_Y, num_betas,
                                 first_direction, zero_is_zero,
                                 interpolation, b_vec, prior_mean, prior_sd)

    new_given_old <- log(get_normalize_const_cnorm(old_samp, proposal_cov, pi, num_betas))

    old_given_new <- log(get_normalize_const_cnorm(new_samp, proposal_cov, pi, num_betas))

    # Calculate Metropolis Hastings ratio
    mh_ratio <- exp(new_post_like + old_given_new - old_post_like - new_given_old)

    if(stats::runif(1) <= mh_ratio){

      # If accepted, update sample
      old_samp <- new_samp
      old_post_like <- new_post_like

      num_accepted <-  num_accepted + 1
    }

    mcmc_samples[i, ] <- old_samp[1, ]

  }
  # Return mcmc samples acceptance rate
  return(list(mcmc_samples = mcmc_samples, acp_prob = num_accepted/n))
}

#' Run multiple parallel chains of Metropolis Hastings for the Diffeomorphism Model
#'
#' @keywords internal
#' @noRd
mcmc_parallel <- function(posterior_modes, n_per_chain, input_X, input_Y, num_betas,
                          first_direction, zero_is_zero, interpolation, b_vec,
                          prior_mean, prior_sd, dim_parameter) {

  `%dopar%` <- foreach::`%dopar%`

  # Make a progress bar for MCMC Chain
  print("MCMC Progress:")
  cl <- parallel::makeCluster(parallel::detectCores())
  doSNOW::registerDoSNOW(cl)
  funcs <- ls(envir = .GlobalEnv)

  pb <- utils::txtProgressBar(max = length(posterior_modes), style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Obtain MCMC Chains with Parallelizations
  mcmc_samples_all <- foreach::foreach(i = 1:length(posterior_modes), .packages = c('mvtnorm', 'pracma',
                                                                       'CompQuadForm', 'MASS', 'stats'),
                              .export = funcs, .options.snow = opts) %dopar% {

                                mcmc_samples <- mcmc_one_chain(n_per_chain, posterior_modes[[i]]$min_par,
                                                               posterior_modes[[i]]$laplace_approx,
                                                               input_X, input_Y, num_betas,
                                                               first_direction, zero_is_zero, interpolation, b_vec,
                                                               prior_mean, prior_sd,
                                                               dim_correction = TRUE,
                                                               dim_parameter = dim_parameter)

                                mcmc_samples

                              }

  close(pb)
  parallel::stopCluster(cl)

  # Get MCMC chains and Acceptance Rates
  mcmc_samples <- lapply(mcmc_samples_all, function(x){x$mcmc_samples})
  acceptance_probs <- do.call(rbind, lapply(mcmc_samples_all, function(x){x$acp_prob}))

  return(list(mcmc_samples = mcmc_samples, acceptance_probs = acceptance_probs))
}













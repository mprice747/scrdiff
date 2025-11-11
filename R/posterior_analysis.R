# Posterior Analysis for Diffeomorphism Method

#' Calculate prob% HPD Intervals for a vector of samples
#'
#' @keywords internal
#' @noRd
compute_HPD <- function(samples, prob = 0.95) {

  n <- length(samples)
  sorted <- sort(samples)
  m <- floor(prob * n)

  intervals <- sapply(1:(n - m), function(i) sorted[i + m] - sorted[i])
  min_idx <- which.min(intervals)

  return(c(sorted[min_idx], sorted[min_idx + m]))
}


#' Make stationary points HPD Intervals (90%, 95%, 99% and 95% Bonferonni Corrected Intervals)
#'
#' @keywords internal
#' @noRd
make_SP_hpd_intervals <- function(stat_points){

  # Have dimension names for matrix
  dimnames_list <- list(unlist(lapply(seq(1, ncol(stat_points)),
                                      function(x){paste0("S", x)})),
                        c("lower", "upper"))

  # Get 90, 95, 99 and Bonferonni Intervals
  intervals_90 <- t(apply(stat_points, 2, compute_HPD, prob = 0.90))
  dimnames(intervals_90) <- dimnames_list

  intervals_95 <- t(apply(stat_points, 2, compute_HPD, prob = 0.95))
  dimnames(intervals_95) <- dimnames_list

  intervals_99 <- t(apply(stat_points, 2, compute_HPD, prob = 0.99))
  dimnames(intervals_99) <- dimnames_list

  bonf_corr <- 1 - 0.05/ncol(stat_points)

  intervals_bonf <- t(apply(stat_points, 2, compute_HPD, prob = bonf_corr))
  dimnames(intervals_bonf) <- dimnames_list

  intervals_list <- list(intervals_90 = intervals_90,
                         intervals_95 = intervals_95,
                         intervals_99 = intervals_99,
                         intervals_bonf = intervals_bonf)

  return(intervals_list)

}

#' Given raw posterior samples for diffeomorphism method, return stationary points
#'
#' @keywords internal
#' @noRd
get_stat_points_predictions <- function(posterior_samples, num_betas,
                                        first_direction, zero_is_zero, interpolation,
                                        b_vec){


  # Transform posterior samples to original parameterization
  post_process_posterior <- t(apply(posterior_samples, 1,
                                    process_post_samples_as_vec,
                                    num_betas = num_betas,
                                    first_direction = first_direction,
                                    zero_is_zero = zero_is_zero))

  total_samples <- nrow(post_process_posterior)
  ncol_post <- ncol(post_process_posterior)

  # Get Betas, Lambdas and B points
  final_Betas <- post_process_posterior[, 1:num_betas]
  final_Lambdas <- post_process_posterior[, (num_betas + 1):(ncol_post - 1)]

  if(is.null(dim(final_Betas))){
    final_Betas <- matrix(final_Betas, nrow = 1)
    final_Lambdas <- matrix(final_Lambdas, nrow = 1)
  }
  if (is.null(b_vec)){
    final_BPoints <- matrix(rep(seq(0, 1, length.out = ncol(final_Lambdas)), total_samples),
                            nrow = total_samples, byrow = TRUE)
  } else{

    final_BPoints <- matrix(rep(b_vec, total_samples),
                            nrow = total_samples, byrow = TRUE)
  }


  # Obtain Posterior predictions for a grid of points
  posterior_predictions <- diffeo_predictions(seq(0, 1, length.out = 500),
                                              final_Betas, final_BPoints,
                                              final_Lambdas,
                                              interpolation = interpolation,
                                              transform_X = FALSE,
                                              return_diffeo = TRUE)
  # Get mean posterior predictive values
  b_vals <- final_BPoints[1, 2:(ncol(final_BPoints) - 1)]

  # Get new stationary points via interpolation
  stat_points_fun <- function(diffeo_res){
    return(pracma::interp1(x = diffeo_res, y = seq(0, 1, length.out = 500),
                   xi = b_vals, method = 'linear'))
  }

  stationary_points <- t(apply(posterior_predictions$diffeomorphism,
                               2, stat_points_fun))

  # Return stationary points and posterior predictions
  return(stationary_points)
}

#' Given posterior modes, sample possible starting proposal means using a softmax
#' of the posterior likelihood
#'
#' @keywords internal
#' @noRd
apply_softmax_sampling <- function(mode_samples, num_sps){

  # Get posterior likelihood and apply softmax
  len_ms <- length(mode_samples)
  num_params <- length(mode_samples[[1]]$min_par)
  log_like_modes <- do.call(rbind, lapply(mode_samples, function(x){-1 * x$min_cost}))

  softmax_post_mass <- softmax(log_like_modes[, 1])

  # Sample posterior modes using softmax probabilities
  subsamp_modes <- sample.int(len_ms, num_sps, TRUE, softmax_post_mass)
  softmax_starting_points <- list()

  for (i in 1:len_ms){
    softmax_starting_points[[i]] <- mode_samples[[subsamp_modes[i]]]
  }

  return(softmax_starting_points)

}

#' Finds MAP of diffeomorphism model. Runs BFGS multiple times and takes value
#' with minimum cost
#'
#' @keywords internal
#' @noRd
find_map <- function(num_searches, input_X, input_Y, num_betas, first_direction,
                     zero_is_zero, interpolation, b_vec,
                     prior_mean, prior_sd) {

  # Find modes and the one with the lowest cost
  modes_found <- find_modes(num_searches, input_X, input_Y, num_betas, first_direction,
                            zero_is_zero, interpolation, b_vec,
                            prior_mean, prior_sd, 1)


  min_list <- modes_found[[which.min(sapply(modes_found, function(x) x$min_cost))]]

  map_stat_points <- get_stat_points_predictions(matrix(min_list$min_par, nrow = 1), num_betas,
                                                 first_direction, zero_is_zero, interpolation,
                                                 b_vec )[1, ]

  return(list(modes_found = modes_found, map_par = min_list$min_par,
              map_stat_points = map_stat_points))

}

#' Sample from Stationary Point Joint Posterior for Diffeomorphism Model
#'
#' Sample from posterior of Diffeomorphism Model using multiple chains of MCMC Model, with
#' Laplace approximations at the posterior modes as the covariance matrices. Low probability
#' modes are filtered out via a softmax correction of the posterior likelihood. The posterior
#' stationary points are outputted
#' @param X length n vector, x points (input data)
#' @param Y length n vector, y points (output data)
#' @param num_betas positive integer, number of weight parameters for the diffeomorphism
#' @param num_stationary positive integer, number of stationary points
#' @param first_direction 1 or -1, 1 if the function is increasing between min(X) and the first
#' stationary point, -1 if the function is decreasing
#' @param zero_is_zero boolean, TRUE if f(min(X)) = 0, FALSE otherwise and f(min(X)) needs to be
#' estimated
#' @param interpolation 'cubic' or 'linear', referring to interpolation type of lambda height vector
#' @param b_vec NULL or vector of length num_stationary + 2, nodes of template function. If NULL, nodes
#' will be equally spaced from eqach other
#' @param prior_mean p dimensional vector of prior mean (Joint Independent Normal)
#' @param prior_sd p dimension vector of prior sds (Joint Independent Normal)
#' @param n_chains positive integer, number of MCMC chains
#' @param n_per_chain positive integer, number of samples per MCMC chain
#' @param cov_scale - positive real number, scaling factor for proposal covariance matrices of MCMC chains
#' @param apply_sm_correction boolean, whether to filter out low probability posterior modes via softmax correction
#' @returns A list with the components
#' \describe{
#'   \item{diff_stat_points}{List of length \code{num_stationary}, posterior samples for each stationary point}
#'   \item{intervals_list_diff}{List of matrices, 90, 95, 99 and Bonferroni Corrected Intervals for each stationary point}
#'   \item{diff_post_cov}{Posterior covariance of stationary points}
#'   \item{diff_map_par}{MAP of unconstrained parameterization}
#'   \item{diff_map_stat_points}{MAP of the stationary points}
#'   \item{map_predictions}{Matrix with X grid and MAP estimate of f(X)}
#'   \item{posterior_modes}{List of posterior modes, likelihoods, inverse Hessians}
#'   \item{acceptance_rates}{Vector of acceptance rates for MCMC chains}
#' }
#' @export
posterior_stat_points_diff_mcmc <- function(X, Y, num_betas, num_stationary, first_direction,
                                            zero_is_zero, interpolation, b_vec,
                                            prior_mean, prior_sd,
                                       n_chains, n_per_chain, cov_scale,
                                       apply_sm_correction = TRUE) {

  # Transform Data to be [0, 1]
  X_trans <- min_max_transform_X(X)
  input_X <- X_trans$new_X
  input_Y <- Y

  # Transform b_vec if necessary
  if(!is.null(b_vec)){
    b_vec <- (b_vec - X_trans$lower_bound)/(X_trans$upper_bound - X_trans$lower_bound)
  } else{
    b_vec <- seq(0, 1, length.out = num_stationary + 2)
  }

  # Get MAP result
  map_result <- find_map(n_chains, input_X, input_Y, num_betas, first_direction,
                         zero_is_zero, interpolation, b_vec,
                         prior_mean, prior_sd)


  # Get posterior modes and the inverse Hessians
  posterior_modes <- map_result$modes_found

  if(apply_sm_correction){
    posterior_modes <- apply_softmax_sampling(posterior_modes, 100)
  }

  # Sample from posterior using standard Metropolis Hastings, using inverse Hessians as
  # proposal covariance matrices
  mcmc_samples <- mcmc_parallel(posterior_modes, n_per_chain, input_X, input_Y, num_betas,
                                first_direction, zero_is_zero, interpolation, b_vec,
                                prior_mean, prior_sd, cov_scale)

  print("Calculating Stationary Points")
  resampled <- do.call(rbind, mcmc_samples$mcmc_samples)

  # Get stationary points from posterior samples
  stat_points_mcmc <- get_stat_points_predictions(resampled, num_betas,
                                                    first_direction, zero_is_zero,
                                                    interpolation, b_vec)

  diff_stat_points_mat <- as.matrix(apply(stat_points_mcmc, 2, reverse_min_max_transform_X,
                            lower_bound = X_trans$lower_bound,
                            upper_bound = X_trans$upper_bound))

  diff_map_stat_points <- unlist(lapply(map_result$map_stat_points,
                                        reverse_min_max_transform_X,
         lower_bound = X_trans$lower_bound,
         upper_bound = X_trans$upper_bound))

  # Make intervals
  intervals_list_diff <- make_SP_hpd_intervals(diff_stat_points_mat)

  # Sampled stationary points
  diff_stat_points <- list()

  for (i in 1:num_stationary){
    diff_stat_points[[i]] <- diff_stat_points_mat[, i]
  }

  # Get predictions at certain X using MAP estimates
  X_for_pred <- seq(0, 1, length.out = 1000)
  map_post_process <- process_post_samples(map_result$map_par, num_betas,
                                           first_direction)
  map_predictions_Y <- diffeo_predictions(X_for_pred,
                                        map_post_process$betas,
                                        b_vec,
                                        map_post_process$lambdas)$predictions[, 1]

  X_for_pred_actual <- reverse_min_max_transform_X(X_for_pred,
                                                   X_trans$lower_bound, X_trans$upper_bound)
  # Bind X and predictions
  map_predictions <- rbind(X_for_pred_actual, map_predictions_Y)


  return(list(diff_stat_points = diff_stat_points,
            intervals_list_diff = intervals_list_diff,
            diff_post_cov = stats::cov(diff_stat_points_mat),
            diff_map_par = map_result$map_par,
            diff_map_stat_points = diff_map_stat_points,
            map_predictions = map_predictions,
            posterior_modes = posterior_modes,
            acceptance_rates = mcmc_samples$acceptance_probs[, 1]))


}








# Functions for DGP Method


#' Given posterior samples of DGP method, obtain MAP estimates and Intervals
#'
#' @keywords internal
#' @noRd
dgp_post_processing <- function(sample_t, num_stationary, grid_t, post_prob_gp){

  # Run K means on sample t's
  dgp_kmeans <- stats::kmeans(sample_t, num_stationary)
  centers_order <- order(dgp_kmeans$centers[, 1])

  dgp_stat_points <- list()
  dgp_post_var <- rep(0, num_stationary)

  # Get sampled points for each cluster
  for (i in 1:num_stationary){
    dgp_stat_points[[i]] <- sample_t[dgp_kmeans$cluster == centers_order[i]]
    dgp_post_var[i] <- stats::var(dgp_stat_points[[i]])
  }

  # Get MAP values for each cluster
  dgp_map <- rep(0, num_stationary)
  for(i in 1:num_stationary){
    samp_post_probs <- pracma::interp1(grid_t, post_prob_gp, dgp_stat_points[[i]])
    dgp_map[i] <- dgp_stat_points[[i]][which.max(samp_post_probs)]
  }


  dimnames_list <- list(unlist(lapply(seq(1, length(dgp_stat_points)),
                                      function(x){paste0("S", x)})),
                        c("lower", "upper"))

  # Calculate 90, 95, 99 and Bonferonni HPD Intervals
  intervals_90 <- do.call(rbind, lapply(dgp_stat_points, compute_HPD, prob = 0.90))
  dimnames(intervals_90) <- dimnames_list

  intervals_95 <- do.call(rbind, lapply(dgp_stat_points, compute_HPD, prob = 0.95))
  dimnames(intervals_95) <- dimnames_list

  intervals_99 <- do.call(rbind, lapply(dgp_stat_points, compute_HPD, prob = 0.99))
  dimnames(intervals_99) <- dimnames_list

  bonf_corr <- 1 - 0.05/num_stationary

  intervals_bonf <- do.call(rbind, lapply(dgp_stat_points, compute_HPD, prob = bonf_corr))
  dimnames(intervals_bonf) <- dimnames_list

  intervals_list <- list(intervals_90 = intervals_90,
                         intervals_95 = intervals_95,
                         intervals_99 = intervals_99,
                         intervals_bonf = intervals_bonf)

  return(list(dgp_stat_points = dgp_stat_points,
              intervals_list_dgp = intervals_list,
              dgp_post_var = dgp_post_var,
              dgp_map_stat_points = dgp_map,
              post_prob_gp = post_prob_gp))
}

#' Create difference matrix (Used for kernel of Gaussian process)
#'
#' @keywords internal
#' @noRd
create_H0_matrix <- function(x_vec) {
  return(outer(x_vec, x_vec,
               FUN = function(x1, x2) (x1 - x2)))
}

#' Sample from Stationary Point Univariate Posterior for DGP Model
#'
#' Sample from the stationary point's posterior for the derivative constrained
#' Gaussian process method
#' @param X length n vector, x points (input data)
#' @param Y length n vector, y points (output data)
#' @param num_stationary positive integer, number of stationary points
#' @param total_samples positive integer, number of posterior samples to obtain
#' @param beta_shape1 positive real number, first shape parameter for beta prior
#' @param beta_shape2 positive real number, second shape parameter for beta prior
#' @returns A list with components
#' \describe{
#'   \item{dgp_stat_points}{List of length \code{num_stationary}, posterior samples for
#'   each stationary point}
#'   \item{dgp_post_var}{Vector of length \code{num_stationary}, posterior variance for
#'   each stationary point}
#'   \item{dgp_map_stat_points}{Vector of length \code{num_stationary}, MAP of stationary
#'   points}
#'   \item{intervals_list_dgp}{List of matrices, 90, 95, 99 and Bonferonni Corrected Intervals for
#'   each stationary point}
#'   \item{post_prob_gp}{Length 5000 vector from of univariate posterior applied to
#'    min(X) to max(X)}
#' }
#' @export
posterior_stat_points_dgp <- function(X, Y, num_stationary, total_samples,
                                      beta_shape1 = 1, beta_shape2 = 1) {

  # Obtain Gaussian process parameters
  n <- length(X)
  H0_mat <- create_H0_matrix(X)

  EB_gp <- Rsolnp::solnp(pars = c(1, 1, 1), fun = dgp::log_mar_lik_gp,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = Y, H0 = H0_mat)$par


  # Gaussian Process parameters
  lambda_gp <- EB_gp[1]^2 / (n * EB_gp[2]^2)
  Kff_gp <- dgp::se_ker(H0 = H0_mat, tau = 1, h = EB_gp[3])
  A_gp <- Kff_gp + diag((n * lambda_gp), n)

  # Obtain log Posterior of stationary point parameter
  grid_t <- seq(min(X), max(X), length.out = 5000)
  log_post_gp <- rep(0, 5000)

  for (j in 1:5000){

    log_post_gp[j] <- dgp::log_post_t_theory(t = grid_t[j], y = Y, x = X,
                                        Kff = Kff_gp, A = A_gp, lambda = lambda_gp,
                                        h = EB_gp[3], sig2 = EB_gp[1]^2,
                                        shape1 = beta_shape1, shape2 = beta_shape2,
                                        a = min(X), b = max(X))
  }

  # Obtain unnormalized posterior values at grid_t
  post_prob_gp <- exp(log_post_gp - max(log_post_gp))

  # Estimate CDF values of posterior
  cdf_post_prob_gp <- pracma::cumtrapz(grid_t, post_prob_gp)
  cdf_post_prob_gp <- cdf_post_prob_gp/cdf_post_prob_gp[5000]

  # Obtain Samples by Inverse transform method
  sample_t <- pracma::interp1(cdf_post_prob_gp[, 1], grid_t, stats::runif(total_samples))

  # Post process samples
  return(dgp_post_processing(sample_t, num_stationary, grid_t, post_prob_gp))

}

#' Apply LogSumExp function to a vector
#'
#' @keywords internal
#' @noRd
log_sum_exp <- function(vec){

  m_val <- max(vec)
  exp_vec <- exp(vec - m_val)

  return(m_val + log(sum(exp_vec)))
}

#' Apply SoftMax function to a vector
#'
#' @keywords internal
#' @noRd
softmax <- function(vec) {

  m_val <- max(vec)
  exp_vec <- exp(vec - m_val)

  return(exp_vec/sum(exp_vec))
}

#' Check if lambdas height vector is valid (peak, valley, peak, valley, etc)
#'
#' @keywords internal
#' @noRd
check_fluctuating_lambdas <- function(lambdas, first_direction) {

  lam_len <- length(lambdas)

  # True Sign of the difference between the consecutive lambdas
  if (first_direction == 1) {
    true_diff_vec <- rep(c(1, -1), length.out = lam_len - 1)
  } else{
    true_diff_vec <- rep(c(-1, 1), length.out = lam_len - 1)
  }

  # See if sign of difference between lambdas matches what it should be
  diff_vec <- sign(lambdas[2:lam_len] - lambdas[1:(lam_len - 1)])

  if(sum(diff_vec == true_diff_vec) == lam_len - 1) {
    return(TRUE)
  } else{
    return(FALSE)
  }
}

#' Return likelihood values for (input_X, input_Y) for specific (betas, lambdas, sigma and b_vec) parameters
#'
#' @keywords internal
#' @noRd
like_value <- function(betas, lambdas, sigma, b_vec, input_X,
                       input_Y, first_direction, interpolation,
                       weight_vector = NULL, return_all = FALSE) {

  # If norm > pi, return infinity (invalid point)
  if (sum(betas^2) >(pi^2) + 1e-05) {
    return(-Inf)
  }

  # Check if lambdas are fluctuating
  if(!check_fluctuating_lambdas(lambdas, first_direction)){
    return(-Inf)
  }

  # Checks is sigma is positive
  if (sigma <= 0) {
    return(-Inf)
  }

  # Use an equally spaced b_vec if none given
  if(is.null(b_vec)){
    b_vec <- seq(0, 1, length.out = length(lambdas))
  }

  # Return diffeomorphism likelihood values for all (X, Y)
  all_log <- diffeo_all(input_X, input_Y, betas, b_vec,
                        lambdas, c(sigma),  interpolation = interpolation,
                        transform_X = FALSE)

  # Return all likelihood values, the sum or
  # Weighted Sum of Weighted Bayesian Bootstrap
  if(return_all){
    return(as.vector(all_log$log_likelihood))
  } else{
    if (is.null(weight_vector)){
      return(sum(all_log$log_likelihood))
    } else {
      return(sum(weight_vector * as.vector(all_log$log_likelihood)))
    }
  }

}

#' Given unconstrained parameterization, transform back to original parameterization of
#' lambdas and sigmas. beta_vec will remain on continuous scale where ||beta|| < pi
#'
#' @keywords internal
#' @noRd
process_post_samples <- function(param_vec, num_betas, first_direction,
                                 zero_is_zero = FALSE){

  total_len <- length(param_vec)

  # Get betas
  betas <- param_vec[1:num_betas]

  # Get the untransformed lambda variables
  lambdas_raw <- param_vec[(num_betas + 1):(total_len - 1)]

  # Depending on if the first lambda value is 0
  if (zero_is_zero){
    correction <- 0
  } else{correction <- -1}

  if(first_direction == 1){
    switch_vec <- rep(c(1, -1),
                      length.out = length(lambdas_raw) + correction)
  } else {
    switch_vec <- rep(c(-1, 1),
                      length.out = length(lambdas_raw) + correction)
  }

  # Transform lambda values into alternating heights
  if(zero_is_zero){
    lambdas <- c(0, cumsum(switch_vec * exp(lambdas_raw)))
  } else{
    lambda_1 <- lambdas_raw[1]
    adding_to <- switch_vec * exp(lambdas_raw[2:length(lambdas_raw)])
    lambdas <- cumsum(c(lambda_1, adding_to))

  }

  # Turn sigma positive
  sigma <- exp(param_vec[total_len])

  return(list(betas = betas, lambdas = lambdas,
              sigma = sigma))

}

#' Given unconstrained parameterization, transform back to original parameterization. Turns
#' result from process_post_samples into vector
#'
#' @keywords internal
#' @noRd
process_post_samples_as_vec <- function(param_vec, num_betas, first_direction,
                                        zero_is_zero = FALSE) {

  process_list <- process_post_samples(param_vec, num_betas, first_direction,
   zero_is_zero)

  return(unlist(process_list, use.names = FALSE))

}


#' Given unconstained parameters (||betas|| < pi and (-infty, infty) for ) find
#' log likelihood value for entire dataset (input_X, input_Y)
#'
#' @keywords internal
#' @noRd
like_value_w_processing <- function(param_vec, input_X, input_Y, num_betas,
                                     first_direction, zero_is_zero,
                                    interpolation, b_vec,
                                     weight_vector = NULL,
                                    return_all = FALSE) {
  # Transform parameters
  process_values <- process_post_samples(param_vec, num_betas,
                                         first_direction,
                                         zero_is_zero)

  # Calculate log likelihood for dataset
  log_like <- like_value(process_values$betas, process_values$lambdas,
                         process_values$sigma, b_vec,
                         input_X, input_Y, first_direction, interpolation,
                         weight_vector, return_all)

  return(log_like)
}


#' Given unconstrained parameters (||betas|| < pi and (-infty, infty)) find log likelihood +
#' log prior for entire dataset. Prior is assumed to be independent multivariate normal
#'
#' @keywords internal
#' @noRd
bayes_value <- function(param_vec, input_X, input_Y, num_betas,
                        first_direction, zero_is_zero,
                        interpolation, b_vec, prior_mean, prior_sd,
                        weight_vector = NULL, w_0 = 1
                       ) {

  # Get log likelihood
  log_like <- like_value_w_processing(param_vec, input_X, input_Y, num_betas,
                                      first_direction, zero_is_zero,
                                      interpolation, b_vec)


  # Get prior likelihood
  if ((is.null(prior_sd)) | (is.null(prior_mean))){
    return(log_like)
  } else{

    log_prior <- sum(w_0 * stats::dnorm(param_vec, mean = prior_mean,
                           sd = prior_sd, log = TRUE))
    return(log_like + log_prior)
  }
}

#' Get negative of bayes_value, used for minimization and finding posterior modes
#'
#' @keywords internal
#' @noRd
neg_bayes_value <- function(param_vec, input_X, input_Y, num_betas,
                            first_direction, zero_is_zero,
                            interpolation, b_vec, prior_mean, prior_sd,
                            weight_vector = NULL, w_0 = 1){

  return(-1 * bayes_value(param_vec, input_X, input_Y, num_betas,
                          first_direction, zero_is_zero,
                          interpolation, b_vec, prior_mean, prior_sd,
                          weight_vector, w_0))
}


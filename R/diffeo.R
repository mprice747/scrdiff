# Functions for applying diffeomorphisms to data

#' Applies interpolation to transformed X via diffeomorphism
#'
#' @keywords internal
#' @noRd
apply_interpolation_per_row <- function(interp_vec, b_length,
                                        lambda_length,
                                        interpolation) {

  # Separates the b-points, lambda values and the
  # diffeomorphism values (points to evaluate interpolation)
  b_vec <- interp_vec[1:b_length]
  lambda_vec <- interp_vec[(b_length + 1):lambda_length]
  diffeomorphism_vec <- interp_vec[(lambda_length + 1):length(interp_vec)]

  # Apply interpolation
  interp_results_vec <- pracma::interp1(b_vec, lambda_vec, diffeomorphism_vec,
                                method = interpolation)

  return(interp_results_vec)

}


#' At predictions (estimation of E(Y|X)), get log likelihood at points Y
#'
#' @keywords internal
#' @noRd
diffeo_log_likelihood <- function(predictions, Y, Sigma) {

  log_likelihood <- stats::dnorm(predictions, mean = Y, sd = Sigma, log = TRUE)

  return(log_likelihood)

}

#' Given X and different Beta (Diffeomorphism Parameters), BPoints (Input Points for Spline Interpolation)
#' and Lambda (Height Vectors) parameters, output model predictions at points X. Option to
#' output diffeomorphism transformation of X
#'
#' @keywords internal
#' @noRd
diffeo_predictions <- function(X, Beta, BPoints, Lambda,
                               interpolation = 'cubic', transform_X = TRUE,
                               return_diffeo = TRUE) {


  # Turn X, Y and sigma into  matrix
  if (is.vector(X)){
    X <- as.matrix(X)
  }
  if ((transform_X) | (min(X) < 0) | (max(X) > 1)) {
    X <- min_max_transform_X(X)$new_X
  }

  # Turn Beta, B Points and Lambda Into Matrix
  Beta <- turn_1d_into_matrix(Beta)
  BPoints <- turn_1d_into_matrix(BPoints)
  Lambda <- turn_1d_into_matrix(Lambda)

  # Get result of diffeomorphism of transformed X space
  diffeomorphism <- gamma_function(X, Beta)
  diffeomorphism[diffeomorphism <= 0] <- 0
  diffeomorphism[diffeomorphism >= 1] <- 1

  # Apply cubic interpolation to transformed diffeomorphism of X space
  b_length <- ncol(BPoints)
  lambda_length <- b_length + ncol(Lambda)

  # Combine BPoints, Lambda and Diffeomorphism into 1 big matrix
  interp_matrix <- cbind(BPoints, Lambda, t(diffeomorphism))

  # Apply interpolations for each row to get predictions for each different
  # set of parameters
  predictions <- apply(interp_matrix, 1, apply_interpolation_per_row,
                       b_length = b_length,
                       lambda_length = lambda_length,
                       interpolation = interpolation)

  if (return_diffeo){
    return(list(predictions = predictions, diffeomorphism = diffeomorphism))
  } else{
    return(predictions)
  }
}

#' At points X, get predictions E(Y|X), the diffeomorphism transformation of X and
#' the log likelihood values for different sets of parameters
#'
#' @keywords internal
#' @noRd
diffeo_all <- function(X, Y, Beta, BPoints, Lambda, Sigma,
                                  interpolation = 'cubic', transform_X = TRUE,
                                  check_compat = TRUE) {


  # Get predictions, and then the log likelihood
  predictions_diffeo <- diffeo_predictions(X, Beta, BPoints, Lambda,
         interpolation = interpolation, transform_X = TRUE, return_diffeo = TRUE)

  log_likelihood <- diffeo_log_likelihood(predictions_diffeo$predictions,
                                          Y, Sigma)
  # List of predictions, diffeomorphism transformations and log likelihood
  return(list(predictions = predictions_diffeo$predictions,
              diffeomorphism = predictions_diffeo$diffeomorphism,
              log_likelihood = log_likelihood))
}


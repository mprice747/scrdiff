# Functions for applying diffeormorphism transformation using diffeomorphic
# parameterization

#' SQRT(2) constant
#' @noRd
SQRT_2 <- sqrt(2)

#' Check if input is matrix. If not, turn into 1 column matrix
#'
#' @keywords internal
#' @noRd
check_matrix <- function(X) {

  if (!is.matrix(X)) {
    X <- matrix(X, ncol = 1)
  }

  return(X)
}

#' Check if input is matrix. If not, turn into 1 row matrix
#'
#' @keywords internal
#' @noRd
turn_1d_into_matrix <- function(X) {

  if (is.null(dim(X))) {
    X <- matrix(X, nrow = 1)
    return(X)
  }
  return(X)
}

#' Inputs X into m different cosine bases of length p, represented as Beta, mxp matrix
#'
#' @keywords internal
#' @noRd
cosine_basis <- function(X, Beta) {
  # Apply cosine, then the dot product with Beta parameters
  Cos_Mat <- cos(pi * (X %*% matrix(seq(1, ncol(Beta)), nrow = 1)))

  result <- tcrossprod(Cos_Mat, (SQRT_2 * Beta))

  return(result)
}

#' Calculates integral from 0 to X for m different cosine bases of length p (Beta)
#'
#' @keywords internal
#' @noRd
int_cosine_basis <- function(X, Beta) {

  # Integral of cos(pi*j*x) is sin(pi * j & x)/(pi * J)
  j_vec <- seq(1, ncol(Beta))

  Sin_Mat <- sin(pi * (X %*% matrix(j_vec, nrow = 1)))

  result <- Sin_Mat %*% (t((SQRT_2 * Beta))/(j_vec * pi))

  return(result)
}

#' Calculates int_{0}^{x}cos^2(y * pi * j)dy
#'
#' @keywords internal
#' @noRd
cos_2_integral <- function(x, j) {

  two_pi_j <- 2 * pi * j
  result <- sin(two_pi_j * x)/(2 * two_pi_j) + x/2

  return(result)

}
#' Calculates int_{0}^{x} cos(y * pi * j1)cos(y * pi * j2)dy where j1 != j2
#'
#' @keywords internal
#' @noRd
cos_cos_integral <- function(x, j1, j2){

  pi_j1_j2 <- pi * (j1 + j2)
  pi_j2_sub_j1 <- pi * (j2 - j1)

  result <- sin(pi_j1_j2 * x)/(2 * pi_j1_j2) +
    sin(pi_j2_sub_j1 * x)/(2 * pi_j2_sub_j1)

  return(result)
}

#' Given vec, output 2 x choose(length(vec), 2) matrix containing every combination of elements within vec
#'
#' @keywords internal
#' @noRd
different_beta_product <- function(vec) {

  vec_comb <- utils::combn(vec, 2)

  return(vec_comb[1, ] * vec_comb[2, ])
}

#' Calculate integral from 0 to X of m cosine basis of length p squared,
#' with Beta as weights (mxp matrix)
#'
#' @keywords internal
#' @noRd
int_cosine_basis_2 <- function(X, Beta) {

  # Produce every possible pair of position
  j_seq <- seq(1, ncol(Beta))
  j_comb <- utils::combn(j_seq, 2)


  # Get cos^2 squared integral of different basis functions
  cos_cos_ints <- turn_1d_into_matrix(apply(X, 1, cos_cos_integral,
                                            j1 = j_comb[1, ], j2 = j_comb[2, ]))
  diff_betas_mult <- turn_1d_into_matrix(apply(Beta, 1, different_beta_product))
  cos_cos_all <- 4 * crossprod(cos_cos_ints, diff_betas_mult)

  # Get cos^2 squared integral of same basis functions
  cos_2_ints <- turn_1d_into_matrix(apply(X, 1, cos_2_integral, j = j_seq))

  # Add all to return result
  cos_2_all <- 2 * t(Beta^2 %*% cos_2_ints)
  result <- cos_cos_all + cos_2_all

  return(result)
}

#' Applies m Gamma function (diffeomorphism), with each being represented with
#' length p cosine basis (Beta (mxp matrix))
#'
#' @keywords internal
#' @noRd
gamma_function <- function(X, Beta) {

  # Find integral of cosine basis and cosine basis squared
  Beta_norms <- sqrt(rowSums(Beta^2))

  int_1 <- t(int_cosine_basis(X, Beta))

  int_2 <- t(int_cosine_basis_2(X, Beta))

  # cosine and sines of norms in exponential map
  cos_beta <- cos(Beta_norms)
  sin_beta_div <- sin(Beta_norms)/Beta_norms
  cos_beta_X <- outer(cos(Beta_norms)^2, X[, 1])

  # Full integral of square exponential map
  exp_map <- cos_beta_X + ((2 * cos_beta *
                              sin_beta_div) * int_1) + sin_beta_div^2 * int_2

  return(t(exp_map))

}


#' Transform X to interval in between 0 and 1 (Required for diffeomorphism)
#'
#' @keywords internal
#' @noRd
min_max_transform_X <- function(X) {

  lower_bound <- min(X)
  upper_bound <- max(X)

  # Min-max scaling
  new_X <- (X - lower_bound)/(upper_bound - lower_bound)

  return(list(new_X = new_X, lower_bound = lower_bound, upper_bound = upper_bound))

}

#' Transform min_max_transform_X output back to original X
#'
#' @keywords internal
#' @noRd
reverse_min_max_transform_X <- function(new_X, lower_bound, upper_bound){
  return(new_X * (upper_bound - lower_bound) + lower_bound)
}





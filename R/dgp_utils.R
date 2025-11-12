# Adapted from dgp (https://github.com/chenghanyustats/dgp)
# Original author: Cheng-Han Yu
# License: GPL-3

#Importing from Github caused issues, this is copy pasted directly from the github link

#' Log dgbeta
#' @noRd
log_dgbeta <- function(x, shape1, shape2, a, b) {
  if (shape1 == 1 && shape2 == 1) {
    if ((b - a) == 1) {
      return(0)
    } else {
      return(log(1 / (b - a)))
    }
  } else {
    return((shape1 - 1) * log(x - a) + (shape2 - 1) * log(b - x) -
             lbeta(shape1, shape2) - (shape1 + shape2 - 1) * log(b - a))
  }
}

#' Kernel function
#' @noRd
se_ker <- function(x1, x2, H0, tau = 1, h = 1) {
  if (missing(x1) || missing(x2)) {
    if(!missing(H0)) {
      ## H0 can be a matrix from outer operation
      return(tau ^ 2 * exp(- (H0 / h) ^ 2 / 2))
    } else {
      stop("must provide either (x1, x2) or H0")
    }
  } else {
    return(tau ^ 2 * exp(- ((x1 - x2) / h) ^ 2 / 2))
  }
}


#' Kernel derivative
#' @noRd
se_ker_der1 <- function(x1, x2, H0, tau = 1, h = 1) {
  if (missing(x1) || missing(x2)) {
    if(!missing(H0)) {
      ### H0 can be a matrix
      return(se_ker(H0, tau, h) * (H0) / (h ^ 2))
    } else {
      stop("must provide either (x1, x2) or H0")
    }
  } else {
    return(se_ker(x1, x2, tau = tau, h = h) * (x1 - x2) / (h ^ 2))
  }
}




#' Compute covariance derivative
#' @noRd
computeCovDer1 <- function(idx1, idx2, ker_der1_fcn = se_ker_der1, ...) {
  # take derivative w.r.t. the 1st argument
  return(outer(idx1, idx2, function(x1, x2) ker_der1_fcn(x1, x2, ...)))
}

#' Computing log posterior for DGP method
#' @noRd
log_post_t_theory <- function(t, y, x, Kff, A, lambda, h, sig2,
                              shape1 = 1, shape2 = 1, a = 0, b = 2) {
  n <- length(x)
  tau2 <- sig2 / (n * lambda)
  Kdf <- computeCovDer1(idx1 = t, idx2 = x,
                        tau = 1, h = h)

  a_mat <- t(Kdf) * h
  mu_t <- t(a_mat) %*% solve(A, y) / h
  dd <- 1 - emulator::quad.form.inv(A, a_mat)
  sig2_t <- tau2 * dd / (h ^ 2)
  cholA <- chol(A)
  log_det_A <- 2 * sum(log(diag(cholA)))
  log_C <- (-1 / 2) * (n * log(2 * pi * tau2) + log_det_A +
                         emulator::quad.form.inv(A, y) / tau2 - log(tau2) +
                         2 * log(h))
  # log_C <- (-n / 2) * log(2 * pi * tau2) -
  #     (1 / 2) * log(det(A)) - quad.form.inv(A, y) / (2 * tau2) +
  #     (1 / 2) * log(tau2) - log(h)
  log_den_t <- (-1 / 2) * (log(sig2_t) + crossprod(mu_t) / sig2_t)
  log_prior <- log_dgbeta(t, shape1, shape2, a, b)

  return(as.vector(log_C + log_den_t + log_prior))
}

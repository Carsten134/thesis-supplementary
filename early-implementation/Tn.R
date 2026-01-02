library(ravetools)
library(rmutil)
library(MASS)

gridProcess <- function(N, M, K, padding = nrow(sigma) %/% 2) {
  N_tilde <- N + 2 * padding
  M_tilde <- M + 2 * padding
  n <- N_tilde * M_tilde
  
  eps <- matrix(rnorm(n),
                nrow = M_tilde,
                ncol = N_tilde)
  
  x <- convolve_image(eps, K)
  # cut padding
  return(x[(padding+1):(M+padding),(padding+1):(N+padding)])
}

# implementation of phi_n
I <- function(x){
  N <- nrow(x)
  M <- ncol(x)
  # use twodim fft
  temp <- apply(x, fft, MARGIN = 2)
  temp <- apply(temp, fft, MARGIN = 1)
  
  # apply scaling
  res <- (1/(pi^2*N*M))*Mod(temp)^2
  
  # fft shift such that 0 is at the center and nyquist is at the borders 
  return(res[
    c(
      ((N %/% 2)+1):N,
      1:(N %/% 2)),
    c(
      ((M %/% 2)+1):M,
      1:(M %/% 2))]
  )
}


tr_k <- function(u, h = 1) ifelse(u/h <= 1 & u / h >= -1, (35 / 32) * (1-(u/h)^2)^3, 0)

I_smooth <- function(I, kernel, hr = 1, hc = 1) {
  N <- nrow(I)
  M <- ncol(I)
  
  omega_M <- 2*pi*((-((M - 2) %/% 2)):(M %/% 2))/M
  omega_N <- 2*pi*((-((N - 2) %/% 2)):(N %/% 2))/N
  function (omega_1, omega_2) {
    K_M <- kernel(omega_M - omega_1, hr)
    K_N <- kernel(omega_N - omega_2, hc)
    
    weights <- K_N %*% t(K_M)
    weighted <- sum(weights * I)
    
    normalized <- (1/(N*M*hr*hc))*weighted
    return(normalized)
  }
}

Tn_f <- function(x,y, kernel, ...) {
  # cache periodogram and fourier frequencies
  I_x <- I(x)
  I_y <- I(y)
  I_tilde <- .5*(I_x + I_y)
  
  I_x_smooth <- I_smooth(I_x - I_tilde, kernel, ...)
  I_y_smooth <- I_smooth(I_y - I_tilde, kernel, ...)
  
  function(omega_1, omega_2) {
    return(
      I_x_smooth(omega_1, omega_2)^2 + I_y_smooth(omega_1, omega_2)^2
    )
  }
}

Tn <- function(x, y, kernel = tr_k, ...) {
  l2 <- Tn_f(x, y, kernel, ...)
  int2(l2, c(-pi, -pi), c(pi, pi))
}

Tn_star_f <- function(x, y, kernel = tr_k, ...) {
  I_x <- I(x)
  I_y <- I(y)
  
  I_tilde <- .5*(I_x + I_y)
  
  # permuting entries
  N <- nrow(I_x)
  M <- ncol(I_x)
  perm <- matrix(runif(N*M), ncol = N)
  perm <- as.numeric(perm > 0.5)
  perm_n <- as.numeric(!perm)
  I_x_rand <- I_x*perm + I_y*perm_n
  I_y_rand <- I_x*perm_n + I_y*perm
  
  I_x_diff <- I_smooth(I_x_rand - I_tilde, kernel, ...)
  I_y_diff <- I_smooth(I_y_rand - I_tilde, kernel, ...)
  
  function(omega_1, omega_2) {
    return(
      I_x_diff(omega_1, omega_2)^2 + I_y_diff(omega_1, omega_2)^2
    )
  }
}

Tn_star <- function(x, y, kernel = tr_k, ...) {
  l2_rand <- Tn_star_f(x, y, kernel, ...)
  int2(l2_rand, c(-pi, -pi), c(pi, pi))
}

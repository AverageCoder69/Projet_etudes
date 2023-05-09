params<-list(
  alpha=c( 0.5, 0.4, 0.1
  ),
  sigma=c( 1, 2,0.2),
  a = 1,
  b = 2,
  v = 0.5,
  m = 5,
  r = 2,
  d = seq(0, 1, length = 95))
alpha <- params$alpha
a <- params$a
b <- params$b
sigma1 <- params$sigma[1]
sigma2 <- params$sigma[2]
sigma3 <- params$sigma[3]
v <- params$v
m <- params$m
r <- params$r
d <- params$d
Alpha<-2/( params$a+params$b)


C <- function(x) {
 C1=(x<=-params$a)*(x>=0)*x
C2=(x>0)*(x<=params$b)*x
return(Alpha*(C1+C2))
}


f2 <- function(x, sigma1, sigma2, v) {
  exp(-x^2 / (2 * sigma1^2)) * (1 + v * cos(2 * pi * x / sigma2))
}


integral <- integrate(f2, lower = -Inf, upper = Inf, sigma1 = sigma1, sigma2 = sigma2, v = v)$value


beta <- 1 / integral


V <- function(x, sigma1, sigma2, v) {
  beta * exp(-x^2 / (2 * sigma1^2)) * (1 + v * cos(2 * pi * x / sigma2))
}


f3 <- function(x, m, sigma3, r) {
  ifelse(x >= m - sigma3 & x <= m + sigma3, (sigma3 - abs(x - m)) / x^r, 0)
}



integral <- integrate(f3, lower = m - sigma3, upper = m + sigma3, m = m, sigma3 = sigma3, r = r)$value


gamma <- 1 / integral


B <- function(x, m, sigma3, r) {
  gamma * ifelse(x >= m - sigma3 & x <= m + sigma3, (sigma3 - abs(x - m)) / x^r, 0)
}


fS <- function(x, a, b, alphaC, alphaV, alphaB, sigma1, sigma2, v, m, sigma3, r) {
  C_x <- C(x)
  V_x <- V(x, sigma1, sigma2, v)
  B_x <- B(x, m, sigma3, r)
  
  alphaC * C_x + alphaV * V_x + alphaB * B_x
}


rs <- function(n, params) {
  a <- params$a
  b <- params$b
  alphaC <- params$alpha[1]
  alphaV <- params$alpha[2]
  alphaB <- params$alpha[3]
  sigma1 <- params$sigma[1]
  sigma2 <- params$sigma[2]
  v <- params$v
  m <- params$m
  sigma3 <- params$sigma[3]
  r <- params$r
  
  XX <- seq(-a, m + sigma3, a / 100)
  envelope_height <- max(fS(XX, a, b, alphaC, alphaV, alphaB, sigma1, sigma2, v, m, sigma3, r))
  
  n_attempts <- ceiling(n * 1.5)
  x_samples <- runif(n_attempts, min = -a, max = max(b, m + sigma3))
  y_samples <- runif(n_attempts)
  
  ps_density <- fS(x_samples, a, b, alphaC, alphaV, alphaB, sigma1, sigma2, v, m, sigma3, r)
  H <- y_samples * x_samples / (a + max(b, m + sigma3))
  
  accepted_samples <- x_samples[H <= ps_density]
  
  if (length(accepted_samples) >= n) {
    return(accepted_samples[1:n])
  } else {
    return(c(accepted_samples, rs(n - length(accepted_samples), params)))
  }
}

esp.TR <- function(n.simul, params) {
  const_1_96 <- 1.96
  sqrt_n_simul <- sqrt(n.simul)
  
  combined <- function(d, params) {
   duree_vie_individus <- rmultinom(1, 1, d)
    position <- match(1, duree_vie_individus)
    lifetime <- position 
     if (lifetime >= 40) {
      lifetime1 <- 39
      T_vector <- numeric(lifetime1)  
      last_value <- 1  
      D <- rs(lifetime1, params)
      T_vector[1] <- last_value * exp(D[1])

      for (i in 2:lifetime1) {
      T_vector[i] <- T_vector[i - 1] * exp(D[i])
}
      derniere <- T_vector[i]
      T <- sum(T_vector)
    } else {
      T_vector <- numeric(lifetime)  
      last_value <- 1  
      D <- rs(lifetime, params)
      T_vector[1] <- last_value * exp(D[1])

      for (i in 2:lifetime) {
      T_vector[i] <- T_vector[i - 1] * exp(D[i])
}
      T <- sum(T_vector)}


    afterlife <- (lifetime - 39)

    if (afterlife > 0) {
      annual_pension <- 0.75 * derniere
      R <- (afterlife * annual_pension)
    } else {
      R <- 0
    }

    return(c(T = T, R = R))

  }

  TR_samples <- sapply(1:n.simul, function(x) combined(d, params))

  std_dev_T <- sd(TR_samples[1, ])
  std_dev_R <- sd(TR_samples[2, ])
  ET <- mean(TR_samples[1, ])
  ER <- mean(TR_samples[2, ])
  demi_largeur_ET <- const_1_96 * (std_dev_T / sqrt_n_simul)
  demi_largeur_ER <- const_1_96 * (std_dev_R / sqrt_n_simul)
  return((list(ET = ET, demi.largeur.ET = demi_largeur_ET, ER = ER, demi.largeur.ER = demi_largeur_ER)))
}

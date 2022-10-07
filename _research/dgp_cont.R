simdata <- function(n = 1000) {
  w1 <- rbinom(n, 1, 0.6)
  w2 <- rbinom(n, 1, 0.3)
  w3 <- rbinom(n, 1, 0.4)
  a <- rbinom(n, 1, 0.5)
  z <- rbinom(n, 1, plogis((-log(1.3)*(w1 + w2 + w3) / 3) + 2*a - 1))
  m <- -log(1.1)*w3 + 2*z - 0.9 + rnorm(n, 2)
  y <- (-log(1.3)*(w1 + w2 + w3) / 3) + z + m + rnorm(n, 5, 2)
  
  data.frame(w1 = w1, 
             w2 = w2, 
             w3 = w3, 
             a = a, 
             z = z, 
             m = m, 
             y = y)
}

truth_dgp_cont <- function() {
  tmp <- expand.grid(w1 = c(1, 0), w2 = c(1, 0), w3 = c(1, 0))
  
  prob_w <- vector("numeric", nrow(tmp))
  for (i in 1:nrow(tmp)) {
    w1 <- tmp[i, "w1"] * 0.6 + (1 - tmp[i, "w1"]) * 0.4
    w2 <- tmp[i, "w2"] * 0.3 + (1 - tmp[i, "w2"]) * 0.7
    w3 <- tmp[i, "w3"] * 0.4 + (1 - tmp[i, "w3"]) * 0.6
    prob_w[i] <- w1 * w2 * w3
  }
  
  gZ <- function(x, a) {
    plogis((-log(1.3)*rowSums(x[, 1:3]) / 3) + 2*a - 1)
  }
  
  beta_0 <- 5
  beta_W <- rep((1/3) * -log(1.3), 3)
  beta_Z <- 1
  beta_M <- 1
  
  alpha_0 <- 1.1
  alpha_W <- -log(1.1)
  alpha_Z <- 2
  
  theta_11_10 <- (beta_0 + as.matrix(tmp[, 1:3]) %*% beta_W + beta_Z + beta_M * 
                    (alpha_0 + alpha_W*tmp$w3 + alpha_Z)) * gZ(tmp, 0)
  
  theta_10_10 <- (beta_0 + as.matrix(tmp[, 1:3]) %*% beta_W + beta_Z + beta_M * 
                    (alpha_0 + alpha_W*tmp$w3)) * (gZ(tmp, 1) - gZ(tmp, 0))
  
  theta_00_10 <- (beta_0 + as.matrix(tmp[, 1:3]) %*% beta_W + beta_M * 
                    (alpha_0 + alpha_W*tmp$w3)) * (1 - gZ(tmp, 1))
  
  theta_10 <- weighted.mean(theta_11_10 + theta_10_10 + theta_00_10, prob_w)
  
  theta_11_00 <- (beta_0 + as.matrix(tmp[, 1:3]) %*% beta_W + beta_Z + beta_M * 
                    (alpha_0 + alpha_W*tmp$w3 + alpha_Z)) * gZ(tmp, 0)
  
  theta_10_00 <- 0
  
  theta_00_00 <- (beta_0 + as.matrix(tmp[, 1:3]) %*% beta_W + beta_M * 
                    (alpha_0 + alpha_W*tmp$w3)) * (1 - gZ(tmp, 0))
  
  theta_00 <- weighted.mean(theta_11_00 + theta_10_00 + theta_00_00, prob_w)
  
  theta_11_11 <- (beta_0 + as.matrix(tmp[, 1:3]) %*% beta_W + beta_Z + beta_M * 
                    (alpha_0 + alpha_W*tmp$w3 + alpha_Z)) * gZ(tmp, 1)
  
  theta_10_11 <- 0
  
  theta_00_11 <- (beta_0 + as.matrix(tmp[, 1:3]) %*% beta_W + beta_M * 
                    (alpha_0 + alpha_W*tmp$w3)) * (1 - gZ(tmp, 1))
  
  theta_11 <- weighted.mean(theta_11_11 + theta_10_11 + theta_00_11, prob_w)
  
  c("nie" = theta_11 - theta_10, 
    "nde" = theta_10 - theta_00)
}

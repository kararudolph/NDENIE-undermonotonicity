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

my_inst <- function(m, z, w) {
  (-log(1.3)*(rowSums(w)) / 3) + z + m + rnorm(nrow(w), 5, 2)
}

pm <- function(m, z, w) {
  dnorm(m, -log(1.1)*w[, "w3"] + 2*z - 0.9, 2)
}

g <- function(a, w) {
  pscore <- .5
  a * pscore + (1 - a) * (1 - pscore)
}

pz <- function(z, a, w) {
  prob1 <- plogis((-log(1.3)*(rowSums(w)) / 3) + 2*a - 1)
  z * prob1 + (1 - z) * (1 - prob1)
}

pzw <- function(z, w) {
  pz(z, 1, w) * g(1, w) + pz(z, 0, w) * g(0, w)
}

e <- function(a, z, w) {
  pz(z, a, w) * g(a, w) / pzw(z, w)
}

pmaw <- function(m, a, w) {
  pm(m, 1, w) * pz(1, a, w) + pm(m, 0, w) * pz(0, a, w)
}

r <- function(z, a, m, w) {
  pm(m, z, w) * pz(z, a, w) / pmaw(m, a, w)
}

truth <- function() {
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


eic <- function(data, aprime, astar) {
  a <- data$a
  z <- data$z
  m <- data$m
  y <- data$y
  z <- data$z
  w <- data[, paste0("w", 1:3)]
  
  beta_0 <- 5
  beta_W <- rep((1/3) * -log(1.3), 3)
  beta_Z <- 1
  beta_M <- 1
  
  alpha_0 <- 1.1
  alpha_W <- -log(1.1)
  alpha_Z <- 2
  
  if (aprime == 1 && astar == 0) {
    muintm_11 <- (beta_0 + as.matrix(w) %*% beta_W + beta_Z + beta_M * 
                    (alpha_0 + alpha_W*data$w3 + alpha_Z))
    
    muintm_10 <- (beta_0 + as.matrix(w) %*% beta_W + beta_Z + beta_M * 
                    (alpha_0 + alpha_W*data$w3))
    
    muintm_00 <- (beta_0 + as.matrix(w) %*% beta_W + beta_M * 
                    (alpha_0 + alpha_W*data$w3))
  }

  if (aprime == 0 && astar == 0) {
    muintm_11 <- (beta_0 + as.matrix(w) %*% beta_W + beta_Z + beta_M * 
                    (alpha_0 + alpha_W*data$w3 + alpha_Z))
    
    muintm_10 <- 0
      
    muintm_00 <- (beta_0 + as.matrix(w) %*% beta_W + beta_M * 
                    (alpha_0 + alpha_W*data$w3))
  }

  if (aprime == 1 && astar == 1) {
    muintm_11 <- (beta_0 + as.matrix(w) %*% beta_W + beta_Z + beta_M * 
                    (alpha_0 + alpha_W*data$w3 + alpha_Z))
    
    muintm_10 <- 0
      
    muintm_00 <- (beta_0 + as.matrix(w) %*% beta_W + beta_M * 
                    (alpha_0 + alpha_W*data$w3))
  }
  
  eif_y_11 <- as.numeric(a == aprime & z == 1) / g(astar, w) * e(astar, 1, w) / e(aprime, 1, w) * (y - my_inst(m, z, w))
  eif_y_10 <- as.numeric(a == aprime & z == 1) / (g(astar, w) * pz(0, astar, w)) * e(astar, 1, w) / e(aprime, 1, w) * r(0, astar, m, w) / r(1, astar, m, w) *  (pz(1, aprime, w) - pz(1, astar, w)) * (y - my_inst(m, z, w))
  eif_y_00 <- as.numeric(a == aprime & z == 0) / g(aprime, w) * e(astar, 0, w) / e(aprime, 0, w) * pz(0, aprime, w) / pz(0, astar, w) * (y - my_inst(m, z, w))
  
  eif_m_11 <- as.numeric(a == astar & z == 1) / g(astar, w) * (my_inst(m, 1, w) - muintm_11)
  eif_m_10 <- as.numeric(a == astar & z == 0) / (g(astar, w) * pz(0, astar, w)) * (pz(1, aprime, w) - pz(1, astar, w)) * (my_inst(m, 1, w) - muintm_10) 
  eif_m_00 <- as.numeric(a == astar & z == 0) / (g(astar, w) * pz(0, astar, w)) * pz(0,aprime, w) * (my_inst(m, 0, w) - muintm_00)
  
  eif_z_11 <- as.numeric(a == astar) / g(astar, w) * muintm_11 * (z - pz(1, a,w))
  eif_z_10 <- (as.numeric(a == aprime) / g(aprime, w) - as.numeric(a == astar) / g(astar,w)) * muintm_10 * (z - pz(1, a, w))
  eif_z_00 <- -as.numeric(a == aprime) / g(aprime,w) * muintm_00 * (z - pz(1, a,w))
  
  eif_w_11 <- muintm_11 * pz(1, astar, w) 
  eif_w_10 <- muintm_10 * (pz(1, aprime, w) - pz(1, astar, w))
  eif_w_00 <- muintm_00 * pz(0, aprime, w) 
  
  eif_11 <- eif_y_11 + eif_m_11 + eif_z_11 + eif_w_11
  eif_10 <- eif_y_10 + eif_m_10 + eif_z_10 + eif_w_10
  eif_00 <- eif_y_00 + eif_m_00 + eif_z_00 + eif_w_00
  (eif_11 + eif_10 + eif_00)[, 1]
}

efficiency_bound <- function(n) {
  tmp <- simdata(n)
  eic_11 <- eic(tmp, 1, 1)
  eic_10 <- eic(tmp, 1, 0)
  eic_00 <- eic(tmp, 0, 0)
  list(nie = var(eic_11 - eic_10), 
       nde = var(eic_10 - eic_00))
}

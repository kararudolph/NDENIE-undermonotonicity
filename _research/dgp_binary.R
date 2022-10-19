suppressPackageStartupMessages(library(tidyverse))

g <- function(a, w) {
  pscore <- .5
  return(a * pscore + (1 - a) * (1 - pscore))
}

pz <- function(z, a, w) {
  prob1 <- plogis((-log(1.3)*(rowSums(w)) / 3) + 2*a - 1)
  return(z * prob1 + (1 - z) * (1 - prob1))
}

pm <- function(m, z, w) {
  prob1 <- plogis(-log(1.1)*w[, 3] + 2*z - 0.9)
  return(m * prob1 + (1 - m) * (1 - prob1))
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

my_inst <- function(m, z, w) {
  plogis((-log(1.3)*(rowSums(w)) / 3) + z + m)
}

intv11 <- function(w, aprime, astar) {
  my_inst(1, 1, w) * pz(1, astar, w) * pm(1, 1, w) + 
    my_inst(0, 1, w) * pz(1, astar, w) * pm(0, 1, w)
}

intv10 <- function(w, aprime, astar) {
  my_inst(1, 1, w) * pm(1, 0, w) * (pz(1, aprime, w) - pz(1,astar, w)) +
    my_inst(0, 1, w) * pm(0, 0, w) * (pz(1, aprime, w) - pz(1, astar, w))
}

intv00 <- function(w, aprime, astar) {
  my_inst(1, 0, w) * pm(1, 0, w) * pz(0, aprime, w) + 
    my_inst(0, 0, w) * pm(0, 0, w) * pz(0, aprime, w) 
}

muintm <- function(z1, z2, w) {
  my_inst(1, z1, w) * pm(1, z2, w) + my_inst(0, z1, w) * pm(0, z2, w)
}

simdata <- function(n = 1000) {
  w1 <- rbinom(n, 1, 0.6)
  w2 <- rbinom(n, 1, 0.3)
  w3 <- rbinom(n, 1, 0.4)
  a <- rbinom(n, 1, 0.5)
  z <- rbinom(n, 1, plogis((-log(1.3)*(w1 + w2 + w3) / 3) + 2*a - 1))
  m <- rbinom(n, 1, plogis(-log(1.1)*w3 + 2*z - 0.9))
  y <- rbinom(n, 1, plogis((-log(1.3)*(w1 + w2 + w3) / 3) + z + m))
  
  data.frame(w1 = w1, 
             w2 = w2, 
             w3 = w3, 
             a = a, 
             z = z, 
             m = m, 
             y = y)
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
  
  w_names <- str_subset(colnames(tmp), "w")
  w <- tmp[, w_names]

  aprime <- astar <- 1

  v11 <- weighted.mean(intv11(w, aprime, astar), prob_w) 
  v10 <- weighted.mean(intv10(w, aprime, astar), prob_w) 
  v00 <- weighted.mean(intv00(w, aprime, astar), prob_w) 
  
  param11 <- v11 + v10 + v00
  
  aprime <- 1
  astar <- 0

  v11 <- weighted.mean(intv11(w, aprime, astar), prob_w) 
  v10 <- weighted.mean(intv10(w, aprime, astar), prob_w) 
  v00 <- weighted.mean(intv00(w, aprime, astar), prob_w) 
  
  param10 <- v11 + v10 + v00
  
  aprime <- astar <- 0
  
  v11 <- weighted.mean(intv11(w, aprime, astar), prob_w) 
  v10 <- weighted.mean(intv10(w, aprime, astar), prob_w) 
  v00 <- weighted.mean(intv00(w, aprime, astar), prob_w) 
  
  param00 <- v11 + v10 + v00

  c("11" = param11, 
    "10" = param10, 
    "00" = param00, 
    "nie" = param11 - param10, 
    "nde" = param10 - param00)
}

eic <- function(data, aprime, astar) {
  a <- data$a
  z <- data$z
  m <- data$m
  y <- data$y
  z <- data$z
  w <- data[, paste0("w", 1:3)]
  
  eif_y_11 <- as.numeric(a == aprime & z == 1) / g(astar, w) * e(astar, 1, w) / e(aprime, 1, w) * (y - my_inst(m, z, w))
  eif_y_10 <- as.numeric(a == aprime & z == 1) / (g(astar, w) * pz(0,astar,w)) * e(astar, 1,  w) / e(aprime, 1, w) * r(0,astar,m,w)/r(1,astar, m,w) *  (pz(1, aprime, w) - pz(1, astar, w)) * (y - my_inst(m, z, w))
  eif_y_00 <- as.numeric(a == aprime & z == 0) / g(aprime, w) * e(astar, 0, w)/e(aprime, 0, w) * pz(0, aprime,w) / pz(0, astar, w) * (y - my_inst(m, z, w))
  eif_m_11 <- as.numeric(a == astar & z == 1) / g(astar,w) * (my_inst(m, 1, w) - muintm(1, 1, w))
  eif_m_10 <- as.numeric(a == astar & z == 0) / (g(astar,w) * pz(0, astar, w)) * (pz(1, aprime, w) - pz(1, astar, w)) * (my_inst(m,1, w) - muintm(1,0,w)) 
  eif_m_00 <- as.numeric(a == astar & z == 0) / (g(astar,w) * pz(0, astar, w)) * pz(0, aprime, w) * (my_inst(m,0,w) - muintm(0,0,w) )
  eif_z_11 <- as.numeric(a == astar) / g(astar, w) * muintm(1, 1, w) * (z - pz(1, a,w))
  eif_z_10 <- (as.numeric(a == aprime) / g(aprime, w) - as.numeric(a == astar) / g(astar,w))  * muintm(1,0,w) * (z - pz(1, a,w))
  eif_z_00 <- -as.numeric(a == aprime) / g(aprime,w) * muintm(0, 0, w) * (z - pz(1, a, w))
  eif_w_11 <- muintm(1, 1, w) * pz(1, astar, w) 
  eif_w_10 <- muintm(1, 0, w)  * (pz(1, aprime, w) - pz(1, astar,w))
  eif_w_00 <- muintm(0, 0, w) * pz(0,aprime, w) 
  
  eif_11 <- eif_y_11 + eif_m_11 + eif_z_11 + eif_w_11
  eif_10 <- eif_y_10 + eif_m_10 + eif_z_10 + eif_w_10
  eif_00 <- eif_y_00 + eif_m_00 + eif_z_00 + eif_w_00
  eif_11 + eif_10 + eif_00
}

efficiency_bound <- function(n) {
  tmp <- simdata(n)
  eic_11 <- eic(tmp, 1, 1)
  eic_10 <- eic(tmp, 1, 0)
  eic_00 <- eic(tmp, 0, 0)
  list(nie = var(eic_11 - eic_10), 
       nde = var(eic_10 - eic_00))
}

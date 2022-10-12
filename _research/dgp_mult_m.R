suppressPackageStartupMessages(library(tidyverse))

g <- function(a, w) {
  pscore <- .5
  return(a * pscore + (1 - a) * (1 - pscore))
}

pz <- function(z, a, w) {
  prob1 <- plogis((-log(1.3)*(rowSums(w)) / 3) + 2*a - 1)
  return(z * prob1 + (1 - z) * (1 - prob1))
}

pm1 <- function(m, z, a, w) {
  prob1 <- plogis(-log(1.1)*w[, 3] + 2*z - 0.9)
  return(m * prob1 + (1 - m) * (1 - prob1))
}

pm2 <- function(m, z, a, w) {
  prob1 <- plogis(log(2)*w[, 3] - 0.5*z - 0.2)
  return(m * prob1 + (1 - m) * (1 - prob1))
}

my_inst <- function(m1, m2, z, w) {
  plogis((-log(1.3)*(rowSums(w)) / 3) + z + m1 - 0.3*m2)
}

intv11 <- function(w, aprime, astar) {
  # m1: 1, m2: 1
  my_inst(1, 1, 1, w) * pz(1, astar, w) * pm1(1, 1, astar, w) * pm2(1, 1, astar, w) + 
    # m1: 0, m2: 1
    my_inst(0, 1, 1, w) * pz(1, astar, w) * pm1(0, 1, astar, w) * pm2(1, 1, astar, w) + 
    # m1: 1, m2: 0
    my_inst(1, 0, 1, w) * pz(1, astar, w) * pm1(1, 1, astar, w) * pm2(0, 1, astar, w) + 
    # m1: 0, m2: 0
    my_inst(0, 0, 1, w) * pz(1, astar, w) * pm1(0, 1, astar, w) * pm2(0, 1, astar, w)
}

intv10 <- function(w, aprime, astar) {
  # m1: 1, m2: 1
  my_inst(1, 1, 1, w) * pm1(1, 0, astar, w) * pm2(1, 0, astar, w) * (pz(1, aprime, w) - pz(1,astar, w)) +
    # m1: 0, m2: 1
    my_inst(0, 1, 1, w) * pm1(0, 0, astar, w) * pm2(1, 0, astar, w) * (pz(1, aprime, w) - pz(1, astar, w)) + 
    # m1: 1, m2: 0
    my_inst(1, 0, 1, w) * pm1(1, 0, astar, w) * pm2(0, 0, astar, w) * (pz(1, aprime, w) - pz(1, astar, w)) + 
    # m1: 0, m2: 0
    my_inst(0, 0, 1, w) * pm1(0, 0, astar, w) * pm2(0, 0, astar, w) * (pz(1, aprime, w) - pz(1, astar, w))
}

intv00 <- function(w, aprime, astar) {
  # m1: 1, m2: 1
  my_inst(1, 1, 0, w) * pm1(1, 0, astar, w) * pm2(1, 0, astar, w) * pz(0, aprime, w) + 
    # m1: 0, m2: 1
    my_inst(0, 1, 0, w) * pm1(0, 0, astar, w) * pm2(1, 0, astar, w) * pz(0, aprime, w) + 
    # m1: 1, m2: 0
    my_inst(1, 0, 0, w) * pm1(1, 0, astar, w) * pm2(0, 0, astar, w) * pz(0, aprime, w) + 
    # m1: 0, m2: 0
    my_inst(0, 0, 0, w) * pm1(0, 0, astar, w) * pm2(0, 0, astar, w) * pz(0, aprime, w)
}

simdata <- function(n = 1000) {
  w1 <- rbinom(n, 1, 0.6)
  w2 <- rbinom(n, 1, 0.3)
  w3 <- rbinom(n, 1, 0.4)
  a <- rbinom(n, 1, 0.5)
  z <- rbinom(n, 1, plogis((-log(1.3)*(w1 + w2 + w3) / 3) + 2*a - 1))
  m1 <- rbinom(n, 1, plogis(-log(1.1)*w3 + 2*z - 0.9))
  m2 <- rbinom(n, 1, plogis(log(2)*w3 - 0.5*z - 0.2))
  y <- rbinom(n, 1, plogis((-log(1.3)*(w1 + w2 + w3) / 3) + z + m1 - 0.3*m2))
  
  data.frame(w1 = w1, 
             w2 = w2, 
             w3 = w3, 
             a = a, 
             z = z, 
             m1 = m1,
             m2 = m2,
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

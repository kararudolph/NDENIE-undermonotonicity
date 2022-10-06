suppressPackageStartupMessages(library(tidyverse))

g <- function(a, w) {
  pscore <- .5
  return(a * pscore + (1 - a) * (1 - pscore))
}

pz <- function(z, a, w) {
  prob1 <- plogis((-log(1.3)*(rowSums(w)) / 3) + 2*a - 1)
  return(z * prob1 + (1 - z) * (1 - prob1))
}

pm <- function(m, z, a, w) {
  prob1 <- plogis(-log(1.1)*w[, 3] + 2*z - 0.9)
  return(m * prob1 + (1 - m) * (1 - prob1))
}

my_inst <- function(m, z, w) {
  plogis((-log(1.3)*(rowSums(w)) / 3) + z + m)
}

pzmw <- function(m, z, w) {
  pm(m, z, 1, w) * g(1, w) + pm(m, z, 0, w) * g(0, w)
}

e <- function(a, z, m, w) {
  pm(m, z, a, w) * g(a, w) / pzmw(m, z,w)
}

pmaw <- function(m, a, w) {
  pm(m, 1, a, w) * pz(1, a, w) + pm(m, 0, a, w) * pz(0, a, w)
}

r <- function(z, a, m, w) {
  pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)
}

muintm <- function(z1, z2, a, w){
  my_inst(1, z1, w) * pm(1, z2, a, w) + my_inst(0, z1, w) * pm(0, z2, a, w)
}

intv11 <- function(w, aprime, astar) {
  my_inst(1, 1, w) * pz(1, astar, w) * pm(1, 1, astar, w) + 
    my_inst(0, 1, w) * pz(1, astar, w) * pm(0, 1, astar, w)
}

intv10 <- function(w, aprime, astar) {
  my_inst(1, 1, w) * pm(1, 0, astar, w) * (pz(1, aprime, w) - pz(1,astar, w)) +
    my_inst(0, 1, w) * pm(0, 0, astar, w)* (pz(1, aprime, w) - pz(1, astar, w))
}

intv00 <- function(w, aprime, astar) {
  my_inst(1, 0, w) * pm(1, 0, astar, w) * pz(0, aprime, w) + 
    my_inst(0, 0, w) * pm(0, 0, astar, w) * pz(0, aprime, w) 
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

truth_monotonicity_mediation <- function(data) {
  w_names <- str_subset(colnames(data), "w")
  m_names <- str_subset(colnames(data), "m")
  w <- data[, w_names]
  a <- data$a
  z <- data$z
  m <- data$m
  y <- data$y

  aprime <- astar <- 1

  v11 <- mean(intv11(w, aprime, astar)) 
  v10 <- mean(intv10(w, aprime, astar)) 
  v00 <- mean(intv00(w, aprime, astar)) 
  
  param11 <- v11 + v10 + v00
  
  aprime <- 1
  astar <- 0

  v11 <- mean(intv11(w, aprime, astar)) 
  v10 <- mean(intv10(w, aprime, astar)) 
  v00 <- mean(intv00(w, aprime, astar)) 
  
  param10 <- v11 + v10 + v00
  
  aprime <- astar <- 0
  
  v11 <- mean(intv11(w, aprime, astar)) 
  v10 <- mean(intv10(w, aprime, astar)) 
  v00 <- mean(intv00(w, aprime, astar)) 
  
  param00 <- v11 + v10 + v00

  c("11" = param11, 
    "10" = param10, 
    "00" = param00, 
    "nie" = param11 - param10, 
    "nde" = param10 - param00)
}

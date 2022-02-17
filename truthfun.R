bound_precision <- function(vals, tol = 1e-6) {
  vals[vals < tol] <- tol
  vals[vals > 1 - tol] <- 1 - tol
  return(vals)
}

bound_propensity <- function(vals, bounds = c(0.001, 0.999)) {
  assertthat::assert_that(!(max(vals) > 1 || min(vals) < 0))
  vals[vals < bounds[1]] <- bounds[1]
  vals[vals > bounds[2]] <- bounds[2]
  return(vals)
}

scale_to_unit <- function(vals) {
  vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals))
  return(vals_scaled)
}

scale_from_unit <- function(scaled_vals, max_orig, min_orig) {
  vals_orig <- (scaled_vals * (max_orig - min_orig)) + min_orig
  return(vals_orig)
}

g <- function(a, w) {
    pscore <- .5
    return(a * pscore + (1 - a) * (1 - pscore))
}

pz <- function(z, a, w) {
    prob1 <- plogis(rowMeans(-log(1.1)* w) + 3*a)
    return(z * prob1 + (1 - z) * (1 - prob1))
}

pm <- function(m, z, a, w) {
    prob1 <- plogis(rowSums(log(3) * w[, -3] + a - z))
    return(m * prob1 + (1 - m) * (1 - prob1))
}

my_inst <- function(m, z, w) {
   plogis(1 / as.numeric(unlist(-1*rowSums(w) + z + m)) -1)
}

pzmw <- function(m, z,w) {
    pm(m,z, 1, w) * g(1, w) +
        pm(m,z, 0, w) * g(0, w)
}

e <- function(a, z, m, w) {
    pm(m, z, a, w) * g(a, w) / pzmw(m, z,w)
}

pmaw <- function(m, a, w) {
    pm(m, 1, a, w) * pz(1, a, w) +
        pm(m, 0, a, w) * pz(0, a, w)
}

r <- function( z, a, m, w) {
    pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)
}

muintm <- function(z1,z2,a,w){
    my_inst(1,z1,w)*pm(1,z2,a,w) + my_inst(0,z1,w)*pm(0,z2,a,w)
}

simdata <- function(n_obs = 1000) {

    ## helper functions for nuisance parameter estimation
    ## baseline covariate -- simple, binary
    w_1 <- rbinom(n_obs, 1, prob = 0.6)
    w_2 <- rbinom(n_obs, 1, prob = 0.3)
    w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
    w <- cbind(w_1, w_2, w_3)
    w_names <- paste("W", seq_len(ncol(w)), sep = "")

    ## exposure/treatment
    a <- as.numeric(rbinom(n_obs, 1, prob = g(1, w)))

    ## mediator-outcome confounder affected by treatment
    z <- rbinom(n_obs, 1, pz(1, a, w))

    ## mediator (possibly multivariate)
    m <- rbinom(n_obs, 1, pm(1, z, a, w))
    m_names <- "M"

    ## outcome
    y <- rbinom(n_obs, 1, my_inst(m, z, w))

    colnames(w) <- w_names
    ## construct data for output
    dat <- as.data.frame(cbind(W = w, A = a, Z = z, M = m, Y = y))
    return(dat)
}


#n<-1000000
#set.seed(63845)

truth_monotonicity_mediation <- function(n){

    library(data.table)
    library(tidyverse)
    data <- simdata(n_obs = n)
    w_names <- str_subset(colnames(data), "W")
    m_names <- str_subset(colnames(data), "M")
    w <- as_tibble(data[, w_names])
    a <- data$A
    z <- data$Z
    m <- data$M
    y <- data$Y

    ## compute quasi-substitution estimator of parameter
    aprime <- astar <- 1
    intv11 <- function(w, aprime, astar) {
    
    #my_inst(m, rep(1,n), w) * pz(1,astar, w) * pm(m, 1, astar, w)
    #KER check that this is ohw to do it
    my_inst(rep(1,n), rep(1,n), w) * pz(1,astar, w) * pm(1, 1, astar, w) + my_inst(rep(0,n), rep(1,n), w) * pz(1,astar, w) * pm(0, 1, astar, w)
}

    v11 <- mean(intv11(w, rep(aprime,n), rep(astar,n))) 
    # 0.5358498

    intv10 <- function(w, aprime, astar) {
    my_inst(rep(1,n), rep(1,n), w) * pm(1, 0, astar, w)* (pz(1,aprime, w) - pz(1,astar,w))  +  my_inst(rep(0,n), rep(1,n), w) * pm(0, 0, astar, w)* (pz(1,aprime, w) - pz(1,astar,w))
}

    v10 <- mean(intv10(w, rep(aprime,n), rep(astar,n))) 

    
    intv00 <- function(w, aprime, astar) {
    my_inst(rep(1,n), rep(0,n), w) * pm(1, 0, astar, w)* pz(0,aprime, w) + my_inst(rep(0,n), rep(0,n), w) * pm(0, 0, astar, w)* pz(0,aprime, w) 
}

    v00 <- mean(intv00(w,  rep(aprime,n), rep(astar,n))) 

    parama1a1<- v11 + v10 + v00
#0.55948

 aprime <- 1
 astar <- 0

    v11 <- mean(intv11(w, rep(aprime,n), rep(astar,n))) 


    v10 <- mean(intv10(w, rep(aprime,n), rep(astar,n))) 


    v00 <- mean(intv00(w,  rep(aprime,n), rep(astar,n))) 
    
    param10 <- v11 + v10 + v00
#0.5390065

    aprime <- 0
    astar <- 0

    v11 <- mean(intv11(w, rep(aprime,n), rep(astar,n))) 


    v10 <- mean(intv10(w, rep(aprime,n), rep(astar,n))) 


    v00 <- mean(intv00(w,  rep(aprime,n), rep(astar,n))) 
  
    param00<- v11 + v10 + v00
    #0.4927015

    nie<- parama1a1 - param10
    nde<- param10 - param00
    #nie:  0.02047349
    #nde: 0.04630506

aprime <- astar <- 1

  eif_y_11 <- as.numeric(a == aprime & z == 1) / g(astar, w) * e(astar, 1, m, w) / e(aprime, 1, m, w) * (y - my_inst(m, z, w))
  eif_y_10 <- as.numeric(a == aprime & z==1) / (g(astar, w) * pz(0,astar,w)) * e(astar, 1, m, w)/e(aprime, 1, m, w) * r(0,astar,m,w)/r(1,astar, m,w) *  (pz(1, aprime, w) - pz(1, astar, w)) * (y - my_inst(m, z, w))
eif_y_00 <- as.numeric(a == aprime & z==0) / g(aprime, w) * e(astar, 0, m, w)/e(aprime, 0, m, w) * pz(0, aprime,w) / pz(0, astar, w) * (y - my_inst(m, z, w))

  eif_m_11 <- as.numeric(a == astar & z==1) / g(astar,w) * (my_inst(m,1,w) - muintm(1, 1, astar, w))

  eif_m_10 <- as.numeric(a == astar & z==0) / (g(astar,w) * pz(0,astar,w)) * (pz(1, aprime, w) - pz(1, astar, w)) * (my_inst(m,1, w) - muintm(1,0,astar,w)) 

  eif_m_00 <- as.numeric(a == astar & z==0) / (g(astar,w) * pz(0,astar,w)) * pz(0,aprime,w) * (my_inst(m,0,w) - muintm(0,0,astar,w) )

  eif_z_11 <- as.numeric(a == astar) / g(astar, w) * muintm(1,1,astar,w) * (z - pz(1, a,w))

  eif_z_10 <- (as.numeric(a == aprime) / g(aprime, w) - as.numeric(a == astar) / g(astar,w))  * muintm(1,0,astar,w) * (z - pz(1, a,w))

  eif_z_00 <- -1* as.numeric(a == aprime) / g(aprime,w) * muintm(0,0,astar,w) * (z - pz(1, a,w))

  eif_w_11 <-  muintm(1,1,astar,w) * pz(1,astar,w) 

  eif_w_10 <- muintm(1,0,astar,w)  * (pz(1, aprime,w) - pz(1, astar,w))

  eif_w_00 <- muintm(0,0,astar,w) * pz(0,aprime,w) 

# un-centered efficient influence function
eif_11 <- eif_y_11 + eif_m_11 + eif_z_11 + eif_w_11
eif_10 <- eif_y_10 + eif_m_10 + eif_z_10 + eif_w_10
eif_00 <- eif_y_00 + eif_m_00 + eif_z_00 + eif_w_00
eifparam11 <- eif_11 + eif_10 + eif_00

aprime <- 1
astar <- 0 

   eif_y_11 <- as.numeric(a == aprime & z == 1) / g(astar, w) * e(astar, 1, m, w) / e(aprime, 1, m, w) * (y - my_inst(m, z, w))
eif_y_10 <- as.numeric(a == aprime & z==1) / (g(astar, w) * pz(0,astar,w)) * e(astar, 1, m, w)/e(aprime, 1, m, w) * r(0,astar,m,w)/r(1,astar, m,w) *  (pz(1, aprime, w) - pz(1, astar, w)) * (y - my_inst(m, z, w))
eif_y_00 <- as.numeric(a == aprime & z==0) / g(aprime, w) * e(astar, 0, m, w)/e(aprime, 0, m, w) * pz(0, aprime,w) / pz(0, astar, w) * (y - my_inst(m, z, w))

  eif_m_11 <- as.numeric(a == astar & z==1) / g(astar,w) * (my_inst(m,1,w) - muintm(1, 1, astar, w))

  eif_m_10 <- as.numeric(a == astar & z==0) / (g(astar,w) * pz(0,astar,w)) * (pz(1, aprime, w) - pz(1, astar, w)) * (my_inst(m,1, w) - muintm(1,0,astar,w)) 

  eif_m_00 <- as.numeric(a == astar & z==0) / (g(astar,w) * pz(0,astar,w)) * pz(0,aprime,w) * (my_inst(m,0,w) - muintm(0,0,astar,w) )

  eif_z_11 <- as.numeric(a == astar) / g(astar, w) * muintm(1,1,astar,w) * (z - pz(1, a,w))

  eif_z_10 <- (as.numeric(a == aprime) / g(aprime, w) - as.numeric(a == astar) / g(astar,w))  * muintm(1,0,astar,w) * (z - pz(1, a,w))

  eif_z_00 <- -1* as.numeric(a == aprime) / g(aprime,w) * muintm(0,0,astar,w) * (z - pz(1, a,w))

  eif_w_11 <-  muintm(1,1,astar,w) * pz(1,astar,w) 

  eif_w_10 <- muintm(1,0,astar,w)  * (pz(1, aprime,w) - pz(1, astar,w))

  eif_w_00 <- muintm(0,0,astar,w) * pz(0,aprime,w) 
# un-centered efficient influence function
eif_11 <- eif_y_11 + eif_m_11 + eif_z_11 + eif_w_11
eif_10 <- eif_y_10 + eif_m_10 + eif_z_10 + eif_w_10
eif_00 <- eif_y_00 + eif_m_00 + eif_z_00 + eif_w_00
eifparam10 <- eif_11 + eif_10 + eif_00

aprime <- 0
astar <- 0 

    eif_y_11 <- as.numeric(a == aprime & z == 1) / g(astar, w) * e(astar, 1, m, w) / e(aprime, 1, m, w) * (y - my_inst(m, z, w))
eif_y_10 <- as.numeric(a == aprime & z==1) / (g(astar, w) * pz(0,astar,w)) * e(astar, 1, m, w)/e(aprime, 1, m, w) * r(0,astar,m,w)/r(1,astar, m,w) *  (pz(1, aprime, w) - pz(1, astar, w)) * (y - my_inst(m, z, w))
eif_y_00 <- as.numeric(a == aprime & z==0) / g(aprime, w) * e(astar, 0, m, w)/e(aprime, 0, m, w) * pz(0, aprime,w) / pz(0, astar, w) * (y - my_inst(m, z, w))

  eif_m_11 <- as.numeric(a == astar & z==1) / g(astar,w) * (my_inst(m,1,w) - muintm(1, 1, astar, w))

  eif_m_10 <- as.numeric(a == astar & z==0) / (g(astar,w) * pz(0,astar,w)) * (pz(1, aprime, w) - pz(1, astar, w)) * (my_inst(m,1, w) - muintm(1,0,astar,w)) 

  eif_m_00 <- as.numeric(a == astar & z==0) / (g(astar,w) * pz(0,astar,w)) * pz(0,aprime,w) * (my_inst(m,0,w) - muintm(0,0,astar,w) )

  eif_z_11 <- as.numeric(a == astar) / g(astar, w) * muintm(1,1,astar,w) * (z - pz(1, a,w))

  eif_z_10 <- (as.numeric(a == aprime) / g(aprime, w) - as.numeric(a == astar) / g(astar,w))  * muintm(1,0,astar,w) * (z - pz(1, a,w))

  eif_z_00 <- -1* as.numeric(a == aprime) / g(aprime,w) * muintm(0,0,astar,w) * (z - pz(1, a,w))

  eif_w_11 <-  muintm(1,1,astar,w) * pz(1,astar,w) 

  eif_w_10 <- muintm(1,0,astar,w)  * (pz(1, aprime,w) - pz(1, astar,w))

  eif_w_00 <- muintm(0,0,astar,w) * pz(0,aprime,w) 
# un-centered efficient influence function
eif_11 <- eif_y_11 + eif_m_11 + eif_z_11 + eif_w_11
eif_10 <- eif_y_10 + eif_m_10 + eif_z_10 + eif_w_10
eif_00 <- eif_y_00 + eif_m_00 + eif_z_00 + eif_w_00
eifparam00 <- eif_11 + eif_10 + eif_00

    var_indirect <- var(eifparam11 - eifparam10)
    var_direct <- var(eifparam10 - eifparam00)
    return(list(estimates = c(nie, nde), eff_bound = c(var_indirect, var_direct)))

}
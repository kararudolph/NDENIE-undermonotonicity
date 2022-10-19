library(purrr)
library(sl3)
suppressPackageStartupMessages(library(tidyverse))
source("_research/Lrnr_glmnet3.R")

devtools::load_all("monomediate")

args <- commandArgs(trailingOnly = TRUE)
dgp <- args[[1]]

if (dgp == "binary") {
  source("_research/dgp_binary.R")
  y_family <- "binomial"
  npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", "m", "y")
}

if (dgp == "cont") {
  source("_research/dgp_cont.R")
  y_family <- "gaussian"
  npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", "m", "y")
}

if (dgp == "mult") {
  source("_research/dgp_mult.R")
  y_family <- "binomial"
  npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", c("m1", "m2"), "y")
}

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

spec <- as.numeric(args[[2]])

specs <- list(
  # All correct
  list(g = Lrnr_mean$new(), 
       q = Lrnr_glmnet3$new(), 
       e = Lrnr_glmnet3$new(), 
       r = Lrnr_glmnet3$new(), 
       mu = Lrnr_glmnet3$new(), 
       rho = Lrnr_glmnet3$new()), 
  # g, q, e, r correct
  list(g = Lrnr_mean$new(), 
       q = Lrnr_glmnet3$new(), 
       e = Lrnr_glmnet3$new(), 
       r = Lrnr_glmnet3$new(), 
       mu = Lrnr_mean$new(), 
       rho = Lrnr_mean$new()), 
  # mu, rho, g correct
  list(g = Lrnr_mean$new(), 
       q = Lrnr_mean$new(), 
       e = Lrnr_mean$new(), 
       r = Lrnr_mean$new(), 
       mu = Lrnr_glmnet3$new(), 
       rho = Lrnr_glmnet3$new()), 
  # mu, rho, q correct
  list(g = Lrnr_mean$new(), 
       q = Lrnr_glmnet3$new(), 
       e = Lrnr_mean$new(), 
       r = Lrnr_mean$new(), 
       mu = Lrnr_glmnet3$new(), 
       rho = Lrnr_glmnet3$new()), 
  # g, q, mu correct
  list(g = Lrnr_mean$new(), 
       q = Lrnr_glmnet3$new(), 
       e = Lrnr_mean$new(), 
       r = Lrnr_mean$new(), 
       mu = Lrnr_glmnet3$new(), 
       rho = Lrnr_mean$new())
)

res <- map_dfr(c(1000, 10000), function(n) {
  tmp <- simdata(n)
  
  folds <- ifelse(n == 1000, 10, 2)
  
  dr <- monomediate(tmp, npsem, specs[[spec]], folds, y_family)
  
  data.frame(estimator = rep("onestep", 2), 
             n = rep(n, 2), 
             spec = rep(spec, 2),
             estimand = c("nie", "nde"), 
             psi = c(dr$nie, dr$nde), 
             var = c(dr$var_nie, dr$var_nde),
             conf_low = c(dr$ci_nie[1], dr$ci_nde[1]), 
             conf_high = c(dr$ci_nie[2], dr$ci_nde[2]))
})

write_csv(res, paste0("_research/data/sim_", args[[1]], "_onestep_", args[[2]], "_", id, ".csv"))

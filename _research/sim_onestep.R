library(boot)
library(purrr)
library(sl3)
source("_research/dgp_all_binary.R")
source("_research/Lrnr_glmnet3.R")

devtools::load_all("monomediate")

npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", "m", "y")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1
args <- commandArgs(trailingOnly = TRUE)
spec <- as.numeric(args[[1]])

specs <- list(
  list(g = Lrnr_mean$new(), 
       q = Lrnr_glmnet3$new(), 
       e = Lrnr_glmnet3$new(), 
       r = Lrnr_glmnet3$new(), 
       mu = Lrnr_glmnet3$new(), 
       rho = Lrnr_glmnet3$new()), 
  list(g = Lrnr_mean$new(), 
       q = Lrnr_glmnet3$new(), 
       e = Lrnr_glmnet3$new(), 
       r = Lrnr_glmnet3$new(), 
       mu = Lrnr_mean$new(), 
       rho = Lrnr_mean$new()), 
  list(g = Lrnr_mean$new(), 
       q = Lrnr_mean$new(), 
       e = Lrnr_mean$new(), 
       r = Lrnr_mean$new(), 
       mu = Lrnr_glmnet3$new(), 
       rho = Lrnr_glmnet3$new()), 
  list(g = Lrnr_mean$new(), 
       q = Lrnr_glmnet3$new(), 
       e = Lrnr_mean$new(), 
       r = Lrnr_mean$new(), 
       mu = Lrnr_glmnet3$new(), 
       rho = Lrnr_glmnet3$new())
)

res <- map_dfr(c(1000, 10000), function(n) {
  tmp <- simdata(n)
  
  folds <- ifelse(n == 1000, 10, 2)
  
  dr <- monomediate(tmp, npsem, specs[[spec]], folds)
  
  data.frame(estimator = rep("onestep", 2), 
             n = rep(n, 2), 
             spec = rep(spec, 2),
             estimand = c("nie", "nde"), 
             psi = c(dr$nie, dr$nde), 
             conf_low = c(dr$ci_nie[1], dr$ci_nde[1]), 
             conf_high = c(dr$ci_nie[2], dr$ci_nde[2]))
})

write_csv(res, paste0("_research/data/sim_binary_onestep_", args[[1]], "_", id, ".csv"))

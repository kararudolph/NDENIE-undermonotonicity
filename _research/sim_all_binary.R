library(boot)
library(purrr)
library(sl3)
source("_research/dgp_all_binary.R")
source("_research/gcomputation.R")
source("_research/Lrnr_glmnet3.R")

devtools::load_all("pkg")

gmboot <- function(data, i) {
  gcomp_mediation_monotinicity(data[i, ])
}

npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", "m", "y")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

map_dfr(c(1000, 10000), function(n) {
  tmp <- simdata(n)
   
  booted <- boot(tmp, statistic = gmboot, R = 500)
  ci_nie <- boot.ci(booted, type = "norm", index = 1)
  ci_nde <- boot.ci(booted, type = "norm", index = 2)
  
  folds <- ifelse(n == 1000, 10, 2)
  
  dr <- monomediate(tmp, npsem, Lrnr_glmnet3$new(), folds)
  
  data.frame(n = rep(n, 2), 
             estimand = c("nie", "nde", "dr_nie", "dr_nde"), 
             psi = c(as.vector(booted$t0), dr$nie, dr$nde), 
             conf_low = c(ci_nie$normal[2], ci_nde$normal[2], dr$ci_nie[1], dr$ci_nde[1]), 
             conf_high = c(ci_nie$normal[3], ci_nde$normal[3], dr$ci_nie[2], dr$ci_nde[2]))
})

library(boot)
library(purrr)
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
dgp <- args[[1]]

if (dgp == "binary") {
  source("_research/dgp_all_binary.R")
  y_family <- "binomial"
}

if (dgp == "cont") {
  source("_research/dgp_cont.R")
  y_family <- "gaussian"
}

devtools::load_all("monomediate")

specs <- list(
  correct = c(y ~ w1 + w2 + w3 + z + m, 
              z ~ w1 + w2 + w3 + a, 
              ~ w3 + z),
  "y incorrect" = c(y ~ 1, 
                    z ~ w1 + w2 + w3 + a, 
                    ~ w3 + z), 
  "z incorrect" = c(y ~ w1 + w2 + w3 + z + m, 
                    z ~ 1, 
                    ~ w3 + z), 
  "m incorrect" = c(y ~ w1 + w2 + w3 + z + m, 
                    z ~ w1 + w2 + w3 + a, 
                    ~ 1)
)

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

spec <- specs[[as.numeric(args[[2]])]]

gmboot <- function(data, i) {
  monomediate_gcomp(data[i, ], 
                    npsem, 
                    spec[[1]], 
                    spec[[2]], 
                    spec[[3]], 
                    y_family)
}

npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", "m", "y")

res <- map_dfr(c(1000, 10000), function(n) {
  tmp <- simdata(n)
  
  booted <- boot(tmp, statistic = gmboot, R = 500)
  ci_nie <- boot.ci(booted, type = "norm", index = 1)
  ci_nde <- boot.ci(booted, type = "norm", index = 2)
  
  data.frame(estimator = rep("gcomp", 2), 
             n = rep(n, 2), 
             spec = rep(names(specs)[as.numeric(args[[2]])], 2),
             estimand = c("nie", "nde"), 
             psi = as.vector(booted$t0), 
             conf_low = c(ci_nie$normal[2], ci_nde$normal[2]), 
             conf_high = c(ci_nie$normal[3], ci_nde$normal[3]))
})

write_csv(res, paste0("_research/data/sim_", args[[1]], "_gcomp_", args[[2]], "_", id, ".csv"))

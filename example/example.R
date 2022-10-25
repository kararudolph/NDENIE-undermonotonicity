# Install the R package from source
# the package path may need to change depending on the analysts 
# machine and working directory
pkg_path <- "monomediate_0.1.0.tar.gz"
install.packages(pkg_path, type = "source", repos = NULL)

library(monomediate)
library(sl3)

# DGP with multiple M
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

tmp <- simdata()

# Create an Npsem object that maps variables in the simulated data
# to the assumed DAG
npsem <- Npsem$new(W = c("w1", "w2", "w3"), 
                   A = "a", 
                   Z = "z", 
                   M = c("m1", "m2"), 
                   Y = "y")

# Define a list of sl3 learners
lrnrs <- list(g = Lrnr_mean$new(), 
              q = Lrnr_glm$new(), 
              e = Lrnr_glm$new(), 
              r = Lrnr_glm$new(), 
              mu = Lrnr_glm$new(), 
              rho = Lrnr_xgboost$new())

# Estimate effects
monomediate(tmp, npsem, lrnrs, 10, "binomial")

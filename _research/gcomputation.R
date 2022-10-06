# Estimates the NDE and NIE using g-computation from Tchetgen Tchetgen and VanderWeele (2014)
# This only works for the specific DGP evaluated in the simulations
# To work with other DGPs, change `y_fit`, `m_fit`, and `z_fit`
# A, Z, M must all be binary
# Std. errors can be calculated usign bootstrap
gcomp_mediation_monotinicity <- function(data) {
  suppressPackageStartupMessages(require("dplyr"))
  
  tmp <- data
  y_fit <- glm(y ~ w1 + w2 + w3 + z + m, family = "binomial", data = tmp) # misspecify for 2
  m_fit <- glm(m ~ w3 + z, family = "binomial", data = tmp) # misspecify for 4
  z_fit <- glm(z ~ w1 + w2 + w3 + a, family = "binomial", data = tmp) # misspecify for 3
  
  pr_z_1 <- predict(z_fit, mutate(tmp, a = 1), type = "response")
  pr_z_0 <- predict(z_fit, mutate(tmp, a = 0), type = "response")
  
  pr_m1_a0_z1 <- predict(m_fit, mutate(tmp, a = 0, z = 1), type = "response")
  pr_m1_a0_z0 <- predict(m_fit, mutate(tmp, a = 0, z = 0), type = "response")
  pr_m1_a1_z1 <- predict(m_fit, mutate(tmp, a = 1, z = 1), type = "response")
  pr_m1_a1_z0 <- predict(m_fit, mutate(tmp, a = 1, z = 0), type = "response")
  
  pr_m0_a0_z1 <- 1 - pr_m1_a0_z1
  pr_m0_a0_z0 <- 1 - pr_m1_a0_z0
  pr_m0_a1_z1 <- 1 - pr_m1_a1_z1
  pr_m0_a1_z0 <- 1 - pr_m1_a1_z0
  
  pr_y_e1_z1_m1 <- predict(y_fit, mutate(tmp, a = 1, z = 1, m = 1), type = "response")
  pr_y_e1_z1_m0 <- predict(y_fit, mutate(tmp, a = 1, z = 1, m = 0), type = "response")
  pr_y_e1_z0_m1 <- predict(y_fit, mutate(tmp, a = 1, z = 0, m = 1), type = "response")
  pr_y_e1_z0_m0 <- predict(y_fit, mutate(tmp, a = 1, z = 0, m = 0), type = "response")
  pr_y_e0_z1_m1 <- predict(y_fit, mutate(tmp, a = 0, z = 1, m = 1), type = "response")
  pr_y_e0_z1_m0 <- predict(y_fit, mutate(tmp, a = 0, z = 1, m = 0), type = "response")
  pr_y_e0_z0_m1 <- predict(y_fit, mutate(tmp, a = 0, z = 0, m = 1), type = "response")
  pr_y_e0_z0_m0 <- predict(y_fit, mutate(tmp, a = 0, z = 0, m = 0), type = "response")
  
  # start with a = 1, a' = 0
  theta_11 <- (pr_y_e1_z1_m0 * pr_m0_a0_z1 * pr_z_0) + 
    (pr_y_e1_z1_m1 * pr_m1_a0_z1 * pr_z_0)
  
  theta_10 <- (pr_y_e1_z1_m0 * pr_m0_a0_z0 * (pr_z_1 - pr_z_0)) + 
    (pr_y_e1_z1_m1 * pr_m1_a0_z0 * (pr_z_1 - pr_z_0))
  
  theta_00 <- (pr_y_e1_z0_m0 * pr_m0_a0_z0 * (1 - pr_z_1)) + 
    (pr_y_e1_z0_m1 * pr_m1_a0_z0 * (1 - pr_z_1))
  
  ey_10 <- mean(theta_11 + theta_10 + theta_00)
  
  # a = 0, a' = 0
  theta_11 <- (pr_y_e0_z1_m0 * pr_m0_a0_z1 * pr_z_0) + 
    (pr_y_e0_z1_m1 * pr_m1_a0_z1 * pr_z_0)
  
  theta_10 <- (pr_y_e0_z1_m0 * pr_m0_a0_z0 * (pr_z_0 - pr_z_0)) + 
    (pr_y_e0_z1_m1 * pr_m1_a0_z0 * (pr_z_0 - pr_z_0))
  
  theta_00 <- (pr_y_e0_z0_m0 * pr_m0_a0_z0 * (1 - pr_z_0)) + 
    (pr_y_e0_z0_m1 * pr_m1_a0_z0 * (1 - pr_z_0))
  
  ey_00 <- mean(theta_11 + theta_10 + theta_00)
  
  # a = 1, a' = 1
  theta_11 <- (pr_y_e1_z1_m0 * pr_m0_a1_z1 * pr_z_1) + 
    (pr_y_e1_z1_m1 * pr_m1_a1_z1 * pr_z_1)
  
  theta_10 <- (pr_y_e1_z1_m0 * pr_m0_a1_z0 * (pr_z_1 - pr_z_1)) + 
    (pr_y_e1_z1_m1 * pr_m1_a1_z0 * (pr_z_1 - pr_z_1))
  
  theta_00 <- (pr_y_e1_z0_m0 * pr_m0_a1_z0 * (1 - pr_z_1)) + 
    (pr_y_e1_z0_m1 * pr_m1_a1_z0 * (1 - pr_z_1))
  
  ey_11 <- mean(theta_11 + theta_10 + theta_00)
  
  c("nie" = ey_11 - ey_10, 
    "nde" = ey_10 - ey_00)
}

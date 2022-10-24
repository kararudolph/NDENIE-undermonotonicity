#' One-step mediation estimator under monotonicity assumption
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} containing all necessary variables
#'  for the estimation problem.
#' @param npsem \[\code{R6(Npsem)}\]\cr
#'  An \code{Npsem} object 
#' @param learners \[\code{list(6)}\]\cr
#'  A named list ("g", "q", "r", "e", "mu", "rho") of \code{sl3} learners.
#' @param folds \[\code{integer(1)}\]\cr
#'  The number of folds to be used for cross-fitting.
#' @param y_family \[\code{character(1)}\]\cr
#'  Outcome variable type (i.e., "gaussian", "binomial").
#'
#' @return Natural (in)direct effect estimates and estimated variance
#' @export
#'
#' @examples
#' library(sl3)
#' npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", c("m1", "m2"), "y")
#' lrnrs <- list(g = Lrnr_mean$new(), 
#'               q = Lrnr_glm$new(), 
#'               e = Lrnr_glm$new(), 
#'               r = Lrnr_glm$new(), 
#'               mu = Lrnr_glm$new(), 
#'               rho = Lrnr_glm$new())
#' 
#' monomediate(multiple_m, npsem, lrnrs, 10, "binomial")
monomediate <- function(data, npsem, learners, folds, y_family) {
  tmp <- data
  folds <- origami::make_folds(tmp, V = folds)

  g <- cv_g(tmp, folds, npsem, learners$g)
  q <- cv_q(tmp, folds, npsem, learners$q)
  r <- cv_r(tmp, folds, npsem, learners$r)
  e <- cv_e(tmp, folds, npsem, learners$e)
  mu <- cv_mu(tmp, folds, npsem, learners$mu, y_family)

  tmp <- cbind(tmp, as.data.frame(mu))
  
  rhos <- list()
  for (as in list(c(1, 1), c(1, 0), c(0, 0))) {
    for (zs in list(c(1, 1), c(1, 0), c(0, 0))) {
      a <- as[1]
      ap <- as[2]
      z <- zs[1]
      zp <- zs[2]
      
      rhos[[paste(a, ap, z, zp, sep = "_")]] <- 
        cv_rho(tmp, folds, npsem, learners$rho, a, ap, z, zp)
    }
  }
  
  Z <- tmp[[npsem$Z]]
  A <- tmp[[npsem$A]]
  Y <- tmp[[npsem$Y]]
  
  eics <- list()
  for (as in list(c(1, 1), c(1, 0), c(0, 0))) {
    a <- as[1]
    ap <- as[2]
    
    `P(a'|W)` <- g[[paste0("a", ap)]]
    `P(a|W)` <- g[[paste0("a", a)]]
    `P(a|M,1,W)` <- e[[paste0("a", a, "_z1")]]
    `P(a'|M,1,W)` <- e[[paste0("a", ap, "_z1")]]
    `P(Z=0|a',W)` <- q[[paste0("z0_a", ap)]]
    `P(Z=0|a,W)` <- q[[paste0("z0_a", a)]]
    `P(Z=1|a',W)` <- q[[paste0("z1_a", ap)]]
    `P(Z=1|a,W)` <- q[[paste0("z1_a", a)]]
    `P(Z=0|M,a',W)` <- r[[paste0("z0_a", ap)]]
    `P(Z=1|M,a',W)` <- r[[paste0("z1_a", ap)]]
    `P(a|M,0,W)` <- e[[paste0("a", a, "_z0")]]
    `P(a'|M,0,W)` <- e[[paste0("a", ap, "_z0")]]
    `E(Y|a,M,1,W)` <- mu[[paste0("a", a, "_z1")]]
    `E(Y|a,M,0,W)` <- mu[[paste0("a", a, "_z0")]]

    H_Y11 <- ((Z == 1 & A == a) / `P(a'|W)`) * (`P(a'|M,1,W)` / `P(a|M,1,W)`)
    H_Y10 <- ((Z == 1 & A == a) / (`P(a'|W)` * `P(Z=0|a',W)`)) * 
      ((`P(a'|M,1,W)` / `P(a|M,1,W)`) * (`P(Z=0|M,a',W)` / `P(Z=1|M,a',W)`)) * 
      (`P(Z=1|a,W)` - `P(Z=1|a',W)`)
    H_Y00 <- ((Z == 0 & A == a) / `P(a|W)`) * (`P(a'|M,0,W)` / `P(a|M,0,W)`) * 
      (`P(Z=0|a,W)` / `P(Z=0|a',W)`)
    
    H_M11 <- (Z == 1 & A == ap) / `P(a'|W)`
    H_M10 <- ((Z == 0 & A == ap) / (`P(a'|W)` * `P(Z=0|a',W)`)) * (`P(Z=1|a,W)` - `P(Z=1|a',W)`)
    H_M00 <- ((Z == 0 & A == ap) / (`P(a'|W)` * `P(Z=0|a',W)`)) * `P(Z=0|a,W)`
    
    H_Z11 <- (A == ap) / `P(a'|W)` * rhos[[paste(a, ap, 1, 1, sep = "_")]]
    H_Z10 <- (((A == a) / `P(a|W)`) - ((A == ap) / `P(a'|W)`)) * 
      rhos[[paste(a, ap, 1, 0, sep = "_")]]
    H_Z00 <- -((A == a) / `P(a|W)`) * rhos[[paste(a, ap, 0, 0, sep = "_")]]
    
    H_W11 <- rhos[[paste(a, ap, 1, 1, sep = "_")]] * `P(Z=1|a',W)`
    H_W10 <- rhos[[paste(a, ap, 1, 0, sep = "_")]] * (`P(Z=1|a,W)` - `P(Z=1|a',W)`)
    H_W00 <- rhos[[paste(a, ap, 0, 0, sep = "_")]] * `P(Z=0|a,W)`
    
    eic_11 <- H_Y11 * (Y - mu$obs) + H_Z11 * (Z - q$obs) + 
      H_M11 * (`E(Y|a,M,1,W)` - rhos[[paste(a, ap, 1, 1, sep = "_")]]) + 
      H_W11
    
    eic_10 <- H_Y10 * (Y - mu$obs) + H_Z10 * (Z - q$obs) + 
      H_M10 * (`E(Y|a,M,1,W)` - rhos[[paste(a, ap, 1, 0, sep = "_")]]) + 
      H_W10
    
    eic_00 <- H_Y00 * (Y - mu$obs) + H_Z00 * (Z - q$obs) + 
      H_M00 * (`E(Y|a,M,0,W)` - rhos[[paste(a, ap, 0, 0, sep = "_")]]) + 
      H_W00
    
    eics[[paste(a, ap, sep = "_")]] <- eic_11 + eic_10 + eic_00
  }
  
  ans <- list(nie = mean(eics$`1_1`) - mean(eics$`1_0`), 
              nde = mean(eics$`1_0`) - mean(eics$`0_0`), 
              var_nie = var(eics$`1_1` - eics$`1_0`) / nrow(tmp), 
              var_nde = var(eics$`1_0` - eics$`0_0`) / nrow(tmp))
  ans$ci_nie <- ans$nie + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_nie)
  ans$ci_nde <- ans$nde + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_nde)
  ans
}
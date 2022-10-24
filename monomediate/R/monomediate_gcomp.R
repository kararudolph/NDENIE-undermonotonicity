#' Parametric mediation estimator under monotonicity assumption
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} containing all necessary variables
#'  for the estimation problem.
#' @param npsem \[\code{R6(Npsem)}\]\cr
#'  An \code{Npsem} object 
#' @param y_form \[\code{formula}\]\cr
#'  A formula for the outcome model
#' @param z_form \[\code{formula}\]\cr
#'  A formula for the "Z" model
#' @param m_form \[\code{formula}\]\cr
#'  A \bold{RIGHT-HAND SIDE ONLY} formula for the mediator model
#' @param y_family \[\code{character(1)}\]\cr
#'  Outcome variable type (i.e., "gaussian", "binomial").
#'
#' @return Natural (in)direct effect estimates
#' @export
#'
#' @examples
#' npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", c("m1", "m2"), "y")
#' y_form <- y ~ w1 + w2 + w3 + z + m1 + m2
#' z_form <- z ~ w1 + w2 + w3 + a
#' m_form <- ~ w3 + z
#' monomediate_gcomp(multiple_m, npsem, y_form, z_form, m_form, "binomial")
monomediate_gcomp <- function(data, npsem, y_form, z_form, m_form, y_family) {
  tmp <- data
  y_fit <- glm(y_form, family = y_family, data = tmp)
  z_fit <- glm(z_form, family = "binomial", data = tmp)

  q <- list()
  q$z1_a1 <- predict(z_fit, dplyr::mutate(tmp, "{npsem$A}" := 1), type = "response")
  q$z0_a1 <- 1 - q$z1_a1
  q$z1_a0 <- predict(z_fit, dplyr::mutate(tmp, "{npsem$A}" := 0), type = "response")
  q$z0_a0 <- 1 - q$z1_a0
  
  tmp$a1_z1 <- predict(y_fit, dplyr::mutate(tmp, "{npsem$A}" := 1, "{npsem$Z}" := 1), type = "response")
  tmp$a1_z0 <- predict(y_fit, dplyr::mutate(tmp, "{npsem$A}" := 1, "{npsem$Z}" := 0), type = "response")
  tmp$a0_z1 <- predict(y_fit, dplyr::mutate(tmp, "{npsem$A}" := 0, "{npsem$Z}" := 1), type = "response")
  tmp$a0_z0 <- predict(y_fit, dplyr::mutate(tmp, "{npsem$A}" := 0, "{npsem$Z}" := 0), type = "response")
  
  rhos <- list()
  for (as in list(c(1, 1), c(1, 0), c(0, 0))) {
    for (zs in list(c(1, 1), c(1, 0), c(0, 0))) {
      a <- as[1]
      ap <- as[2]
      z <- zs[1]
      zp <- zs[2]
      
      m_form_new <- update.formula(as.formula(paste0(paste0("a", a, "_z", z), "~ .")), m_form)
      rho_fit <- glm(m_form_new, family = "gaussian", data = tmp)
      
      rhos[[paste(a, ap, z, zp, sep = "_")]] <- 
        predict(rho_fit, dplyr::mutate(tmp, "{npsem$A}" := ap, "{npsem$Z}" := zp), type = "response")
    }
  }
  
  Z <- tmp[[npsem$Z]]
  A <- tmp[[npsem$A]]
  Y <- tmp[[npsem$Y]]

  comp <- list()
  for (as in list(c(1, 1), c(1, 0), c(0, 0))) {
    a <- as[1]
    ap <- as[2]

    `P(Z=0|a,W)` <- q[[paste0("z0_a", a)]]
    `P(Z=1|a',W)` <- q[[paste0("z1_a", ap)]]
    `P(Z=1|a,W)` <- q[[paste0("z1_a", a)]]
    
    H_W11 <- rhos[[paste(a, ap, 1, 1, sep = "_")]] * `P(Z=1|a',W)`
    H_W10 <- rhos[[paste(a, ap, 1, 0, sep = "_")]] * (`P(Z=1|a,W)` - `P(Z=1|a',W)`)
    H_W00 <- rhos[[paste(a, ap, 0, 0, sep = "_")]] * `P(Z=0|a,W)`
    
    comp[[paste(a, ap, sep = "_")]] <- H_W11 + H_W10 + H_W00
  }
  
  list(nie = mean(comp$`1_1`) - mean(comp$`1_0`), 
       nde = mean(comp$`1_0`) - mean(comp$`0_0`))
}

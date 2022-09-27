mediation_monotinicity <- function(data, a, z, m, y, w, family) {
  y_form <- as.formula(paste(y, paste(paste(paste(w, collapse = "+"), a, z, m, sep = "+"), 
                                      paste(z, m, sep = "*"), sep = "+"), 
                             sep = "~"))
  
  m_form <- reformulate(c(w, a, z), response = m)
  z_form <- reformulate(c(w, a), response = z)
  
  y_fit <- glm(y_form, data = data, family = family)
  m_fit <- glm(m_form, data = data)
  z_fit <- glm(z_form, data = data, family = "binomial")
  
  data_1 <- data_0 <- data
  data_1[[a]] <- 1
  data_0[[a]] <- 0
  
  w_c <- predict(z_fit, data_1, type = "response") - 
    predict(z_fit, data_0, type = "response")
  
  alpha_e <- coef(y_fit)[a]
  alpha_n <- coef(y_fit)[z]
  theta_0 <- coef(m_fit)["(Intercept)"]
  theta_c <- 0
  alpha_mn <- coef(y_fit)[paste(z, m, sep = ":")]
  
  alpha_e + mean((alpha_mn*theta_0 + alpha_mn*theta_c*0 + alpha_n)*w_c)
}

mediation_monotinicity <- function(data, a, z_fit, m_fit, y_fit, w) {
  data_1 <- data_0 <- data
  data_1[[a]] <- 1
  data_0[[a]] <- 0
  
  w_c <- predict(z_fit, data_1, type = "response") - 
    predict(z_fit, data_0, type = "response")
  
  alpha_e <- coef(y_fit)[a]
  alpha_n <- coef(y_fit)[z]
  theta_0 <- coef(m_fit)["(Intercept)"]
  theta_c <- 0
  alpha_mn <- coef(y_fit)[paste(z, m, sep = ":")]
  
  alpha_e + mean((alpha_mn*theta_0 + alpha_mn*theta_c*0 + alpha_n)*w_c)
}



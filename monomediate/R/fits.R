cv_g <- function(data, folds, npsem, learners) {
  fit <- function(fold, data, npsem, learners) {
    train <- origami::training(data) 
    valid <- origami::validation(data)
    
    task <- sl3::sl3_Task$new(
      data = train,
      covariates = npsem$W,
      outcome = npsem$A,
      outcome_type = "binomial"
    )
    
    est <- learners$train(task)
    
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = npsem$W,
      outcome = npsem$A,
      outcome_type = "binomial"
    )
    
    list(a1 = est$predict(task), 
         a0 = 1 - est$predict(task))
  }
  
  ans <- origami::cross_validate(fit,
                                 folds = folds,
                                 data = data,
                                 npsem = npsem,
                                 learners = learners,
                                 use_future = FALSE)
  ans$errors <- NULL
  i <- Reduce(c, lapply(folds, function(x) x$validation_set))
  lapply(ans, function(x) x[order(i)])
}

cv_e <- function(data, folds, npsem, learners) {
  fit <- function(fold, data, npsem, learners) {
    train <- origami::training(data) 
    valid <- origami::validation(data)
    
    task <- sl3::sl3_Task$new(
      data = train,
      covariates = c(npsem$M, npsem$Z, npsem$W),
      outcome = npsem$A,
      outcome_type = "binomial"
    )
    
    est <- learners$train(task)
    
    valid[, npsem$Z] <- 1
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$Z, npsem$W),
      outcome = npsem$A,
      outcome_type = "binomial"
    )
    
    ans <- list(a1_z1 = est$predict(task), 
                a0_z1 = 1 - est$predict(task))
    
    valid[, npsem$Z] <- 0
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$Z, npsem$W),
      outcome = npsem$A,
      outcome_type = "binomial"
    )
    
    ans$a1_z0 <- est$predict(task)
    ans$a0_z0 <- 1 - ans$a1_z0
    ans
  }
  
  ans <- origami::cross_validate(fit,
                                 folds = folds,
                                 data = data,
                                 npsem = npsem,
                                 learners = learners,
                                 use_future = FALSE)
  ans$errors <- NULL
  i <- Reduce(c, lapply(folds, function(x) x$validation_set))
  lapply(ans, function(x) x[order(i)])
}

cv_q <- function(data, folds, npsem, learners) {
  fit <- function(fold, data, npsem, learners) {
    train <- origami::training(data) 
    valid <- origami::validation(data)
    
    task <- sl3::sl3_Task$new(
      data = train,
      covariates = c(npsem$A, npsem$W),
      outcome = npsem$Z,
      outcome_type = "binomial"
    )
    
    est <- learners$train(task)
    
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$A, npsem$W),
      outcome = npsem$Z,
      outcome_type = "binomial"
    )
    
    ans <- list(obs = est$predict(task))
    
    valid[, npsem$A] <- 1
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$A, npsem$W),
      outcome = npsem$Z,
      outcome_type = "binomial"
    )
    
    ans$z1_a1 <- est$predict(task)
    ans$z0_a1 = 1 - ans$z1_a1
    
    valid[, npsem$A] <- 0
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$A, npsem$W),
      outcome = npsem$Z,
      outcome_type = "binomial"
    )
    
    ans$z1_a0 <- est$predict(task)
    ans$z0_a0 <- 1 - ans$z1_a0
    ans
  }
  
  ans <- origami::cross_validate(fit,
                                 folds = folds,
                                 data = data,
                                 npsem = npsem,
                                 learners = learners,
                                 use_future = FALSE)
  ans$errors <- NULL
  i <- Reduce(c, lapply(folds, function(x) x$validation_set))
  lapply(ans, function(x) x[order(i)])
}

cv_r <- function(data, folds, npsem, learners) {
  fit <- function(fold, data, npsem, learners) {
    train <- origami::training(data) 
    valid <- origami::validation(data)
    
    task <- sl3::sl3_Task$new(
      data = train,
      covariates = c(npsem$M, npsem$A, npsem$W),
      outcome = npsem$Z,
      outcome_type = "binomial"
    )
    
    est <- learners$train(task)
    
    valid[, npsem$A] <- 1
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$A, npsem$W),
      outcome = npsem$Z,
      outcome_type = "binomial"
    )
    
    ans <- list(z1_a1 = est$predict(task), 
                z0_a1 = 1 - est$predict(task))
    
    valid[, npsem$A] <- 0
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$A, npsem$W),
      outcome = npsem$Z,
      outcome_type = "binomial"
    )
    
    ans$z1_a0 <- est$predict(task)
    ans$z0_a0 <- 1 - ans$z1_a0
    ans
  }
  
  ans <- origami::cross_validate(fit,
                                 folds = folds,
                                 data = data,
                                 npsem = npsem,
                                 learners = learners,
                                 use_future = FALSE)
  ans$errors <- NULL
  i <- Reduce(c, lapply(folds, function(x) x$validation_set))
  lapply(ans, function(x) x[order(i)])
}

cv_mu <- function(data, folds, npsem, learners, type) {
  fit <- function(fold, data, npsem, learners) {
    train <- origami::training(data) 
    valid <- origami::validation(data)
    
    task <- sl3::sl3_Task$new(
      data = train,
      covariates = c(npsem$M, npsem$Z, npsem$A, npsem$W),
      outcome = npsem$Y,
      outcome_type = type
    )
    
    est <- learners$train(task)
    
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$Z, npsem$A, npsem$W),
      outcome = npsem$Y,
      outcome_type = type
    )
    
    ans <- list(obs = est$predict(task))
    
    valid[, npsem$A] <- 1
    valid[, npsem$Z] <- 1
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$Z, npsem$A, npsem$W),
      outcome = npsem$Y,
      outcome_type = type
    )
    
    ans$a1_z1 <- est$predict(task)
    
    valid[, npsem$A] <- 0
    valid[, npsem$Z] <- 1
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$Z, npsem$A, npsem$W),
      outcome = npsem$Y,
      outcome_type = type
    )
    
    ans$a0_z1 <- est$predict(task)
    
    valid[, npsem$A] <- 1
    valid[, npsem$Z] <- 0
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$Z, npsem$A, npsem$W),
      outcome = npsem$Y,
      outcome_type = type
    )
    
    ans$a1_z0 <- est$predict(task)
    
    valid[, npsem$A] <- 0
    valid[, npsem$Z] <- 0
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$M, npsem$Z, npsem$A, npsem$W),
      outcome = npsem$Y,
      outcome_type = type
    )
    
    ans$a0_z0 <- est$predict(task)
    ans
  }
  
  ans <- origami::cross_validate(fit,
                                 folds = folds,
                                 data = data,
                                 npsem = npsem,
                                 learners = learners,
                                 use_future = FALSE)
  ans$errors <- NULL
  i <- Reduce(c, lapply(folds, function(x) x$validation_set))
  lapply(ans, function(x) x[order(i)])
}

cv_rho <- function(data, folds, npsem, learners, a, ap, z, zp) {
  fit <- function(fold, data, npsem, learners, a, ap, z, zp) {
    train <- origami::training(data) 
    valid <- origami::validation(data)
    
    task <- sl3::sl3_Task$new(
      data = train,
      covariates = c(npsem$A, npsem$Z, npsem$W),
      outcome = paste0("a", a, "_z", z),
      outcome_type = "continuous"
    )
    
    est <- learners$train(task)
    
    valid[, npsem$Z] <- zp
    valid[, npsem$A] <- ap
    task <- sl3::sl3_Task$new(
      data = valid,
      covariates = c(npsem$A, npsem$Z, npsem$W),
      outcome = paste0("a", a, "_z", z),
      outcome_type = "continuous"
    )
    
    list(res = est$predict(task))
  }
  
  ans <- origami::cross_validate(fit,
                                 folds = folds,
                                 data = data,
                                 npsem = npsem,
                                 learners = learners,
                                 a = a, ap = ap, 
                                 z = z, zp = zp,
                                 use_future = FALSE)$res
  i <- Reduce(c, lapply(folds, function(x) x$validation_set))
  ans[order(i)]
}

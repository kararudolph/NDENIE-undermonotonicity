#' @export
Npsem <- R6::R6Class(
  "Npsem",
  public = list(
    W = NULL,
    A = NULL,
    Z = NULL,
    M = NULL,
    Y = NULL,
    initialize = function(W, A, Z, M, Y) {
      self$W <- W
      self$A <- A
      self$Z <- Z
      self$M <- M
      self$Y <- Y
    },
    get = function(.data, var = c("W", "A", "Z", "M", "Y")) {
      .data[[self[[match.arg(var)]]]]
    },
    dt = function(.data) {
      data.table::as.data.table(.data[, c(self$W, self$A, self$Z, self$M, self$Y)])
    }
  )
)
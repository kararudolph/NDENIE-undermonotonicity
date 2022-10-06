Lrnr_constant <- R6::R6Class(
  classname = "Lrnr_constant", inherit = sl3::Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(val = 0.5, ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),
    
    .train = function(task) {
      args <- self$params
      fit_object <- args$val
      return(fit_object)
    },
    .predict = function(task = NULL) {
      rep(self$fit_object, nrow(task$X))
    }
  )
)

#' Validate inputs for MPM functions
#'
#' Checks that all arguments passed to the main MPM functions are of the
#' correct type and structure before any computation begins.
#'
#' @param fys Numeric scalar. First-year survival (age 0 to age 1).
#'   In a pre-breeding census model this is folded into the fertility term.
#' @param js Data frame of juvenile survival rates, one row per juvenile
#'   stage. Either one column (`mean`) for constant survival, or two columns
#'   (`mean`, `min`) when a survival ramp is requested.
#' @param as Numeric scalar. Adult (final-stage) survival rate.
#' @param f Numeric scalar. Annual fecundity — female offspring produced per
#'   adult female per year. Does **not** include first-year survival adjustment
#'   (that is handled via `fys`).
#' @param t Data frame of stage durations, one row per juvenile stage. Either
#'   one column (`dur`) for fixed duration, or three columns (`dur`, `min`,
#'   `max`) for variable durations.
#'
#' @return Invisibly returns `NULL`. Called for its side-effect of stopping
#'   with an informative error message if any input is invalid.
#'
#' @keywords internal

init_input_check <- function(fys,js,as,f,t){
  if (length(fys)>1) {
    stop("first-year survival must be a single value")
  }
  if (!is.data.frame(js)) {
    stop("Argument 'js' (juv surv) must be data frame")
  }
  if (!all(sapply(js,is.numeric))){
    stop("All columns of data frame 'js' must be numeric")
  }
  if (length(as)>1) {
    stop("adult survival must be a single value")
  }
  if (!is.numeric(as)){
    stop("adult survival must be numeric")
  }
  if (length(f)>1) {
    stop("fecundity must be a single value")
  }
  if (!is.numeric(f)){
    stop("fecundity must be numeric")
  }
  if (!is.data.frame(t)) {
    stop("Argument 't' must be a data frame")
  }
  if (!all(sapply(t,is.numeric))){
    stop("All columns of data frame 't' must be numeric")
  }
  if (!ncol(t)%in%c(1,3)){
    stop("Data frame 't' have either 1 column or 3 columns (mean, min, and max duration)")
  }
}

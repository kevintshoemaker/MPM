#' Build a scenario list for MPM functions
#'
#' A convenience wrapper that packages life-history parameters into a named
#' list accepted by [do_aas()], [do_fas()], and [do_unroll()].
#'
#' @param fysurv Numeric scalar. First-year survival probability.
#' @param jsurv Either a numeric vector of mean juvenile survival rates (one
#'   per stage) or a data frame with columns `mean` and optionally `min` when
#'   a survival ramp is desired.
#' @param asurv Numeric scalar. Adult survival probability.
#' @param fec Numeric scalar. Annual fecundity (female offspring per adult
#'   female per year, excluding first-year survival).
#' @param dur Either a numeric vector of stage durations (one per juvenile
#'   stage) or a data frame with columns `dur`, `min`, and `max` for variable
#'   stage durations and number of rows equal to the number of stages.
#' @param ramp Logical. If `TRUE`, juvenile survival increases monotonically
#'   within each stage (a "ramp"), constrained to match the specified mean.
#'   Only used when `jsurv` is supplied as a vector; if `jsurv` is already a
#'   data frame this argument is ignored.
#'
#' @return A named list with elements `fysurv`, `jsurv`, `asurv`, `fec`,
#'   and `dur`.
#'
#' @examples
#' # Single juvenile stage, fixed duration, no ramp
#' scen <- gen_scen(
#'   fysurv = 0.45,
#'   jsurv  = 0.75,
#'   asurv  = 0.96,
#'   fec    = 1.5,
#'   dur    = 9,
#'   ramp   = FALSE
#' )
#'
#' # Two juvenile stages with survival ramp
#' scen2 <- gen_scen(
#'   fysurv = 0.45,
#'   jsurv  = data.frame(mean = c(0.5, 0.8), min = c(0.25, 0.65)),
#'   asurv  = 0.95,
#'   fec    = 2.29,
#'   dur    = c(3, 4),
#'   ramp   = TRUE
#' )
#'
#' @export

gen_scen <- function(fysurv,jsurv,asurv,fec,dur,ramp){
  if(is.vector(dur)) dur = data.frame(dur=dur)
  nst = length(jsurv) + 1
  if(is.vector(jsurv)){
    if(ramp){
      jsurv = data.frame(mean=jsurv,min=as.numeric(NA))   # the NA means that the ramp by default starts at the previous max value
    }else{
      jsurv = data.frame(mean=jsurv)
    }
  }
  scen <- list()
  scen$fysurv = fysurv
  scen$jsurv = jsurv
  scen$asurv = asurv
  scen$fec = fec  # annual fecundity (assume pre-breeding census- this represents the number of new one-year-olds entering the population next year)
  scen$dur = dur  # duration of juvenile stage (mean, min, max)
  scen
}





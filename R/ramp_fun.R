#' Generate a within-stage survival ramp
#'
#' Computes annual survival rates across the years of a single stage such
#' that survival increases (or decreases) monotonically from a starting value
#' to an ending value while matching a specified mean. The shape follows a
#' scaled exponential curve whose curvature is found by one-dimensional
#' optimization.
#'
#' @param sm Numeric scalar. Target **mean** survival across the included
#'   years of the stage.
#' @param sa Numeric scalar. **Starting** (minimum) survival value at the
#'   beginning of the stage.
#' @param sz Numeric scalar. **Ending** (maximum) survival value at the end
#'   of the stage (typically the minimum survival of the next stage).
#' @param y Integer. Total number of years in the stage.
#' @param include_first Logical. Whether to include the first year when
#'   computing and constraining the mean. Default `TRUE`.
#' @param include_last Logical. Whether to include the last year when
#'   computing and constraining the mean. Default `FALSE` (the last value
#'   belongs to the next stage).
#'
#' @return Numeric vector of length `sum(include_first, rep(TRUE, y-2), include_last)`
#'   giving annual survival rates for the included years of the stage.
#'
#' @examples
#' # Survival ramp from 0.5 to 0.95 over 10 years, mean ~0.8
#' r <- ramp_fun(sm=0.8, sa=0.5, sz=0.95, y=10)
#' plot(r, type="l", ylab="Survival", xlab="Year within stage")
#' mean(r)   # should be close to 0.8
#'
#' @export
ramp_fun = function(sm,sa,sz,y,include_first=T,include_last=F){
  include = rep(T,y)
  include[1] = include_first; include[y]=include_last
  opt_func = function(k){
    thisseq = sa + (sz-sa) * ( (1-exp(-k * (seq(1,y,length=y)-1) ) ) / (1-exp(-k * (y - 1) ) ) )  # sequence of survival in each year
    msurv = mean(thisseq[include])
    (msurv-sm)^2
  }
  opt = optimize(opt_func,c(-10,10))
  thisk = opt$minimum
  (sa + (sz-sa) * ( (1-exp(-thisk * ( 1:y - 1 ) ) ) / (1 - exp(-thisk * (y - 1) ) ) ))[include]
}

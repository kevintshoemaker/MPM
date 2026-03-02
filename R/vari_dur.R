#' Negative entropy (objective function for max-entropy optimization)
#'
#' @param p Numeric vector of probabilities summing to 1.
#' @return Numeric scalar (negative Shannon entropy).
#' @keywords internal
negent <- function(p) {
  sum(p * log(p))
}

#' Generate equality constraint function for stage duration optimization
#'
#' Returns a function that evaluates three equality constraints for use with
#' [Rsolnp::solnp()]: (1) probabilities sum to 1, (2) mean matches target,
#' (3) variance matches a default spread based on the range of ages.
#'
#' @param ages Integer vector of possible stage-exit ages.
#' @param mage Numeric scalar. Target mean age at stage exit.
#'
#' @return A function `f(p)` returning a length-3 numeric vector of
#'   constraint residuals.
#'
#' @keywords internal
generate_eq = function(ages, mage){
  sig=(max(ages)-min(ages))/5
  function(p){
    z1 = sum(p)-1       # constraint 1: sum to 1
    thismean = sum(p*ages)
    z2 = thismean - mage   # constraint 2: match user defined mean
    this2cm = sum((ages-mage)^2 * p)
    z3 = this2cm-sig^2  # constraint 3: match user defined second central moment
    c(z1,z2,z3)  # z4
  }
}

#' Inequality constraint for log-concavity
#'
#' Ensures the second difference of `log(p)` is non-positive, producing a
#' unimodal (log-concave) distribution.
#'
#' @param p Numeric vector of probabilities.
#' @return Numeric vector of second differences of `log(p)`.
#' @keywords internal
ineq <- function(p) {
  diff(log(p), differences=2)
}

#' Compute maximum-entropy distribution for variable stage duration
#'
#' For a stage with a mean, minimum, and maximum duration, finds the
#' probability distribution over integer ages at stage exit that maximizes
#' entropy subject to matching the mean and a variance that approximately matches
#' the min/max constraints. Results are
#' cached with [memoise::memoise()] to avoid redundant optimization.
#'
#' @param meand Numeric. Mean stage duration (years).
#' @param mind Numeric. Minimum stage duration (years).
#' @param maxd Numeric. Maximum stage duration (years).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{`dur`}{Integer stage-exit ages considered.}
#'     \item{`prob`}{Probability of exiting the stage at each age.}
#'   }
#'
#' @examples
#' d <- vari_dur(meand=9, mind=5, maxd=12)
#' plot(d$dur, d$prob, type="h", xlab="Age at exit", ylab="Probability")
#'
#' @export
vari_dur <- function(meand, mind, maxd){
  dur <- floor(mind):ceiling(maxd)
  p0 = dnorm(dur,meand,(maxd-mind)/4) ; p0= p0/sum(p0)
  bounds = cbind(rep(1e-9,length(dur)),rep(1-1e-9,length(dur)))
  eq = generate_eq(dur,meand)
  s = Rsolnp::solnp(pars=p0, fun = negent,
            eqfun = eq, eqB = c(0,0,0),
            ineqfun = ineq, ineqLB = rep(-1e9,length(dur)-2), ineqUB=rep(0,length(dur)-2),
            LB=bounds[,1],UB=bounds[,2], control = list(trace=0))   #
  data.frame(
    dur=dur,
    prob=s$pars
    # cprob=pmax(0,pmin(1,c(s$pars[1], s$pars[2:length(dur)] / (1-cumsum(s$pars))[1:(length(dur)-1)]) ))
  )
}
vari_dur = memoise::memoize(vari_dur)

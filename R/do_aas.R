#' Compute stage transition probabilities using the
#' “asymptotic age-within-stage structure” (AAS) method of
#' Kendall et al. (2019)
#'
#' Internal helper. Computes `gamma` (the fraction of survivors that
#' transition out of a stage) using the 'AAS' method of Kendall et al. (2019).
#'
#' @param s Numeric scalar. Survival rate for the focal stage.
#' @param l Numeric scalar. Current estimate of the population growth rate
#'   (lambda).
#' @param t Numeric scalar. Duration of the focal stage (years).
#'
#' @return Numeric scalar. The fraction of survivors that transition to the
#'   next stage.
#'
#' @references
#' Kendall, B.E., Fujiwara, M., Diaz-Lopez, J., Schneider, S., Voigt, J. and Wiesner, S., 2019.
#' Persistent problems in the construction of matrix population models.
#' *Ecological modelling*, 406, pp.33-43.
#'
#' @keywords internal

gamma_aas <- function(s,l,t){
  sl = s/l
  if(sl==1) 0.5 else (sl^(t-1)*(1-sl)) / (1-sl^t)  # note that the formula involves the sum of finite geometric series so can be simplified
}

#' Stage transition sub-diagonal entries (AAS method)
#'
#' @param s Numeric scalar. Survival of the focal stage.
#' @param l Numeric scalar. Current lambda estimate.
#' @param t Numeric scalar. Duration of the focal stage.
#'
#' @return Numeric vector of length 2: `c(stay, transition)` probabilities
#'   (both multiplied by `s`).
#'
#' @keywords internal
thistrans_aas <- function(s,l,t){
  gam = gamma_aas(s,l,t)   # compute gamma (fraction of survivors transitioning)
  c(s*(1-gam),s*gam)    # allocate the transitions
}

#' Construct a stage-based MPM using the “asymptotic age-within-stage structure” (AAS) method of
#' Kendall et al. (2019)
#'
#' Builds a stage-structured population matrix using the “asymptotic age-within-stage structure” (AAS) method of
#' Kendall et al. (2019). Lambda is solved iteratively. The method
#' assumes constant survival within each stage and exact (fixed) stage
#' durations; a warning is issued if variable durations or a survival ramp
#' are supplied.
#'
#' @inheritParams init_input_check
#'
#' @return A square numeric matrix of dimension `(number of stages) x
#'   (number of stages)`. The dominant eigenvalue of this matrix equals the
#'   asymptotic population growth rate lambda.
#'
#' @references
#' Kendall, B.E., Fujiwara, M., Diaz-Lopez, J., Schneider, S., Voigt, J. and Wiesner, S., 2019.
#' Persistent problems in the construction of matrix population models.
#' *Ecological modelling*, 406, pp.33-43.
#'
#' @seealso [do_unroll()], [do_fas()]
#'
#' @examples
#' scen <- gen_scen(fysurv=0.45, jsurv=0.75, asurv=0.96, fec=1.5,
#'                  dur=9, ramp=FALSE)
#' mat <- do_aas(scen$fysurv, scen$jsurv, scen$asurv, scen$fec, scen$dur)
#' popbio::lambda(mat)
#'
#' @export
do_aas <- function(fys,js,as,f,t){
  init_input_check(fys,js,as,f,t)
  ss = c(js[[1]],as)   # join juv and adult survival rates
  n=length(ss); m0 <- matrix(0,n,n) ; m0[n,n] <- ss[n] ;  m0[1,n] <- f*fys    # construct init matrix
  lam = 1  # take an initial guess at lambda
  dif=Inf ; tol=1e-6
  while(dif>tol){
    l = popbio::lambda(m0)
    for(g in 1:(n-1)) m0[g:(g+1),g] = thistrans_aas(ss[g],l,t[[1]][g])
    dif = abs(popbio::lambda(m0) - l)
  }
  if(ncol(t)>1) warning("this method is only valid for exact stage durations: variable stage durations not supported. Only first column of input 't' was used in this analysis")
  if(ncol(js)>1) warning("this method is only valid for constant survival within stage: survival 'ramp' not supported. Only first column of input 'js' was used in this analysis")
  m0
}


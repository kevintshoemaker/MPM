#' Construct a stage-based MPM using the "Flat Age-Within-Stage" (FAS) method (not recommended)
#'
#' Implements the "Fixed Age at Stage" (FAS) approximation, in which exactly
#' `1/t` of each stage is assumed to transition each year regardless of
#' lambda. This method is known to be incorrect in most cases and is included
#' for comparison purposes only.
#'
#' @inheritParams init_input_check
#'
#' @return A square numeric stage-based matrix.
#'
#' @section Warning:
#' This method systematically overestimates population growth rates and
#' **should not be used** for management or conservation applications. It is
#' retained here for benchmarking against legacy analyses (e.g., the Gopher
#' Tortoise SSA model).
#' @references
#' Kendall, B.E., Fujiwara, M., Diaz-Lopez, J., Schneider, S., Voigt, J. and Wiesner, S., 2019.
#' Persistent problems in the construction of matrix population models.
#' *Ecological modelling*, 406, pp.33-43.
#'
#' @seealso [do_aas()], [do_unroll()]
#'
#' @examples
#' scen <- gen_scen(fysurv=0.45, jsurv=0.75, asurv=0.96, fec=1.5,
#'                  dur=9, ramp=FALSE)
#' mat <- do_fas(scen$fysurv, scen$jsurv, scen$asurv, scen$fec, scen$dur)
#' popbio::lambda(mat)   # compare with do_aas()
#'
#' @export
do_fas <- function(fys,js,as,f,t){
  init_input_check(fys,js,as,f,t)
  ss = c(js[[1]],as)
  n=length(ss); m0 <- matrix(0,n,n) ; m0[n,n] <- ss[n] ;  m0[1,n] <- f*fys    # construct init matrix
  gam <- 1/t[[1]]    # fraction of stage that is in the transition age
  for(g in 1:(n-1)) m0[g:(g+1),g] <- c(ss[g]*(1-gam[g]), ss[g]*gam[g] )
  warning("this method is incorrect in most cases and can be highly misleading in many real-world species")
  if(ncol(t)>1) warning("variable stage durations not supported. Only first column of input 't' was used in this analysis")
  if(ncol(js)>1) warning("this method is only valid for constant survival within stage: survival 'ramp' not supported. Only first column of input 'js' was used in this analysis")
  m0
}

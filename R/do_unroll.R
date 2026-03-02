#' Construct an age-based MPM by unrolling stage durations
#'
#' Converts life-history parameters into a full age-based (Leslie-type)
#' matrix by "unrolling" each juvenile stage into its constituent annual age
#' classes. Supports fixed or variable stage durations and optional
#' within-stage survival ramps. This is generally the recommended method for
#' modeling any age-structured population.
#'
#' @inheritParams init_input_check
#'
#' @details
#' **Pre-breeding census model.** All fertility terms implicitly assume that
#' the census occurs just before the breeding event, so the fertility
#' row of the matrix is `f * fys` (fecundity × first-year survival). The
#' matrix dimension equals `sum(stage durations) + 1` (the `+1` is for the
#' adult stage).
#'
#' **Variable stage durations.** When `t` has three columns (`dur`, `min`,
#' `max`), the probability distribution of ages at stage transition is
#' computed via maximum-entropy optimization ([vari_dur()]) and convolved
#' across stages.
#'
#' **Survival ramp.** When `js` has a `min` column, annual survival within
#' each stage follows a monotonically increasing or decreasing curve
#' (computed by [ramp_fun()]) that starts at the specified minimum (or the
#' previous stage's maximum if `min = NA`) and ends at the minimum survival
#' of the subsequent stage, while matching the specified mean.
#'
#' @return A square numeric age-based matrix of dimension
#'   `(total ages) x (total ages)`. The dominant eigenvalue equals lambda.
#' @references
#' Kendall, B.E., Fujiwara, M., Diaz-Lopez, J., Schneider, S., Voigt, J. and Wiesner, S., 2019.
#' Persistent problems in the construction of matrix population models.
#' *Ecological modelling*, 406, pp.33-43.
#' @seealso [do_aas()], [do_fas()], [ramp_fun()], [vari_dur()]
#'
#' @examples
#' # Fixed duration, no ramp
#' scen <- gen_scen(fysurv=0.45, jsurv=0.75, asurv=0.96,
#'                  fec=1.5, dur=9, ramp=FALSE)
#' mat <- do_unroll(scen$fysurv, scen$jsurv, scen$asurv, scen$fec, scen$dur)
#' popbio::lambda(mat)
#'
#' # Variable stage duration with survival ramp
#' scen7 <- gen_scen(
#'   fysurv = 0.45,
#'   jsurv  = data.frame(mean = 0.75, min = 0.5),
#'   asurv  = 0.95,
#'   fec    = 2.29,
#'   dur    = data.frame(dur=9, min=6, max=13),
#'   ramp   = TRUE
#' )
#' mat7 <- do_unroll(scen7$fysurv, scen7$jsurv, scen7$asurv,
#'                   scen7$fec, scen7$dur)
#' popbio::lambda(mat7)
#'
#' @export
do_unroll <- function(fys,js,as,f,t){
  init_input_check(fys,js,as,f,t)
  na = sum(t[[ncol(t)]])+1; ns = nrow(t)   # na is number of age classes, ns is number of juv stage classes
  m0 <- matrix(0,na,na) ;  m0[na,na] <- as    # construct init matrix
  stnames = paste0("st_", 0:(ns+1))   # kts added zero age/stage for completeness [but actually probably should take it out again to simplify]
  survnames = paste0("surv_", 0:(ns+1))
  fystage = stnames[1]
  adstage = stnames[ns+2]
  juvstages = setdiff(stnames,c(fystage,adstage))
  agedf = data.frame(    # make data frame to do preliminary computations
    age = 0:na
  )
  agedf[stnames] <- lapply(1:(ns+1), function(x) 0)  # add cols for fraction in stage by age
  agedf[survnames] <- lapply(1:(ns+1), function(x) 0)  # add cols for survival of each stage by age
  agedf$st_0[1] = 1; agedf$st_1[2] = 1   # set all newborns in stage zero, and temporarily set all age 1 individuals in stage 1.
  agedf[[survnames[1]]] = fys     # set first and last survival rates
  agedf[[tail(survnames,1)]] = as

  s=1
  for(s in 1:ns){
    thisst = juvstages[s]; thissurv=gsub("st_","surv_",thisst)

    # fill in survival for juvenile stages
    if(ncol(js)>1){  # if survival ramp
      ss = c(fys,js$mean,as)   # put together all the survival rates for all stages
      names(ss) = stnames
      startage = ifelse(s==1,1,cumsum(t$dur[1:s])+1)
      theseages = startage:(startage+t$dur[s]-1)
      thisramp = ramp_fun(
        js$mean[s],
        ifelse(is.na(js$min[s]),ss[s],js$min[s]),
        ss[s+2],
        t$dur[s]+ifelse(is.na(js$min[s]),2,1),
        include_first = !is.na(js$min[s])
      )
      agedf[[thissurv]][agedf$age%in%theseages] = thisramp
      agedf[[thissurv]][agedf$age<min(theseages)] = thisramp[1]
      agedf[[thissurv]][agedf$age>max(theseages)] = tail(thisramp,1)
    }else{   # if not ramp
      agedf[[thissurv]]=js$mean[s]
    }
  }  # end loop through stages

  # fill in stage probabilities by age for juvenile stages
  if(ncol(t)>1){
    td = lapply(1:nrow(t),function(z) vari_dur(t[z,1],t[z,2],t[z,3] )   )  # variable ages at transition
    j=2
    for(j in 1:ns){   # determine probability densities of transitioning at each age given variable durations
      thiscol = stnames[1+j]
      nextcol = stnames[1+j+1]
      thistd = td[[j]]
      if(j==1){
        agedf[[nextcol]][match(thistd$dur,agedf$age-1)] = thistd$prob
      }else{
        nonzero_ages = agedf$age[agedf[[thiscol]]>1e-15 ]
        a=2
        for(a in 1:length(nonzero_ages)){  # loop through all ages of the prior stage and redistribute them
          thisa = nonzero_ages[a]
          agedf[[nextcol]][agedf$age%in%(thisa+thistd$dur)] = agedf[[nextcol]][agedf$age%in%(thisa+thistd$dur)] +
            agedf[[thiscol]][agedf$age==thisa] * thistd$prob
        }
      }
    }# end loop through stages
    colSums(agedf[,stnames])   # make sure cols sum to 1...

    stprob = agedf[,stnames]
    stprob$st_1[agedf$age>0] = 1  # start with everyone in first stage
    s=1
    for(s in 1:ns){
      start = stprob[,s+1]   # all who did not transition to the next stage stay
      pstay = pmax(0,pmin(1,1-cumsum(agedf[,sprintf("st_%s",s+1)])))
      stprob[,s+1] = start * pstay   # need to be in stage currently AND stay in stage
      stprob[,s+2] = start * (1 - pstay)
    }
    rowSums(stprob)   # make sure rows sum to 1
    agedf[,stnames] = stprob

  }else{  # if no variable-age stages
    prev=2; s=1
    for(s in 1:ns){
      agedf[prev:(prev+t$dur[s]-1),juvstages[s]] = 1
      prev=prev+t$dur[s]
    }
    agedf[[adstage]][nrow(agedf)] = 1
  } # end if no variable-age stages

  agedf$surv = rowSums( agedf[,survnames]*agedf[,stnames] )   # survival by age
  agedf$f = agedf[,adstage] * f * fys

  m1 <- diag(agedf$surv[agedf$age%in%c(1:(na-1))])
  m0[2:na,1:(na-1)] = m1
  m0[1,] = agedf$f[-1]
  m0
}


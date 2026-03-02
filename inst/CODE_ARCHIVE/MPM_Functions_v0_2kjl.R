
# FUNCTIONS for matrix population modeling in R ----------


## Initial error checking function ------------


# s is stage specific survival, l is lambda, T is focal stage duration

init_input_check <- function(js,as,f,t){
  if (!is.data.frame(js)) {
    stop("Argument 'js' (juv surv) must be data frame")
  }
  #kjl
  if (!all(sapply(js[,setdiff(names(js),"startbefore")],is.numeric))){
    stop("All columns (except startbefore) of data frame 'js' must be numeric")
  }
  #kjl
  if(length(js$startbefore)==nrow(js)) { #if we have a startbefore column...
    if (!all(sapply(js[,"startbefore"],is.logical))){ #and if it is not logical...
     stop("Column 'startbefore' of data frame 'js' must be logical if it exists")
    }
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

## Scenario generating function -------------

gen_scen <- function(jsurv,asurv,fec,dur){
  if(is.vector(dur)) dur = data.frame(dur=dur)
  nst = length(jsurv) + 1
  if(is.vector(jsurv)) jsurv = data.frame(mean=jsurv)
  scen <- list()
  scen$jsurv = jsurv 
  scen$asurv = asurv 
  scen$fec = fec  # annual fecundity (assume pre-breeding census- this represents the number of new one-year-olds entering the population next year)  
  scen$dur = dur  # duration of juvenile stage (mean, min, max)
  scen
}

## Functions for AAS method ---------------

#  Assume pre-breeding census
#  s is vector of survival, f is fecundity (scalar), t is duration of juvenile stage
#  l is lambda, mat is a stage matrix...

gamma_aas <- function(s,l,t){         
  sl = s/l
  if(sl==1) 0.5 else (sl^(t-1)*(1-sl)) / (1-sl^t)  # note that the formula involves the sum of finite geometric series so can be simplified
}

# s is survival of focal stage, l is current lambda estimate, t is duration of focal stage
thistrans_aas <- function(s,l,t){
  gam = gamma_aas(s,l,t)   # compute gamma (fraction of survivors transitioning)
  c(s*(1-gam),s*gam)    # allocate the transitions
}

# s is survival of all stages, f is fecundity of final stage, t is duration of juvenile stages
do_aas <- function(js,as,f,t){
  init_input_check(js,as,f,t)
  ss = c(js[[1]],as)
  n=length(ss); m0 <- matrix(0,n,n) ; m0[n,n] <- ss[n] ;  m0[1,n] <- f    # construct init matrix
  lam = 1  # take an initial guess at lambda
  dif=Inf ; tol=1e-6
  while(dif>tol){
    l = popbio::lambda(m0)
    for(g in 1:(n-1)) m0[g:(g+1),g] = thistrans_aas(ss[g],l,t[[1]][g]) 
    dif = abs(lambda(m0) - l)
  }
  if(ncol(t)>1) warning("this method is only valid for exact stage durations: variable stage durations not supported. Only first column of input 't' was used in this analysis")
  if(ncol(js)>1) warning("this method is only valid for constant survival within stage: survival 'ramp' not supported. Only first column of input 'js' was used in this analysis")
  m0
}

## Function for incorrect "FAS" method ------------

# function for implementing the incorrect "FAS" method (used in Gopher Tortoise SSA model)
do_fas <- function(js,as,f,t){
  init_input_check(js,as,f,t)
  ss = c(js[[1]],as)
  n=length(ss); m0 <- matrix(0,n,n) ; m0[n,n] <- ss[n] ;  m0[1,n] <- f    # construct init matrix
  gam <- 1/t[[1]]    # fraction of stage that is in the transition age
  for(g in 1:(n-1)) m0[g:(g+1),g] <- c(ss[g]*(1-gam[g]), ss[g]*gam[g] )
  warning("this method is incorrect in most cases and can be highly misleading in many real-world species")
  if(ncol(t)>1) warning("variable stage durations not supported. Only first column of input 't' was used in this analysis")
  if(ncol(js)>1) warning("this method is only valid for constant survival within stage: survival 'ramp' not supported. Only first column of input 'js' was used in this analysis")
  m0
}

## Variable age at maturity functions ------------

### Function to compute the objective of max entropy (minimize the negative entropy) for discrete distribution
negent = function(p){
  sum(p*log(p))
}

### Closure for generating equality constraint function- assume equality constraints are zeros a.la lagrange multiplier method
# ages are the possible ages for this stage, 
# mean is the point estimate for maturation, 
# sd is the standard deviation from the mean,
# skew is the skewness (third central moment, standardized by dividing by sigma cubed)
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
# ages=1:10;mage=5.5    #;skew=0

ineq = function(p){
  diff(diff(log(p)))
}

# sum(s$pars)
# s$pars
# plot(s$pars)
# sum(ages*s$pars)
# sqrt(sum((ages-m)^2*s$pars))
# sum((ages-m)^3*s$pars)/(sig^3)

### Function for allocating stage durations in variable stage duration model
meand=9;mind=5;maxd=12
vari_dur <- function(meand, mind, maxd){
  dur <- floor(mind):ceiling(maxd)
  p0 = dnorm(dur,meand,(maxd-mind)/4) ; p0= p0/sum(p0)
  bounds = cbind(rep(1e-9,length(dur)),rep(1-1e-9,length(dur)))
  eq = generate_eq(dur,meand)
  s = solnp(pars=p0, fun = negent, 
            eqfun = eq, eqB = c(0,0,0),
            ineqfun = ineq, ineqLB = rep(-1e9,length(dur)-2), ineqUB=rep(0,length(dur)-2),
            LB=bounds[,1],UB=bounds[,2], control = list(trace=0))   # 
  data.frame(
    dur=dur,
    prob=s$pars
    # cprob=pmax(0,pmin(1,c(s$pars[1], s$pars[2:length(dur)] / (1-cumsum(s$pars))[1:(length(dur)-1)]) ))
  )
}

#vari_dur = memoize(vari_dur)

## Survival ramp functions ------------
# Survival ramp function (with optional “startbefore” hatchling year)
#
# sm = target MEAN survival across JUVENILE years returned (length y)
# sa = start survival:
#      - startbefore=FALSE: sa is FIRST juvenile-year survival (age 1–2)
#      - startbefore=TRUE : sa is HATCHLING survival (age 0–1) and will be dropped
# sz = survival at END of juvenile stage (approaching next stage survival)
# y  = number of juvenile years to RETURN (integer)
# startbefore = if TRUE, build ramp one year longer (include hatchling), then drop it
#
# Returns: always a vector of length y (juvenile survivals only)
#
# Adds a warning if the returned juvenile mean is not close to sm.
#
ramp_fun <- function(sm, sa, sz, y, startbefore = FALSE,
                     warn_tol = 0.01){
  
  # Internal ramp length:
  # - normal: y juvenile years
  # - startbefore: y juvenile years + 1 hatchling year at the front
  y_internal <- if (startbefore) y + 1 else y
  
  # Build internal ramp of length y_internal for a given k
  ramp_from_k <- function(k){
    sa + (sz - sa) *
      ((1 - exp(-k * ((1:y_internal) - 1))) /
         (1 - exp(-k * (y_internal - 1))))
  }
  
  # Extract juvenile years (the ones we want mean(sm) to apply to)
  juvenile_slice <- function(r){
    if (startbefore) r[-1] else r
  }
  
  # Objective: match mean over juvenile years only
  opt_func <- function(k){
    r <- ramp_from_k(k)
    juv <- juvenile_slice(r)
    abs(mean(juv) - sm)
  }
  
  # Optimize k to minimize absolute mean difference
  opt   <- optimize(opt_func, c(-10, 10))
  thisk <- opt$minimum
  
  # Build final ramp
  out_internal <- ramp_from_k(thisk)
  
  # Always return juvenile years only
  out <- juvenile_slice(out_internal)
  
  # Warn if mean isn't matching target (report magnitude)
  mean_out <- mean(out)
  diff     <- mean_out - sm
  if (abs(diff) > warn_tol) {
    warning(sprintf(
      "ramp_fun(): juvenile mean does not match target sm. mean(out)=%.10f, sm=%.10f, diff=%.10f (abs diff=%.3e).",
      mean_out, sm, diff, abs(diff)
    ), call. = FALSE)
  }
  
  out
}

#try some examples to verify that the startbefore logic works
if(F){
  x <- 1:7
  y <- ramp_fun(0.75, 0.13, 0.97, 7, startbefore = TRUE)
  mean(y)
  plot(x, y,
       type = "s",
       ylim = c(0, 1),
       xlab = "Age (years)",
       ylab = "Survival",
       main = "Juvenile survival ramp")
}




## Function for generating variable age at maturity from a given percentage of ma  ------------------
#make a dur dataframe vor variable age at maturity
#inputs: 
# - ma
# - the percentage of ma that represents the width of the range of possible ages at maturity
# - the hard minimum age for maturation (youngest ever...7yo from Jon Moore's data)

make_var_dur<-function(ma,range,hardmin=7){
  data.frame(dur=ma-1, 
             min=floor(max(ma-(range/2*ma),hardmin)-1),
             max=ceiling(ma+(range/2*ma)-1))
}

## Function for "UNROLL" method  ------------------
#example t with variable duration:
#dur min max
#   9   8  10

# KJL notes 1/27/26: 
#  - added annotated by chatgpt
#  - updated by me to include the new ramp function (see #KJL)
#  - updated by me to fix behavior when variable duration produces incorrect dimensions when filling in post-ramp surv (when t$max = t$dur+1) (see #KJL)
  
# Unroll a stage-based specification into an age-based (Leslie-like) matrix.
#
# Notes / assumptions:
#  - Pre-breeding census formulation.
#  - Fecundity term `f` is assumed to ALREADY include neonate / first-year survival.
#  - Survival terms represent survival while in each stage (juvenile stages + adult).
#  - Supports:
#      (A) fixed stage durations (t is a vector / 1-col object)
#      (B) variable stage durations (t has dur/min/max columns)
#      (C) optional survival “ramps” within juvenile stages (js has mean + min)
#
#js=data.frame(mean=0.75,min=0.13,startbefore=T); as=0.97; f=0.1; t=make_var_dur(ma=8,range=0.2,hardmin = 7)

# example problem case: 
# t=make_var_dur(ma=8,range=0.2,hardmin = 7)
# dur min max
#    7   6   8


do_unroll <- function(js, as, f, t){
  
  #----------------------------
  # 1) Basic input checks
  #----------------------------
  init_input_check(js, as, f, t)
  
  #----------------------------
  # 2) Dimensions and bookkeeping
  #----------------------------
  
  # na = number of age classes in the final age-based matrix.
  #      For variable durations, this uses the *maximum* possible total juvenile duration
  #      (i.e., the last column of t, usually "max"), then +1 for the adult age class.
  #
  # ns = number of juvenile STAGES (rows of t). (In your “single juvenile stage” case, ns = 1.)
  na = sum(t[[ncol(t)]]) + 1 #summing last column of t gives total juv stage duration (whether it is a matrix or single row matrix or single value).  add one for adults
  ns = nrow(t)
  
  #----------------------------
  # 3) Initialize the projection matrix "m0"
  #----------------------------
  
  # Start with an na x na matrix of zeros
  m0 <- matrix(0, na, na)
  
  # Put fecundity in top row, adult column:
  #  - Adults (last age class) produce newborns (row 1).
  m0[1, na] <- f
  
  # Put adult survival in bottom-right corner:
  #  - Adults can survive and remain adults.
  m0[na, na] <- as
  
  #----------------------------
  # 4) Stage names used in intermediate computations
  #----------------------------
  
  # Stage names: juvenile stages st_1 ... st_ns plus adult stage st_(ns+1)
  stnames = paste0("st_", 1:(ns + 1))
  
  # Name of the adult stage column (last one)
  adstage = stnames[ns + 1]
  
  #----------------------------
  # 5) Build an "age dataframe" used for intermediate calculations
  #----------------------------
  
  # Agedf has one row per age interval (1 to na-1).
  # These are the survival transitions between consecutive ages.
  agedf = data.frame(
    age = 1:(na - 1)
  )
  
  # Add one column per stage to store stage-related probabilities by age.
  # (These are used later to compute weighted-average survival by age.)
  agedf[stnames] <- lapply(1:(ns + 1), function(x) 0)
  
  # Placeholders for:
  #  - agedf$surv : final survival per age interval after weighting by stage occupancy
  #  - agedf$f    : fecundity per age interval (non-zero only when adult)
  agedf$surv = NA
  agedf$f    = NA
  
  # Initialize cohort to begin entirely in the first juvenile stage at the first age interval
  agedf$st_1 = 0
  agedf$st_1[1] = 1
  
  #======================================================================
  # CASE A: VARIABLE STAGE DURATIONS  (t has >1 column; usually dur/min/max)
  #======================================================================
  if(ncol(t) > 1){
    
    #----------------------------------------------------------
    # 6) Build survmat: survival-at-age conditional on stage
    #----------------------------------------------------------
    #
    # survmat will be an (na-1) x (ns+1) matrix:
    #   - rows: ages 1..(na-1)
    #   - cols: each stage st_1..st_(ns+1)
    # Each column says: "If you are in this stage, what is your survival at each age?"
    #
    survmat = NULL
    
    #----------------------------
    # 6A) If js has >1 column: we interpret as "ramp" specification
    #     (typically columns mean + min)
    #----------------------------
    if(ncol(js) > 1){
      
      # ss is stage “targets” for survival at stage ends:
      #  - juvenile stage means (js$mean) plus adult survival appended (as)
      ss = c(js$mean, as)
      
      # startage is intended as the age index where the current juvenile stage begins.
      # NOTE: in this original code, startage is *never updated* inside the loop,
      #       so if ns > 1 this will not shift stages forward correctly.
      #       If ns = 1 (your case), this does not matter.
      startage = 1
      
      # Loop over juvenile stages and create the survival profile for each stage
      for(sss in 1:ns){
        
        # Build a ramp of survival values across the stage.
        # Length = t$dur[sss] + 1
        # Interpretation: “dur” intervals plus an endpoint value (so dur+1 points).
        thisramp = ramp_fun(
          js$mean[sss],     # desired mean survival in this stage
          js$min[sss],      # minimum survival at youngest age in stage
          ss[sss + 1],      # “max”/endpoint survival aligned to next stage target
          t$dur[sss] + 1,    # ramp length (dur+1 points) 
        startbefore = js$startbefore[sss] #KJL do we start the ramp prior to the juvenile age?
          )
        
        # Indices of ages where this stage applies.
        # Length of thisages is also dur+1, matching thisramp.
        thisages = startage:(startage + t$dur[sss])
        
        # Initialize stage-conditional survival vector for all ages (length na-1)
        thissurv = numeric(na - 1)
        
        # Assign ramp values over the ages for this stage
        thissurv[thisages] = thisramp
        
        # For ages BEFORE the stage begins, peg survival at the first ramp value
        thissurv[1:startage] = thisramp[1]
        
        # For ages AFTER the stage ends, peg survival at the last ramp value
        #
        # WARNING: In certain edge cases, (startage + length(thisramp)) can be > (na-1),
        #          and the ":" operator can generate indices that include "na", which
        #          will EXTEND thissurv beyond length (na-1) and cause cbind warnings.
        #
        #kjl this happens when the variable stage duration means that there is just one more age after the mean...
        #kjl do a check to make sure we actually need this!
        if((startage + length(thisramp))<=(na - 1)){ #KJL fix
          thissurv[(startage + length(thisramp)):(na - 1)] = thisramp[length(thisramp)]
        } #KJL
        # Add this stage’s survival profile as a new column of survmat
        survmat = cbind(survmat, thissurv)
        
        # NOTE: If ns > 1, you would normally advance startage here:
        # startage = startage + t$dur[sss]
      }
      
      # Append adult survival column (constant adult survival for all ages)
      survmat = cbind(survmat, rep(as, na - 1))
      
      # Name columns
      colnames(survmat) = stnames
      
    } else {
      #----------------------------
      # 6B) No ramp: juvenile survival is constant within each stage
      #----------------------------
      
      # For each juvenile stage, make a column repeating that stage survival for all ages
      for(sss in 1:ns){
        survmat = cbind(survmat, rep(js[[1]][sss], na - 1))
      }
      
      # Append adult survival column (constant)
      survmat = cbind(survmat, rep(as, na - 1))
      
      # Name columns
      colnames(survmat) = stnames
    }
    
    #----------------------------------------------------------
    # 7) Compute variable “ages at transition” distributions
    #----------------------------------------------------------
    #
    # For each juvenile stage row of t (dur/min/max), build a discrete
    # distribution over possible durations using vari_dur().
    #
    td = lapply(1:nrow(t), function(z){
      vari_dur(t[z,1], t[z,2], t[z,3])
    })
    
    #----------------------------------------------------------
    # 8) Fill in the probability of transitioning to the next stage by age
    #----------------------------------------------------------
    #
    # agedf[[st_(j+1)]] is used as the probability mass of transitioning
    # from stage j to stage j+1 at each age.
    #
    for(j in 1:ns){
      
      thisst = j
      nextst = j + 1
      
      thiscol = sprintf("st_%s", thisst)
      nextcol = sprintf("st_%s", nextst)
      
      thistd = td[[j]]  # data.frame with columns dur and prob
      
      if(j == 1){
        # For the first juvenile stage, transition probability depends directly on age:
        # match durations (thistd$dur) to ages (agedf$age) and assign probabilities.
        agedf[[nextcol]][match(thistd$dur, agedf$age)] = thistd$prob
        
      } else {
        # For later stages, transition timing depends on when individuals entered this stage.
        # Find ages where there is non-negligible probability of being in this stage.
        nonzero_ages = agedf$age[agedf[[thiscol]] > 1e-15]
        
        # For each age of entry into this stage, redistribute probability forward by possible durations
        for(a in 1:length(nonzero_ages)){
          thisa = nonzero_ages[a]
          
          agedf[[nextcol]][thisa + thistd$dur] =
            agedf[[nextcol]][thisa + thistd$dur] +
            agedf[[thiscol]][thisa] * thistd$prob
        }
      }
    } # end loop through stages
    
    # Optional diagnostic: column sums of stage-transition mass
    colSums(agedf[, stnames])
    
    #----------------------------------------------------------
    # 9) Convert transition probabilities into stage-occupancy probabilities (stprob)
    #----------------------------------------------------------
    #
    # stprob[age, stage] = probability an individual is in that stage at that age.
    #
    stprob = agedf[, stnames]
    
    # Everyone starts in first juvenile stage at all ages until they transition
    stprob[,1] = 1
    
    for(s in 1:ns){
      
      # start = current probability mass in stage s (before applying “stay”)
      start = stprob[, s]
      
      # pstay = probability of NOT having transitioned to stage (s+1) yet by that age
      pstay = pmin(1, 1 - cumsum(agedf[, sprintf("st_%s", s + 1)]))
      
      # Avoid tiny negative / floating artifacts
      pstay[pstay < 1e-10] = 0
      
      # Probability of being in stage s at age = probability you were “start” AND stayed
      stprob[, s] = start * pstay
      
      # Probability of being in next stage = probability you were “start” AND did NOT stay
      stprob[, s + 1] = start * (1 - pstay)
    }
    
    # Optional diagnostic: rows should sum to 1 (probability of being in *some* stage)
    rowSums(stprob)
    
    #----------------------------------------------------------
    # 10) Compute final age-specific survival and fecundity
    #----------------------------------------------------------
    
    # Age-specific survival is a weighted mean of stage-conditional survival, weighted by stprob
    agedf$surv = rowSums(survmat * stprob)
    
    # Age-specific fecundity: only adult stage contributes
    agedf$f = stprob[, adstage] * f
    
  } else {
    
    #======================================================================
    # CASE B: FIXED STAGE DURATIONS  (t is effectively 1-column)
    #======================================================================
    
    if(ncol(js) > 1){
      
      #----------------------------------------------------------
      # 11) Fixed durations + ramp case
      #----------------------------------------------------------
      #
      # This part fills agedf$surv from the oldest stage backward to youngest.
      # It uses 'prev' as the stage index (starting at adult = ns+1),
      # and 'ndx' as a pointer to the current end of the agedf$surv vector.
      #
      prev = ns + 1
      ndx  = nrow(agedf)
      
      # Stage survival targets: juvenile means + adult survival
      ss = c(js$mean, as)
      
      # "thismax" is the target (end) survival for the current stage
      thismax = ss[prev]
      
      for(z in 1:ns){
        
        # Duration (years) of the current juvenile stage (working backward)
        thisdur = t[[1]][prev - 1]
        
        # Build a ramp of length thisdur+1 (includes endpoint)
        thisramp = ramp_fun(
          js$mean[prev - 1],
          js$min[prev - 1],
          thismax,
          thisdur + 1
        )
        
        # Fill in the appropriate slice of agedf$surv.
        # The [-length(thisramp)] drops the last element, because that endpoint
        # corresponds to the next (older) stage’s survival target.
        agedf$surv[(ndx - length(thisramp) + 2):ndx] = thisramp[-length(thisramp)]
        
        # Move backward one stage
        prev = prev - 1
        
        # Update "thismax" to the minimum of this ramp, for continuity with the next younger stage
        thismax = min(thisramp)
        
        # Move the index pointer backward by the number of values we just filled
        ndx = ndx - length(thisramp) + 1
        
        # Optional plotting for debugging:
        # plot(agedf[,c("age","surv")])
      }
      
    } else {
      #----------------------------------------------------------
      # 12) Fixed durations + no ramp:
      #     Repeat each stage survival value for its duration
      #----------------------------------------------------------
      agedf$surv = rep(js[[1]], t[[1]])
    }
    
    # For now assume only adults reproduce in fixed-duration formulation
    agedf$f = 0
  }
  
  #======================================================================
  # 13) Convert age-specific survival and fecundity into the final matrix
  #======================================================================
  
  # Put survival values on the sub-diagonal via a diagonal matrix
  m1 <- diag(agedf$surv)
  
  # Insert survival into m0 so that age i transitions to age i+1
  m0[2:na, 1:(na - 1)] = m1
  
  # Put age-specific fecundity in the top row
  m0[1, 1:(na - 1)] = agedf$f
  
  # Return the full age-based projection matrix
  m0
}
# END Functions Script -----------



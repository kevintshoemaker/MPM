
# NOTES ---------

# v6 implements "unroll" with a ramp in juvenile survival up to adulthood
#   to do this, we need separate repro and nonrepro stages for ages affected by variable maturity
#   this will allow individuals to reproduce while juvenile while still varying in survival rates.


# Load packages -----------

library(popbio)  # package for matrix population models
# library(nloptr)  
library(Rsolnp)   # package for nonlinear constrained optimization (eg. for finding maxent probability distributions)
library(memoise)

# Load functions ---------

# s is stage specific survival, l is lambda, T is focal stage duration

## Initial error checking function ------------

init_input_check <- function(js,as,f,t){
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

vari_dur = memoize(vari_dur)

## Survival ramp functions ------------

# sm = 0.75; sa = 0.5; sz = 0.95; y= 5    # note: t should be integer
ramp_fun = function(sm,sa,sz,y){
  opt_func = function(k){
    msurv = mean( sa + (sz-sa) * ( (1-exp(-k * (seq(1,y,length=20)-1) ) ) / (1-exp(-k * (y - 1) ) ) ) )
    abs(msurv-sm)
  }
  opt = optimize(opt_func,c(-10,10))
  thisk = opt$minimum
  sa + (sz-sa) * ( (1-exp(-thisk * ( 1:y - 1 ) ) ) / (1 - exp(-thisk * (y - 1) ) ) )
} 

# plot(1:10,ramp_fun(0.8,0.5,0.96,10),type="l")


# plot(age_at_mat(9,5,12))

## Function for "UNROLL" method  ------------------
  # function for "unrolling" a stage-based matrix with a fixed duration stage into an age-based matrix  [update: allow variable stage durations]
  # currently only implemented for pre-breeding census model. Assumes fertility term already incorporates neonate/first year survival
  # survival terms represent survival rates only while in each stage.

do_unroll <- function(js,as,f,t){
  init_input_check(js,as,f,t)
  na = sum(t[[ncol(t)]])+1; ns = nrow(t)   # na is number of age classes, ns is number of juv stage classes
  stnames = paste0("st_", 1:(ns+1))
  adstage = stnames[ns+1]
  
  if(ncol(t)>1){
    na = sum(t[[3]]) + 1; ns = nrow(t)   # na is number of age classes, ns is number of juv stage classes
    if(ncol(js)>1){
      
    }else{
      survdf=NULL
      for(ss in 1:(ns)){
        survdf = cbind(survdf,rep(js[[1]][ss],na-1))
      }
      survdf = cbind(survdf, rep(as,na-1))
      colnames(survdf)=stnames    # keep as matrix?
      # survdf=as.data.frame(survdf)
    }
      
    m0 <- matrix(0,na,na) ;  m0[1,na] <- f;  m0[na,na] <- as    # construct init matrix
    td = lapply(1:nrow(t),function(z) vari_dur(t[z,1],t[z,2],t[z,3] )   )  # variable ages at transition
    
    agedf = data.frame(    # make data frame to do preliminary computations
      age = 1:(na-1)
    )
    agedf[stnames] <- lapply(1:(ns+1), function(x) 0)
    agedf$surv = NA; agedf$f = NA; agedf$st_1=0; agedf$st_1[1] = 1 
    j=1
    for(j in 1:ns){   # determine probability densities of transitioning at each age given variable durations
      thisst = j ; nextst = j+1
      thiscol = sprintf("st_%s",thisst)
      nextcol = sprintf("st_%s",nextst)
      thistd = td[[j]]
      if(j==1){
        agedf[[nextcol]][match(thistd$dur,agedf$age)] = thistd$prob
      }else{
        nonzero_ages = agedf$age[agedf[[thiscol]]>1e-15 ]
        a=1
        for(a in 1:length(nonzero_ages)){  # loop through all ages of the prior stage and redistribute them 
          thisa = nonzero_ages[a]
          agedf[[nextcol]][thisa+thistd$dur] = agedf[[nextcol]][thisa+thistd$dur] + 
            agedf[[thiscol]][thisa] * thistd$prob
        }
      }
    }# end loop through stages
    
    colSums(agedf[,stnames])   # make sure cols sum to 1...
    
    stprob = agedf[,stnames]
    stprob[,1] = 1  # start with everyone in first stage
    s=2
    for(s in 1:ns){
      start = stprob[,s]   # all who did not transition to the next stage stay
      pstay = pmin(1,1-cumsum(agedf[,sprintf("st_%s",s+1)])) ; pstay[pstay<1e-10] = 0
      stprob[,s] = start * pstay   # need to be in stage currently AND stay in stage
      stprob[,s+1] = start * (1 - pstay) 
    }
    rowSums(stprob)   # make sure rows sum to 1
    
    agedf$surv = rowSums( survdf*stprob )   # survival by age
    agedf$f = stprob[,adstage] * f
    
    # for matrix, transition rates must be conditional on already being in the stage...
      # assume that all individuals must pass through all stages. 
    
    # put matrix together
  
    m1 <- diag(agedf$surv)
    m0[2:na,1:(na-1)] = m1
    m0[1,1:(na-1)] = agedf$f
    
    # agedf[[thiscol]][match(thistd$dur,agedf$age)] = pmax(0,pmin(1,1-cumsum(thistd$prob)) )
    # agedf[[thiscol]][max(thistd$dur):(na-1)] = 0; agedf[[nextcol]][max(thistd$dur):(na-1)] = 1
    
    
    # next: loop through all the stages and mix together all the ages .
    # compute the prob of survival as a weighted average
    # compute fecundity as a weighted average 
    
    for(r in 1:nrow(t)){
      m1 <- diag(ss[r]*(1-pr2))
      m1[na-1,] <- s[r]*pr2
      m0[2:na,1:(na-1)] = m1
    }
  }else{
    if(ncol(js)>1){   # if ramp and no variable aam
      
    }else{    # if no ramp and no variable aam
      m0 <- matrix(0,na,na) ;  m0[1,na] <- f;  m0[na,na] <- as    # construct init matrix
      m1 <- diag(rep(js,t[[1]]))
      m0[2:na,1:(na-1)] = m1
    }
  }
  m0
}

# Generate scenarios -------

# jsurv=0.75; asurv=0.96; fec=0.5; dur=9
scen1 <- gen_scen(jsurv=0.75, asurv=0.96, fec=0.5, dur=9)

scen2 <- gen_scen(jsurv=c(0.75,0.96),asurv=0.96,fec=0.5,dur=c(6,3))

scen3 <- gen_scen(jsurv=c(0.75),asurv=0.96,fec=0.5,dur=data.frame(dur=9,min=6,max=13))   # variable stage duration

scen4 <- gen_scen(jsurv=c(0.6,0.8),asurv=0.96,fec=0.5,dur=data.frame(dur=c(4,3),min=c(2,2),max=c(5,5)))   # variable stage duration

scen5 <- gen_scen(jsurv=data.frame(mean=0.75,min=0.5),asurv=0.95, fec=1.29, dur=data.frame(dur=9,min=6,max=13))

scen6 <- gen_scen(jsurv=data.frame(mean=c(0.6,0.8),min=c(0.3,0.6)),asurv=0.95, fec=1.29, dur=data.frame(dur=c(3,5),min=c(1,3),max=c(3,8) ))

# Run tests --------

## Single juvenile stage  ------------

# js=scen1$jsurv;as=scen1$asurv;f=scen1$fec;t=scen1$dur
mat <- do_aas(js=scen1$jsurv,as=scen1$asurv,f=scen1$fec,t=scen1$dur)    # AAS
mat
lambda(mat)

# js=scen1$jsurv;as=scen1$asurv;f=scen1$fec;t=scen1$dur
mat <- do_unroll(scen1$jsurv,scen1$asurv,scen1$fec,scen1$dur)   # unroll
mat
lambda(mat)       # gives the same result as AAS- as it should!

mat <- do_fas(scen1$jsurv,scen1$asurv,scen1$fec,scen1$dur)  # FAS
mat
lambda(mat)  # 1.065  --- 6.5% growth rate per year   - much higher growth rate!!!  

## Two pre-maturation stages ----------

# js=scen2$jsurv;as=scen2$asurv;f=scen2$fec;t=scen2$dur
mat <- do_aas(scen2$jsurv,scen2$asurv,scen2$fec,scen2$dur)    # AAS
mat
lambda(mat)

mat <- do_unroll(scen2$jsurv,scen2$asurv,scen2$fec,scen2$dur)   # unroll
mat
lambda(mat)       # gives the same result as AAS- as it should!

mat <- do_fas(scen2$jsurv,scen2$asurv,scen2$fec,scen2$dur)  # FAS
mat
lambda(mat)  # still higher, but not as wrong  


## variable juvenile stage duration ---------
js=scen3$jsurv; as=scen3$asurv; f=scen3$fec; t=scen3$dur
mat <- do_aas(scen4$jsurv,scen4$asurv,scen4$fec,scen4$dur)   # aas- gives warning message
mat
lambda(mat) 

js=scen3$jsurv;as=scen3$asurv;f=scen3$fec;t=scen3$dur
mat <- do_unroll(scen3$surv,scen3$fec,scen3$dur)   # very minor difference from fixed duration!!
mat
lambda(mat) 


## multiple variable juvenile stages -----

js=scen4$jsurv;as=scen4$asurv;f=scen4$fec;t=scen4$dur
mat <- do_unroll(scen3$surv,scen3$fec,scen3$dur)   # very minor difference from fixed duration!!
mat
lambda(mat) 

# CHECK AND DEBUG ------------------

## DEBUG STEP: Check out difference between "unroll" method and AAS for 2-year stage...

# m = do_aas(c(0.5,0.95),0.72,2)
# m2 = do_unroll(c(0.5,0.95),0.72,2)
# m;m2
# lambda(m);lambda(m2)   # the unroll method gives much higher - NOW FIXED

## Compare with Kendall's package...
## library(mpmtools)

# install.packages("devtools")    # code to install package from github from lead author of Kendall et al paper. 
# devtools::install_github("BruceKendall/mpmtools")

# check result against kendall's package
# stab = data.frame(stage_name=paste0("s",c(1:length(scen2$surv) )),survival=c(NULL,scen2$surv),maternity=c(NULL,0,0,scen2$fec),duration=c(NULL,scen2$dur,Inf))
# mat2 <- make_stage4age_matrix(stab,approx_method = "AAS",model="post")
# mat2[1,2]  = fec
# mat2          # why does kendall's function increase the adult survival rate for "pre" model?  Also includes a term for adult shrinking in post model. What am I missing?
# mat2[2,2] = surv[2]
# lambda(mat2)   # not totally analogous- not sure what Kendall's function is doing exactly..
# mat2 <- make_stage4age_matrix(stab,approx_method = "unrolled",model="post")
# mat2[1,9]  = fec
# lambda(mat2)     # exactly the same...


# mat2 <- make_stage4age_matrix(stab,approx_method = "FAS",model="pre")
# mat2[1,2]  = fec     # Kendall's function is increasing adult survival to over 100%. Not correct??
# lambda(mat2)

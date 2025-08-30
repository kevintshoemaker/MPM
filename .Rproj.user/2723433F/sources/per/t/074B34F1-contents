

# Load packages -----------

library(popbio)  # package for matrix population models
library(nloptr)  # package for nonlinear constrained optimization (eg. for finding maxent probability distributions)
library(Rsolnp)  

# Load functions ---------

# s is stage specific survival, l is lambda, T is focal stage duration


## Function for AAS method ---------------

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

# n is number of stages (there are n-1 juvenile stages), s is survival of all stages, f is fecundity of final stage, t is duration of juvenile stages
do_aas <- function(s,f,t){
  n=length(s); m0 <- matrix(0,n,n) ; m0[n,n] <- s[n] ;  m0[1,n] <- f    # construct init matrix
  lam = 1  # take an initial guess at lambda
  dif=Inf ; tol=1e-6
  while(dif>tol){
    l = lambda(m0)
    for(g in 1:(n-1)) m0[g:(g+1),g] = thistrans_aas(s[g],l,t[g]) 
    dif = abs(lambda(m0) - l)
  }
  m0
}

## Function for "UNROLL" method  ------------------

# function for "unrolling" a stage-based matrix with a fixed duration stage into an age-based matrix

do_unroll <- function(s,f,t){
  na = sum(t)+1; ns = length(t)   # ny is number of age classes, ns is number of juv stage classes
  m0 <- matrix(0,na,na) ; m0[na,na] <- s[ns+1] ;  m0[1,na] <- f    # construct init matrix
  m1 <- diag(rep(s[1:ns],t))
  m0[2:na,1:(na-1)] = m1
  m0
}

## Function for incorrect "FAS" method ------------

# function for implementing the incorrect "FAS" method (used in Gopher Tortoise SSA model)
do_fas <- function(s,f,t){
  n=length(s); m0 <- matrix(0,n,n) ; m0[n,n] <- s[n] ;  m0[1,n] <- f    # construct init matrix
  gam <- 1/t    # fraction of stage that is in the transition age
  for(g in 1:(n-1)) m0[g:(g+1),g] <- c(s[g]*(1-gam[g]), s[g]*gam[g] )
  m0
}

## Variable age at maturity ------------

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
  sig=(max(ages)-min(ages))/4
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


### Function for allocating age at maturity in variable-maturity age-structured model 
meanam=9;minam=5;maxam=12
age_at_mat <- function(meanam, minam, maxam){
  aam <- minam:maxam
  p0 = dnorm(aam,meanam,(maxam-minam)/4) ; p0= p0/sum(p0)
  bounds = cbind(rep(1e-9,length(aam)),rep(1-1e-9,length(aam)))
  eq = generate_eq(aam,meanam)
  sink("new");
  s = solnp(pars=p0, fun = negent, 
            eqfun = eq, eqB = c(0,0,0),
            ineqfun = ineq, ineqLB = rep(-1e9,length(aam)-2), ineqUB=rep(0,length(aam)-2),
            LB=bounds[,1],UB=bounds[,2])   # 
  sink()
  data.frame(
    ages=aam,
    prob=s$pars
  )
}
# plot(age_at_mat(9,5,12))


# Generate scenarios -------

scen1 <- list()
scen1$surv = c(0.75,0.96) # one survival rate per stage
scen1$fec = 0.5  # annual fecundity (assume pre-breeding census- this represents the number of new one-year-olds entering the population next year)  
scen1$dur = 9  # mean duration of juvenile stage

scen2 <- list()
scen2$surv = c(0.75,0.96,0.96) # one survival rate per stage
scen2$fec = 0.5  # annual fecundity (assume pre-breeding census- this represents the number of new one-year-olds entering the population next year)  
scen2$dur = c(6,3)  # mean duration of juvenile stages

# Run tests --------

## Single juvenile stage  ------------

mat <- do_aas(scen1$surv,scen1$fec,scen1$dur)    # AAS
mat
lambda(mat)

mat <- do_unroll(scen1$surv,scen1$fec,scen1$dur)   # unroll
mat
lambda(mat)       # gives the same result as AAS- as it should!

mat <- do_fas(surv,fec,ma)  # FAS
mat
lambda(mat)  # 1.065  --- 6.5% growth rate per year   - much higher growth rate!!!  

## Two pre-maturation stages ----------

mat <- do_aas(scen2$surv,scen2$fec,scen2$dur)    # AAS
mat
lambda(mat)

mat <- do_unroll(scen2$surv,scen2$fec,scen2$dur)   # unroll
mat
lambda(mat)       # gives the same result as AAS- as it should!

mat <- do_fas(scen2$surv,scen2$fec,scen2$dur)  # FAS
mat
lambda(mat)  # still higher, but not as wrong  


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

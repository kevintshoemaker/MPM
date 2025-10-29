

# Load packages -----------

library(popbio)  # package for matrix population models

# Load functions ---------

# s is stage specific survival, l is lambda, T is stage duration
gamma_aas <- function(s,l,T){
  (s/l)^(T-1) / sum((s/l)^(1:(T-1)))
}

thismat_aas <- function(mat,s,l,T){
  gam = gamma_aas(s,l,T)
  mat[,1] <- c(s*(1-gam),s*gam) 
  mat
}

# function for implementing AAS method 
# assume that there is only one multiyear stage of interest- the first stage (juveniles) 
# s is vector of survival, f is fecundity (scalar), t is duration of juvenile stage
do_aas <- function(n,s,f,t){
  m0 <- matrix(0,n,n) ; m0[n,n] <- s[2] ;  m0[1,n] <- f    # construct init matrix
  lam = 1  # take an initial guess at lambda
  dif=Inf ; tol=1e-6
  while(dif>tol){
    l = lambda(m0)
    m0 = thismat_aas(m0,s[1],l,t) 
    dif = abs(lambda(m0) - l)
  }
  m0
}

# function for "unrolling" a stage-based matrix with a fixed duration stage into an age-based matrix
do_unroll <- function(n,s,f,t){
  m0 <- matrix(0,t,t) ; m0[t,t] <- s[2] ;  m0[1,t] <- f    # construct init matrix
  for(i in 1:(t-1)) m0[i+1,i] <- s[1]
  m0
}

# function for implementing the incorrect "FAS" method (used in Gopher Tortoise SSA model)
do_fas <- function(n,s,f,t){
  m0 <- matrix(0,n,n) ; m0[n,n] <- s[2] ;  m0[1,n] <- f    # construct init matrix
  gam <- 1/(t-1)
  m0[,1] <- c(s[1]*(1-gam), s[1]*gam )
  m0
}

# function for divvying up age at maturity in age-structured model
meanam=9.5;minam=7;maxam=15
age_at_mat <- function(meanam, minam, maxam){
  aam <- minam:maxam
  obj = function(par) {
    c1 = (sum(aam*dnorm(aam,par[1],par[2]))-meanam)^2 
    c2 = ((sum(aam^2*dnorm(aam,par[1],par[2])) - meanam^2) - 0.5*diff(range(aam)))^2
    (c1+c2)
  }
  optpars <- suppressWarnings( optim(obj,par=c(10,2),method="BFGS") )
  list(
    ages=aam,
    dist=dnorm(aam, optpars$par[1],optpars$par[2])
  )
}


# Declare pop params -------

nstag = 2  # number of stages: juv and adult

surv = c(0.75,0.96) # one survival rate per stage

fec = 0.5  # annual fecundity (assume pre-breeding census- this represents the number of new one-year-olds entering the population next year)  

T = 9  # mean age of transition from one stage to the next


# Kendall method "AAS": asymptotic age within stage (correct method 1) ------------

   # formula: [surv/lam)^(T-1)] / Sum({i=0:T-1} (survi / lam)^i)

mat <- do_aas(nstag,surv,fec,T)
lambda(mat)


# Kendall "unroll" method (the other 'correct' method)  --------------

mat <- do_unroll(nstag,surv,fec,T)
lambda(mat)     

# NEXT: write up a "flat age within stage" routine
    #     

mat <- do_fas(nstag,surv,fec,T)
mat

lambda(mat)  # 1.07  --- 7% growth rate per year   - much higher growth rate!!!  


# develop a function that does the AAS method with 2 or more juvenile stages
    # I think all we need to do is implement each stage one at a time from youngest to oldest
    #  or is there a way to do both stages in one optimization step?


# develop a method that allows for variable maturation age 
     # mean stage duration AND a stdev stage duration  
  # probably easiest to start with the unrolled formulation
  # then try to figure out a way to incorporate variability in the AAS formula.. 





# first of all, what does variation in stage duration look like? Should it be normally distributed? 
   # most likely what we might know is that the average individual will mature after 10 years, and the earliest known 
     # maturation occurs at 7 years old, and some individuals may take up to 15 years. Something like that. 
     # not sure if standard deviation makes sense... 
     # for unrolled model, maybe put in mean, min, and max? Then fit a lognormal distribution where the min is the 2.5% quantile from lognormal, 
         # or maybe use a truncated normal distribution with the defined mean, min, and max (this is the maximum entropy solution... )
         # use an iterative procedure to find a discrete distribution close to the specified mean, min and max?





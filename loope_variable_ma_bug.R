
rm(list=ls())

# Load packages -----------

library(popbio)  # package for matrix population models
library(Rsolnp)   # package for nonlinear constrained optimization (eg. for finding maxent probability distributions)
library(memoise)

# Load functions ---------

source("MPM_Functions_v0_1.R")

#scenario: age at maturity is 15 with a range of +/- 2.  
#expectation: reproduction starts at age 13, experiences largest jump from 14 to 15, and is 1.0 (with adult survival rate) at age 17
#incorrect matrix when using "durations" of age-1 for the dur dataframe
scen <- gen_scen(jsurv=c(0.75),asurv=0.96,fec=1,dur=data.frame(dur=14,min=12,max=16))   # variable stage duration
mat <- do_unroll(scen$jsurv,scen$asurv,scen$fec,scen$dur)  
mat
plot(1:ncol(mat),mat[1,],type="s")
abline(v=c(13,15,17),lt="dotted")

#reproduction here starts at 12 and is 1.0 at 16

#correct matrix when using "durations" one longer htan intended
scen <- gen_scen(jsurv=c(0.75),asurv=0.96,fec=1,dur=data.frame(dur=15,min=13,max=17))   # variable stage duration
mat <- do_unroll(scen$jsurv,scen$asurv,scen$fec,scen$dur) 
mat
plot(1:ncol(mat),mat[1,],type="s")
abline(v=c(13,15,17),lt="dotted")

#or am i thinking about the timing of reproduction incorrectly?  I am all tangled up in my brain about this!

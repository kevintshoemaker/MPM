
rm(list=ls())


# Load functions ---------

source("MPM_Functions_v0_2.R")

#scenario: age at maturity is 14 with a range of +/- 2.  
#expectation: reproduction starts at age 12, experiences largest jump from 13 to 14, and is 1.0 (with adult survival rate) at age 16
#incorrect matrix when using "durations" of age-1 for the dur dataframe
scen <- gen_scen(# variable stage duration
  fysurv=0.5,
  jsurv=c(0.75),
  asurv=0.96,
  fec=2,
  dur=data.frame(dur=14,min=12,max=16),
  ramp=F
)   

mat <- do_unroll(scen$fysurv,scen$jsurv,scen$asurv,scen$fec,scen$dur)  
mat
plot(1:ncol(mat),mat[1,],type="s")
abline(v=c(13,15,17),lt="dotted")

#reproduction here starts at 13 and is 1.0 at 17

#  new version seems to work okay??
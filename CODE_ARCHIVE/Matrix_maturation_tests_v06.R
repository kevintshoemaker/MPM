
# NOTES ---------

# v6 implements "unroll" with a ramp in juvenile survival up to adulthood
#   tested for up to 2 stages with variable stage duration
#   still assumes only one reproductive stage (adult) but this could be relaxed in the future


# Load packages -----------

library(popbio)  # package for matrix population models
# library(nloptr)  
library(Rsolnp)   # package for nonlinear constrained optimization (eg. for finding maxent probability distributions)
library(memoise)

# Load functions ---------

# Generate scenarios -------

# jsurv=0.75; asurv=0.96; fec=0.5; dur=9
scen1 <- gen_scen(jsurv=0.75, asurv=0.96, fec=0.5, dur=9)

scen2 <- gen_scen(jsurv=c(0.75,0.96),asurv=0.96,fec=0.5,dur=c(6,3))

scen3 <- gen_scen(jsurv=c(0.75),asurv=0.96,fec=0.5,dur=data.frame(dur=9,min=6,max=13))   # variable stage duration

scen4 <- gen_scen(jsurv=c(0.6,0.8),asurv=0.96,fec=0.5,dur=data.frame(dur=c(4,3),min=c(2,2),max=c(5,5)))   # variable stage duration

scen5 <- gen_scen(jsurv=data.frame(mean=0.75,min=0.5),asurv=0.95, fec=1.29, dur=9)

scen6 <- gen_scen(jsurv=data.frame(mean=c(0.5,0.8),min=c(0.25,0.65)),asurv=0.95, fec=1.29, dur=c(3,4))

scen7 <- gen_scen(jsurv=data.frame(mean=0.75,min=0.5),asurv=0.95, fec=1.29, dur=data.frame(dur=9,min=6,max=13))

scen8 <- gen_scen(jsurv=data.frame(mean=c(0.6,0.8),min=c(0.3,0.6)),asurv=0.95, fec=1.29, dur=data.frame(dur=c(3,5),min=c(1,3),max=c(3,8) ))

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
mat <- do_unroll(scen3$jsurv,scen3$asurv,scen3$fec,scen3$dur)   # very minor difference from fixed duration!!
mat
lambda(mat) 


## multiple variable juvenile stages -----

js=scen4$jsurv;as=scen4$asurv;f=scen4$fec;t=scen4$dur
mat <- do_unroll(scen5$jsurv,scen5$asurv,scen5$fec,scen5$dur)   # very minor difference from fixed duration!!
mat
lambda(mat) 


## ramp with fixed stage duration- one stage -----

js=scen5$jsurv;as=scen5$asurv;f=scen5$fec;t=scen5$dur
mat <- do_unroll(scen2$jsurv,scen2$asurv,scen2$fec,scen2$dur)  
mat
lambda(mat) 


## ramp with fixed stage duration- two stage -----

js=scen6$jsurv;as=scen6$asurv;f=scen6$fec;t=scen6$dur
mat <- do_unroll(scen6$jsurv,scen6$asurv,scen6$fec,scen6$dur)   
mat
lambda(mat) 

## ramp with variable stage duration- one stage -----

js=scen7$jsurv; as=scen7$asurv; f=scen7$fec; t=scen7$dur
mat <- do_unroll(scen6$jsurv,scen6$asurv,scen6$fec,scen6$dur)   
mat
lambda(mat) 

## ramp with variable stage duration- two stage -----

js=scen8$jsurv; as=scen8$asurv; f=scen8$fec; t=scen8$dur
mat <- do_unroll(scen8$jsurv,scen8$asurv,scen8$fec,scen8$dur)   
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

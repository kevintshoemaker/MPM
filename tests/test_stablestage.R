

library(MPM)
library(popbio)
library(ggplot2)


# ?ramp_fun

?gen_scen

basescen = gen_scen(
  0.3,    # fy surv
  0.85,   # juv surv
  0.97,   # ad surv
  1.75,    # fec
  14,     # years in juvenile stage
  T       # ramp in juvenile survival
)

# ?do_unroll
basemat = do_unroll(basescen$fysurv, basescen$jsurv, basescen$asurv, basescen$fec, basescen$dur)
lambda(basemat)

plot(diag(basemat[-1,]),type="l")  # juvenile ramp
mean(diag(basemat[-1,]),type="l")


## examine stable stage distribution for adults

# lambda=0.98;adultsurv=0.97;ma = 15
plot_adult_frac = function(lambda,adultsurv,ma){
  if(lambda<adultsurv) stop("lambda can't be less than adult survival")
  seq = 0:499
  this1 = (adultsurv/lambda)^seq
  adultages = ma+seq
  bins = c(15:199,rep(200,length(seq) - length(15:199) ) )

  adultfracs_num = tapply(this1,bins,sum)
  adultfract_den = sum(adultfracs_num)
  adultfracs = adultfracs_num/adultfract_den
  sum(adultfracs)

  d = data.frame(
    age = as.numeric(names(adultfracs_num)),
    fraction = adultfracs
  )

  print(ggplot(d,aes(age,fraction)) + geom_path(lwd=1.1) + theme_classic())

  return(d)

}

d=plot_adult_frac(lambda=1,adultsurv=0.99,ma=15)

d$fraction[d$age==200]  # fraction of adults greater than 200 years old


sum(d$fraction[d$age>100])


## compute mean survival with senescence

survs = c(rep(0.95,60),seq(0.95,0.1,length=20))
prod(survs)^(1/length(survs))

sum(survs)/(length(survs))



## try running model with senescence.

senescence = gen_scen(
  fysurv = 0.45,
  jsurv  = data.frame(mean = c(0.75), min = c(0.5)),
  asurv  = data.frame(mean = 0.95, old_age = 80, max_age = 100),
  fec    = 1.29,
  dur    = data.frame(dur=9, min=6, max=13),
  ramp   = TRUE
)

nosenescence = gen_scen(
  fysurv = 0.45,
  jsurv  = data.frame(mean = c(0.75), min = c(0.5)),
  asurv  = 0.95,
  fec    = 1.29,
  dur    = data.frame(dur=9, min=6, max=13),
  ramp   = TRUE
)

# ?do_unroll
fys = senescence$fysurv
js = senescence$jsurv
as = senescence$asurv
f = senescence$fec
t = senescence$dur

m=do_unroll(senescence$fysurv,senescence$jsurv,senescence$asurv,senescence$fec,senescence$dur)
popbio::lambda(m)

m2=do_unroll(nosenescence$fysurv,nosenescence$jsurv,nosenescence$asurv,nosenescence$fec,nosenescence$dur)
popbio::lambda(m2)

scen10 <- gen_scen(
  fysurv=0.45,
  jsurv=0.75,
  asurv=data.frame(mean=0.96,old_age=80,max_age=100),
  fec=1.5,
  dur=9,
  ramp=F
)

# fys = scen10$fysurv
# js = scen10$jsurv
# as = scen10$asurv
# f = scen10$fec
# t = scen10$dur
m3=do_unroll(scen10$fysurv,scen10$jsurv,scen10$asurv,scen10$fec,scen10$dur)
popbio::lambda(m3)

m4=do_unroll(scen10$fysurv,scen10$jsurv,data.frame(mean=0.96),scen10$fec,scen10$dur)
popbio::lambda(m4)

scen11 <- gen_scen(
  fysurv=0.45,
  jsurv=data.frame(
    mean=0.75,
    min=0.5
  ),
  asurv=data.frame(mean=0.96,old_age=80,max_age=100),
  fec=2.29,
  dur=data.frame(
    dur=9,
    min=6,
    max=13
  ),
  ramp=T
)

scenx1 <- gen_scen(
  fysurv=0.45,
  jsurv=0.75,
  asurv=data.frame(mean=0.99,old_age=80,max_age=100),
  fec=1.5,
  dur=15,
  ramp=F
)
mx1=do_unroll(scenx1$fysurv,scenx1$jsurv,scenx1$asurv,scenx1$fec,scenx1$dur)
lambda(mx1)

mx2=do_unroll(scenx1$fysurv,scenx1$jsurv,data.frame(mean=0.99),scenx1$fec,scenx1$dur)
lambda(mx2)

scenx2 <- gen_scen(
  fysurv=0.45,
  jsurv=0.75,
  asurv=data.frame(mean=0.96),
  fec=2.0,
  dur=8,
  ramp=F
)
mx2=do_unroll(scenx2$fysurv,scenx2$jsurv,scenx2$asurv,scenx2$fec,scenx2$dur)
lambda(mx2)
lambda(mx1)



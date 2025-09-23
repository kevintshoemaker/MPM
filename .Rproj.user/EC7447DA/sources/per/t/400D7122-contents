

stages = c(1,2)
fecs = seq(0.5,10,length=5)
survs = seq(0.1,0.8,length=5)
Ts = seq(1,10,1)

df1 <- expand.grid(fecs,survs,Ts)
names(df1) = c("fec","surv","T")

df2 <- expand.grid(fecs,survs,survs,Ts,Ts)
names(df2) = c("fec","surv1","surv2","T1","T2")

df1$lam_AAS <- NA
df1$lam_FAS <- NA
df2$lam_AAS <- NA
df2$lam_FAS <- NA
# loop through scenarios and determine lamb

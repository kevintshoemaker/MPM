library(popbio)
library(purrr)
library(dplyr)
library(tidyr)

gamma_aas <- function(s,l,t){         
  sl = s/l
  if(sl==1) 0.5 else (sl^(t-1)*(1-sl)) / (1-sl^t)
}
thistrans_aas <- function(s,l,t){
  gam = gamma_aas(s,l,t)  
  c(s*(1-gam),s*gam)
}

do_aas <- function(s,f,t){
  n=length(s); m0 <- matrix(0,n,n) ; m0[n,n] <- s[n] ;  m0[1,n] <- f
  lam = 1
  dif=Inf ; tol=1e-6
  while(dif>tol){
    l = lambda(m0)
    for(g in 1:(n-1)) m0[g:(g+1),g] = thistrans_aas(s[g],l,t[g]) 
    dif = abs(lambda(m0) - l)
  }
  m0
}

do_fas <- function(s,f,t){
  n=length(s); m0 <- matrix(0,n,n) ; m0[n,n] <- s[n] ;  m0[1,n] <- f
  gam <- 1/t
  for(g in 1:(n-1)) m0[g:(g+1),g] <- c(s[g]*(1-gam[g]), s[g]*gam[g])
  m0
}

fecs  <- seq(0.5, 10,length=25)
jvsurvs <- seq(0.1, 0.8, length=12)
adsurvs <- seq(0.5,0.99,length=12)
TT    <- seq(1,10,1)
jvsurv1 <- seq(0.1, 0.8, length=12)
jvsurv2 <- seq(0.1, 0.8, length=12)

df1 <- expand.grid(fec=fecs, jvsurvs=jvsurvs, adsurvs=adsurvs, TT=TT, stringsAsFactors = FALSE)
df2 <- expand.grid(fec=fecs, jvsurv1=jvsurv1, jvsurv2=jvsurv2, adsurvs=adsurvs, TT=TT, stringsAsFactors = FALSE)

#Vectorized lambda calculations for df1
df1 <- df1 %>%
  mutate(
    lambda_aas = pmap_dbl(list(fec, jvsurvs, adsurvs, TT), function(fec, jvsurvs, adsurvs, TT){
      s <- c(jvsurvs, adsurvs)           # juvenile + adult
      m <- do_aas(s, fec, t=TT)
      lambda(m)
    }),
    lambda_fas = pmap_dbl(list(fec, jvsurvs, adsurvs, TT), function(fec, jvsurvs, adsurvs, TT){
      s <- c(jvsurvs, adsurvs)
      m <- do_fas(s, fec, t=TT)
      lambda(m)
    })
  )

plan(multisession)  # parallel execution

# making functions vectorized
gamma_aas <- function(s, l, t){         
  sl = s/l
  if(abs(sl - 1) < 1e-12) 0.5 else (sl^(t-1)*(1-sl)) / (1 - sl^t)
}

thistrans_aas <- function(s, l, t){
  gam = gamma_aas(s,l,t)  
  c(s*(1-gam), s*gam)
}

# Faster hopefully
do_aas <- function(s, f, t){
  n = length(s)
  m0 <- matrix(0, n, n)
  m0[n,n] <- s[n]
  m0[1,n] <- f
  
  f_obj <- function(l){
    m <- m0
    for(g in 1:(n-1)) m[g:(g+1), g] = thistrans_aas(s[g], l, t[g])
    lambda(m) - l
  }
  
  l_star <- uniroot(f_obj, c(0.01, 10))$root
  for(g in 1:(n-1)) m0[g:(g+1), g] = thistrans_aas(s[g], l_star, t[g])
  m0
}

do_fas <- function(s,f,t){
  n=length(s)
  m0 <- matrix(0,n,n)
  m0[n,n] <- s[n]
  m0[1,n] <- f
  gam <- 1/t
  for(g in 1:(n-1)) m0[g:(g+1),g] <- c(s[g]*(1-gam[g]), s[g]*gam[g])
  m0
}

df1 <- expand.grid(fec=fecs, jvsurvs=jvsurvs, adsurvs=adsurvs, TT=TT, stringsAsFactors = FALSE)

df1 <- df1 %>%
  mutate(
    results = furrr::future_pmap(
      list(fec, jvsurvs, adsurvs, TT),
      function(fec, jvsurvs, adsurvs, TT){
        s <- c(jvsurvs, adsurvs)
        m_aas <- do_aas(s, fec, t=TT)
        m_fas <- do_fas(s, fec, t=TT)
        list(lambda_aas=lambda(m_aas), lambda_fas=lambda(m_fas))
      }
    )
  )
  tidyr::unnest_wider(results)

df1 <- df1 %>%
    mutate(
      lambda_aas = vapply(results, function(x) x$lambda_aas, numeric(1)),
      lambda_fas = vapply(results, function(x) x$lambda_fas, numeric(1))
    ) %>%
    select(-results)
  
#Vectorized lambda calculations for df2
df2 <- df2 %>%
  mutate(
    lambda_aas = pmap_dbl(list(fec, jvsurv1, jvsurv2, TT, adsurvs), function(fec, surv1, surv2, TT, adsurvs){
      s <- c(jvsurv1, jvsurv2, adsurvs)   # two juveniles + adult
      m <- do_aas(s, fec, t=c(TT,TT,1))
      lambda(m)
    }),
    lambda_fas = pmap_dbl(list(fec, surv1, surv2, TT), function(fec, surv1, surv2, TT){
      s <- c(surv1, surv2, adsurvs)
      m <- do_fas(s, fec, t=c(TT,TT,1))
      lambda(m)
    })
  )

#differences
df1 <- df1 %>%
  mutate(lambda_diff = lambda_fas - lambda_aas)
df2 <- df2 %>%
  mutate(lambda_diff = lambda_fas - lambda_aas)

#an attempt at making graphs
library(ggplot2) #recommended on stack exchange for big dataframes

#Subplots for df1 plots
sub1 = df1 %>% filter(fec==fec[3]&adsurvs==adsurvs[8])
ggplot(sub1, aes(x = TT, y = lambda_diff)) +
  geom_point(aes(colour = jvsurvs))+
labs(title = "Difference in λ vs Time to Maturation (T)", x = "TT", y = "λ_diff")

sub2 = df1 %>% filter(TT==2&adsurvs==adsurvs[8])
ggplot(sub2, aes(x = fec, y = lambda_diff)) +
  geom_point(aes(colour = jvsurvs))+
  labs(title = "Difference in λ vs Time to Maturation (T)", x = "Fecundity", y = "λ_diff")

sub3 = df1 %>% filter(fec==fecs[3]&adsurvs==adsurvs[8])
ggplot(sub3, aes(x = jvsurvs, y = lambda_diff)) +
  geom_point(aes(colour = TT))
  labs(title = "Difference in λ vs Time to Maturation (T)", x = "Survival", y = "λ_diff")

sub4 = df1 %>% filter(jvsurvs==jvsurvs[6]&adsurvs==adsurvs[12])
ggplot(sub4, aes(x = TT, y = lambda_diff)) +
  geom_point(aes(colour = fec))
  labs(title = "Difference in λ vs Time to Maturation (T)", x = "Survival", y = "λ_diff")
  
sub5 = df1 %>% filter(TT==2&adsurvs==adsurvs[8])
ggplot(sub5, aes(x = jvsurvs, y = lambda_diff)) +
  geom_point(aes(colour = fec))
  labs(title = "Difference in λ vs Time to Maturation (T)", x = "Survival", y = "λ_diff")
  
sub6 = df1 %>% filter(jvsurvs==0.1&adsurvs==adsurvs[8])
ggplot(sub6, aes(x = fec, y = lambda_diff)) +
  geom_point(aes(colour = TT))
  labs(title = "Difference in λ vs Time to Maturation (T)", x = "Survival", y = "λ_diff")
  
#Heatmaps for df1
subh1 = df1 %>% filter(TT==2&adsurvs==adsurvs[8])
ggplot(subh1, aes(x = fec, y = jvsurvs, fill = lambda_diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    name = expression(lambda[diff])
  ) +
  labs(
    x = "Fecundity",
    y = "Juvenile survival",
    title = expression("Heatmap of " ~ lambda[diff])
  ) 
  theme_minimal(base_size = 14)

subh2 = df1 %>% filter(fec==fecs[3]&adsurvs==adsurvs[8])
ggplot(subh2, aes(x = TT, y = jvsurvs, fill = lambda_diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    name = expression(lambda[diff])
  ) +
  labs(
    x = "Time in Juvenile Stage",
    y = "Juvenile survival",
    title = expression("Heatmap of " ~ lambda[diff])
  ) 
  theme_minimal(base_size = 14)
  
subh3 = df1 %>% filter(jvsurvs==0.1&adsurvs==adsurvs[8])
ggplot(subh3, aes(x = TT, y = fec, fill = lambda_diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    name = expression(lambda[diff])
  ) +
  labs(
    x = "Time in stage",
    y = "Fecundity",
    title = expression("Heatmap of " ~ lambda[diff])
  ) 
  theme_minimal(base_size = 14)
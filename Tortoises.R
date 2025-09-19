library(popbio)
library(purrr)
library(dplyr)
library(tidyr)
library(lhs)

#Tortoise specific parameters
clutchsize <- seq(4,12,1) #multiply by parameters below for fecundity
nestsurv <- seq(0.15,0.55,0.05) #nest success
propfem <- seq(0.5,0.8,0.05) #proportion female
hatchsurv <- seq(0.05,0.4,0.05) #hatchling survival
hatchrate <- seq(0.75,0.95,0.05) #how many eggs become tortoises

jvsurv <- seq(0.63, 0.91,0.01)
adsurv <- seq(0.92,0.99,0.005)
TT <- seq(8,20,1)
varimature <- seq(0,25) #variability in maturation, (0,coefficient of maturation aka fraction of the mean to get to stddev)
#age of maturation and fecundity driven by latitude (so temperature)
#

#latin heighten cube sampling???
?lhs::randomLHS #gets help file for this function
randomLHS(10, 8)

n_params = 8
n_samples = 100
lhs_raw <- randomLHS(n_samples,n_params)
#convert into actual parameter samples
lhs <- lhs_raw

colnames(lhs) <- c("clutchsize","nestsurv", "propfem", "hatchrate", "hatchsurv", "juvsurv", "adsurv", "TT")
lhs <- as.data.frame(lhs)
lhs$clutchsize <-min(clutchsize) + lhs_raw[,1]*diff(range(clutchsize)) #redo this process for the rest of the problems

#lhs$fec <- everything multipled by each other so lhs$clutchsize*lhs$nestsurv*lhs$propfem*lhs$hatchrate*lhs$hatchsurv
#then run scenarios and see what parameters produce growing populations?

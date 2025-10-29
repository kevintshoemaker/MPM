


library(lhs)


?lhs::randomLHS

randomLHS(10,5)


n_params = 8


n_samples = 100


lhs_raw <- randomLHS(n_samples,n_params)


# convert into parameter samples

lhs <- lhs_raw

colnames(lhs) <- c(
  "clutchsize",
  "nestsurv",
  "propfem",
  "hatchrate",
  "hatchsurv",
  "juvsurv",
  "adsurv",
  "TT"
)
lhs <- as.data.frame(lhs)

clutchsize = c(4,12)

lhs$clutchsize <- 4 + lhs_raw[,1] * diff(range(clutchsize))

lhs$fecundity <- lhs$clutchsize*
                 lhs$nestsurv * lhs$propfem *
                 lhs$hatchrate * lhs$hatchsurv



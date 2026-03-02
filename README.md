
# MPM

<!-- badges: start -->
<!-- badges: end -->

MPM provides tools for constructing and comparing age-structured matrix 
population models (MPMs). It implements three methods for handling 
age-structured stage duration — AAS, FAS, and "Unroll" (see Kendall et al. 2019 for details).

MPM supports variable stage durations and within-stage survival ramps.

## Installation

You can install the development version of MPM from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("kevintshoemaker/MPM")
```

Or:

``` r
# install.packages("devtools")
devtools::install_github("kevintshoemaker/MPM")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MPM)
# Single juvenile stage, fixed duration
scen <- gen_scen(fysurv=0.45, jsurv=0.75, asurv=0.96,
                     fec=1.5, dur=9, ramp=FALSE)
    
mat <- do_unroll(scen$fysurv, scen$jsurv, scen$asurv, scen$fec, scen$dur)
popbio::lambda(mat)
```

## Methods

- **Unroll**: recommended method, converts stage parameters into a full 
  age-based Leslie matrix
- **AAS**: Asymptotic Age within Stage (Kendall et al. 2019)
- **FAS**: Flat Age within Stage — included for legacy comparisons only, 
  not recommended

## Citations

Kendall, B.E., Fujiwara, M., Diaz-Lopez, J., Schneider, S., Voigt, J. and Wiesner, S., 2019. 
Persistent problems in the construction of matrix population models. 
Ecological modelling, 406, pp.33-43.  



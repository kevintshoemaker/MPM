# # Code for building R package --------
#
# # Load packages -----------
#
# # install.packages(c("devtools", "roxygen2", "usethis", "testthat"))
#
# library(devtools)
# library(roxygen2)
# library(usethis)
# library(testthat)
#
# # Dependencies ------
#
# usethis::use_package(package = "popbio")
# usethis::use_package(package = "Rsolnp")
# usethis::use_package(package = "memoise")
#
# # Add functions to R directory
#
# ## Initial error checking function ------------
#
# usethis::use_r("init_input_check")
#
# ## Scenario generating function -------------
#
# usethis::use_r("gen_scen")
#
# ## Functions for AAS method ---------------
#
# usethis::use_r("do_aas")
#
#
# ## Function for incorrect "FAS" method ------------
# # function for implementing the incorrect "FAS" method (used in Gopher Tortoise SSA model)
#
# usethis::use_r("do_fas")
#
#
# ## Variable age at maturity functions ------------
#
#
# ### Function for allocating stage durations in variable stage duration model
#
# usethis::use_r("vari_dur")
#
# ## Survival ramp functions ------------
#
# usethis::use_r("ramp_fun")
#
#
# ## Function for "UNROLL" method  ------------------
#
# usethis::use_r("do_unroll")
#
#
## Unit tests etc--------------

# usethis::use_testthat()
# usethis::use_test("do_unroll")

## Generate documentation ------

# devtools::document()

# # Check for problems (run this often!) --------
# devtools::check()
#
# # Install locally to test  ----------
# devtools::install()
#
# # Load without installing (fastest during development)  --------
# devtools::load_all()


#
# # END Functions Script -----------

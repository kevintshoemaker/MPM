test_that("unroll matches AAS for single fixed-duration stage", {
  scen <- gen_scen(fysurv=0.45, jsurv=0.75, asurv=0.96,
                   fec=1.5, dur=9, ramp=FALSE)
  mat_aas    <- do_aas(scen$fysurv, scen$jsurv, scen$asurv, scen$fec, scen$dur)
  mat_unroll <- do_unroll(scen$fysurv, scen$jsurv, scen$asurv, scen$fec, scen$dur)
  expect_equal(popbio::lambda(mat_aas), popbio::lambda(mat_unroll), tolerance=1e-5)
})

test_that("FAS gives higher lambda than AAS in typical application", {
  scen <- gen_scen(fysurv=0.45, jsurv=0.75, asurv=0.96,
                   fec=1.5, dur=9, ramp=FALSE)
  lam_aas <- popbio::lambda(do_aas(scen$fysurv, scen$jsurv, scen$asurv, scen$fec, scen$dur))
  lam_fas <- popbio::lambda(do_fas(scen$fysurv, scen$jsurv, scen$asurv, scen$fec, scen$dur))
  expect_gt(lam_fas, lam_aas)
})

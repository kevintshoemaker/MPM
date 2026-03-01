# mpmR 0.2.0

* Reimagined `do_unroll()` to allow survival ramp to start at the previous 
  stage's survival rate OR at a user-specified minimum value.
* Added `ramp_fun()` as an exported function with full documentation.

# mpmR 0.1.0

* Initial implementation of `do_unroll()` with ramp using maximum entropy 
  distribution, monotonically increasing or decreasing between fixed 
  min and max survival with fixed mean.
* Implemented `do_aas()`, `do_fas()`, and supporting functions.

# Notes
This package currently assumes only one reproductive stage (adult) but this could be relaxed in the future.
For now, we assume a pre-breeding census model. Therefore "fecundity" terms must always include survival of newborns to age 1.
I may implement post-breeding census model but I'm reluctant because it will inevitably add confusion!
However, I think if we can design a way to input parameters that is clear in either model, then we can potentially include an argument for constructing pre- vs post-breeding census matrices



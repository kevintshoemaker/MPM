#' MPM: Matrix Population Models for Stage-Structured Populations in R
#'
#' Provides tools for constructing and comparing age- and stage-based matrix
#' population models (MPMs). Implements three methods for handling stage
#' duration, when stages are defined in terms of age:
#'
#' - **AAS** (Asymptotic Age-within Stage): iterative method described by Kendall et al. 2019.
#' - **FAS** (Fixed Age within Stage): simple but biased method, included for
#'   legacy comparisons.
#' - **Unroll**: converts age-based stages into a full age-based Leslie
#'   matrix, supporting variable stage durations and within-stage survival
#'   ramps.
#'
#'#' @references
#' Kendall, B.E., Fujiwara, M., Diaz-Lopez, J., Schneider, S., Voigt, J. and Wiesner, S., 2019.
#' Persistent problems in the construction of matrix population models.
#'
#' All methods assume a **pre-breeding census** model.
#'
#' @keywords internal
#' @importFrom stats optimize
#' @importFrom utils tail
#' @keywords internal
"_PACKAGE"


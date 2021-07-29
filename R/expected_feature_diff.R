#' Expected number of not optimally balanced features as defined by a threshold
#'
#' @description Calculate expected number of features showing balance difference greater than a threshold
#' 
#' @md
#' @param n.feature (`numeric`) total number of features
#' @param n.arm1 (`numeric`) number of patients in comparison arm.
#' @param n.arm2  (`numeric`) number of patients in control arm
#' @param threshold (`numeric`) positive number(s) for threshold to compare to.
#' @note The output number indicates when running a randomized
#' trial with n.arm1 and n.arm2 samples in two arms and n.feature features
#' are of interest, the expected number of
#' features showing balance difference greater than threshold.
#' p = Prob(|Y|>threshold) is calculated from t distribution.
#' With n.feature features in total,
#' expected number of features with abs value > threshold can be calculated from Binomial(n.feature, p)
#' 
#' @return return a `numeric` vector for expected number of unbalanced features
#' @examples
#' expected_feature_diff(n.feature = 10, n.arm1 = 240, n.arm2 = 300, threshold = 0.2)
#' expected_feature_diff(n.feature = 10, n.arm1 = 240, n.arm2 = 300, threshold = c(0.1, 0.25))
#' @importFrom stats pt
#' @export

expected_feature_diff <- function(n.feature, n.arm1, n.arm2, threshold) {
  stopifnot(all(threshold > 0))
  threshold.t <- threshold / sqrt((1 / n.arm1) + 1 / (n.arm2))
  df <- n.arm1 + n.arm2 - 2
  stopifnot(threshold.t >= 0)
  p <- (1 - pt(threshold.t, df)) * 2
  res <- n.feature * p
  names(res) <- paste("Expected # features >", threshold)
  res
}

#' A simulated meta-analytic data set
#'
#' Data simuated from a true model with a three-way interaction between three
#' moderators: m1, m2 and m3. If the values of the three moderators are all "B"s
#' the ture effect size will be 0.80. Otherwise, the true effect size is 0.
#'
#' @docType data
#'
#' @usage data(SimData)
#'
#' @keywords datasets
#' @format A data frame of 120 studies with 5 moderators
#' \itemize{
#'   \item efk: The effect size of each study
#'   \item vark: The sampling variance of each study
#'   \item m1 to m5: Five randomly generated moderators. m1 and m2 have two levels (A and B),
#'   whereas m3, m4 and m5 have three levels (A, B and C)
#'   }
"SimData"

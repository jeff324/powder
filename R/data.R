#' Simulated example dataset from the LBA (varies drift rate across conditions)
#'
#' A simulated dataset from the LBA with 10 subjects, 2 conditions, and 300 trials per condition.
#' Simulates correct mean drift rate varying across conditions.
#' The parameters used to generate the null dataset were the following:
#' * starting point variability, A: 1
#' * relative threshold, b: 0.4
#' * correct drift rate mean condition 1, vc.1: 2
#' * correct drift rate mean condition 2, vc.2: 3.5
#' * error drift rate mean, ve: 1
#' * drift rate variability, svc, sve: 1
#' * non-decision time, t0: 0.3
#' @md
#' @format A list with 10 elements, each representing a simulated subject. Each element contains the following:
#' \describe{
#'   \item{Cond}{condition number}
#'   \item{Time}{response time}
#'   \item{Correct}{whether the response was error (2) or correct (1)}
#' }
"drift"

#' Simulated example dataset from the LBA (no parameter varied across conditions)
#'
#' A simulated dataset from the LBA with 10 subjects, 2 conditions, and 300 trials per condition.
#' Simulates no parameters varying across conditions.
#' The parameters used to generate the null dataset were the following:
#' * starting point variability, A: 1
#' * relative threshold, b: 0.4
#' * correct drift rate mean vc: 3
#' * error drift rate mean, ve: 1
#' * drift rate variability, svc, sve: 1
#' * non-decision time, t0: 0.3
#' @md
#' @format A list with 10 elements, each representing a simulated subject. Each element contains the following:
#' \describe{
#'   \item{Cond}{condition number}
#'   \item{Time}{response time}
#'   \item{Correct}{whether the response was error (2) or correct (1)}
#' }
"null"

#' The dataset from Rae et al. (2014)
#'
#' @format A list with 34 elements, each representing a subject. Each element contains the following:
#' \describe{
#'   \item{Cond}{condition number}
#'   \item{Correct}{whether the response was error (2) or correct (1)}
#'   \item{Time}{response time}
#' }
"rae"

#' The "Simple" dataset from Evans and Brown (2017)
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{Cond}{condition number}
#'   \item{Correct}{whether the response was error (2) or correct (1)}
#'   \item{Time}{response time}
#' }
"individual"

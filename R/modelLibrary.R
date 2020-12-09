#' Theophylline dataset, as included in the NONMEM distribution.
#'
#' This is comprised of 12 subjects with each 12 concentrations.
#'
#' @format A data frame with 12x13 rows: 1 administration and 12 observation rows
#' \describe{
#'   \item{ID}{Subject ID}
#'   \item{AMT}{Theophylline dose, in mg}
#'   \item{TIME}{Time of observation, in days after first dose}
#'   \item{CONC}{Concentration, in mg/L}
#'   \item{EVID}{0 for observation, 1 for administration}
#'   \item{WT}{Weight in kg}
#' }
"theopp"


#' NlmixrUI object describing theophylline 1-compartment model with oral absorption
#'
#' The model was estimated by nlmixr SAEM
#'
#' @format A fitted Nlmixr model
"theopp_nlmixr"


#' Phenobarbitol dataset, as included in the NONMEM distribution.
#'
#' This is comprised of 59 subjects with sparse sampled concentrations.
#' There are 155 measured concentrations and 589 administrations.
#'
#' @format A data frame with 744 rows
#' \describe{
#'   \item{ID}{Subject ID}
#'   \item{TIME}{Time of observation, in days after first dose}
#'   \item{AMT}{phenobarbitol dose, in mg}
#'   \item{WT}{Weight in kg}
#'   \item{APGR}{APGAR score, a measure for sedation}
#'   \item{DV}{Concentration, in mg/L}
#'   \item{EVID}{1 for administration, 0 for observation}
#'   \item{MDV}{1 for administration, 0 for observation}
#' }
"pheno"


#' NlmixrUI object describing Phenobarbitol 1-compartment model with IV administration
#' and allometric scaling on weight
#'
#' The model was estimated by nlmixr SAEM
#' @format A fitted Nlmixr model
"pheno_nlmixr"

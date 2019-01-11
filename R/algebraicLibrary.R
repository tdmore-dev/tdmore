## Library of functions
#' A tdmore model to predict concentration in a 1cpt PK model with bolus administration
#'
#' @param TVV Typical value for volume
#' @param TVCL Typical value for clearance
#' @param OV Omega for volume
#' @param OCL Omega for clearance
#' @param res_var Residual error model
#'
#' @return a tdmore model
#' @export
#'
#' @examples
#' m1 <- pk1cptivbolusVCL(
#'   TVV=10, TVCL=5,
#'   OV=0.20, OCL=0.30,
#'   errorModel(prop=0.1)
#' )
#' summary(m1)
pk1cptivbolusVCL <- function(TVV=10, TVCL=5, OV=0.20, OCL=0.30, res_var) {
  m1 <- algebraic({
    function(t, TIME, AMT, EV, ECL) {
      V = TVV * exp(EV)
      CL = TVCL * exp(ECL)
      k = CL / V
      tD = TIME
      D = AMT
      ifelse(t >= tD, D / V * exp(-k*(t-tD)), 0)
    }
  })
  tdmore(m1,
          res_var=res_var,
          omega=c(EV=OV, ECL=OCL))
}


#' A tdmore model to predict concentration in a 1cpt PK model with bolus administration
#'
#' @param TVKA Typical value for absorption constant
#' @param TVV Typical value for volume
#' @param TVCL Typical value for clearance
#' @param OKA Typical value for clearance
#' @param OV Omega for volume
#' @param OCL Omega for clearance
#' @param res_var Residual error model
#'
#' @return a tdmore model
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' m1 <- pk1cptoralbolusVCL(
#'   TVKA=0.5, TVV=10, TVCL=5,
#'   OKA=0.70, OV=0.20, OCL=0.30,
#'   errorModel(prop=0.1)
#' )
#' summary(m1)
pk1cptoralbolusVCL <- function(TVKA=0.5, TVV=10, TVCL=5,
                               OKA=0.70, OV=0.20, OCL=0.30,
                               res_var) {
  m1 <- algebraic({
    function(t, TIME, AMT, EKA, EV, ECL) {
      V = TVV * exp(EV)
      CL = TVCL * exp(ECL)
      Ka = TVKA * exp(EKA)
      k = CL / V
      tD = TIME
      D = AMT
      ifelse(t >= tD,
             D/V * Ka/(Ka-k) * (exp(-k*(t-tD))-exp(-Ka*(t-tD)))
             ,
             0)
    }
  })
  tdmore(m1,
         res_var=res_var,
         omega=c(EV=OV, ECL=OCL, EKA=OKA))
}

#' This base function can be used to develop an algebraic model
#' for a 2cpt pk model with oral absorption.
#'
#' The caller should have access to the following values:
#' t, TIME, AMT
#' as well as model parameters
#' CL, V1, Q, V2
#'
#' The following alternative parametrization is also permitted
#' K12, K21, K
#'
#' Optionally, TLAG can be defined.
#'
#' @examples
#' m1 <- algebraic(function(t, TIME, AMT, EV, ECL) {
#'     KA <- 0.3
#'     CL <- 3 * exp(ECL)
#'     V1 <- 10 * exp(EV)
#'     Q <- 3
#'     V2 <- 800
#'     pk2cptoralbolus()
#' })
#' @export
pk2cptivbolus <- function() {
  env <- parent.frame()
  with(env, {
    if(!exists('K12')) K12 <- Q/V1
    if(!exists('K21')) K21 <- Q/V2
    if(!exists('K')) K <- CL/V1
    if(!exists('TLAG')) TLAG <- 0

    Beta = 1/2*(K12+K21+K-sqrt((K12+K21+K)^2 - 4*K21*K) )
    Alpha=K21*K/Beta
    A = (1/V1) * (Alpha - K21) / (Alpha - Beta)
    B = (1/V1) * (Beta - K21) / (Beta - Alpha)

    D=AMT
    tD=TIME
    browser()
    D*ifelse(t-tD <= TLAG,
             0,
             A*exp(-Alpha*(t-tD-TLAG))+B*exp(-Beta*(t-tD-TLAG))
    )
  })
}

#' This base function can be used to develop an algebraic model
#' for a 2cpt pk model with oral absorption.
#'
#' The caller should have access to the following values:
#' t, TIME, AMT
#' as well as model parameters
#' KA, CL, V1, Q, V2
#'
#' The following alternative parametrization is also permitted
#' KA, K12, K21, K
#'
#' Optionally, TLAG can be defined.
#'
#' @examples
#' m1 <- algebraic(function(t, TIME, AMT, EV, ECL) {
#'     KA <- 0.3
#'     CL <- 3 * exp(ECL)
#'     V1 <- 10 * exp(EV)
#'     Q <- 3
#'     V2 <- 800
#'     pk2cptoralbolus()
#' })
#' @export
pk2cptoralbolus <- function() {
  env <- parent.frame()
  with(env, {
    if(!exists('K12')) K12 <- Q/V1
    if(!exists('K21')) K21 <- Q/V2
    if(!exists('K')) K <- CL/V1
    if(!exists('TLAG')) TLAG <- 0

    Beta = 1/2*(K12+K21+K-sqrt((K12+K21+K)^2 - 4*K21*K) )
    Alpha=K21*K/Beta
    A = (KA/V1) * (K21-Alpha) / ((KA-Alpha)*(Beta-Alpha))
    B = (KA/V1) * (K21-Beta)/((KA-Beta)*(Alpha-Beta))

    D=AMT
    tD=TIME
    D*ifelse(t-tD <= TLAG,
             0,
             A*exp(-Alpha*(t-tD-TLAG))+B*exp(-Beta*(t-tD-TLAG))
             -(A+B)*exp(-KA*(t-tD-TLAG))
    )
  })
}

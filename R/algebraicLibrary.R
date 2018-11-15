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


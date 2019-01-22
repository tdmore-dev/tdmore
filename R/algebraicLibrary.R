## Library of functions


## These equations come from PFIM / INSERM / Monolix PKPD library
# Core algebraic equations ----------------------------------------------------

#'
#' Algebraic equations for PK models
#'
#' @param t sampling times
#' @param tD dosing time
#' @param D dosing amount
#' @param V volume of distribution
#' @param K elimination constant
#'
#' @return vector of same size as `t`, with the concentrations
#' @name pk_
NULL

prepare1cpt_ <- function(env=parent.frame()) {
  assign('tD', env$TIME, envir=env)
  assign('D', env$AMT, envir=env)
  assign('Tau', env$II, envir=env)
  assign('TInf', env$AMT / env$RATE, envir=env)
  if( length(env$TInf) > 0 && env$TInf > env$Tau ) warning("Infusion time larger than interdose interval")
}
#' @rdname pk_
#' @export
pk1cptiv_ <- function(t, TIME, AMT, V, K, SS=0, II=Inf) {
  prepare1cpt_()
  if(SS==0)
    D / V * exp(-K*(t-tD))
  else
    D / V * exp(-K*(t-tD)) / (1 - exp(-K*Tau) )
}

#' @rdname pk_
#' @export
pk1cptinfusion_ <- function(t, TIME, AMT, RATE, V, K, SS=0, II=Inf) {
  prepare1cpt_()
  if(SS==0)
    D/TInf * 1/(K*V)* ifelse(t-tD <= TInf,
           (1 - exp(-K*(t-tD)) ),
           (1 - exp(-K*TInf) )*exp(-K*(t-tD-TInf))
    )
  else
    D/TInf * 1/(K*V) * ifelse(t-tD <= TInf,
                              (1-exp(-K*(t-tD))) +   exp(-K*Tau)*(1-exp(-K*TInf))*exp(-K*(t-tD-TInf))/(1-exp(-K*Tau)),
                              (1-exp(-K*TInf)) * exp(-K*(t-tD-TInf)) / (1-exp(-K*Tau))
    )
}
#' @rdname pk_
#' @export
pk1cptoral_ <- function(t, TIME, AMT, V, K, KA, SS=0, II=Inf) {
  prepare1cpt_()
  if(SS==0)
    D/V * KA/(KA-K) * ( exp(-K*(t-tD)) - exp(-KA*(t-tD))    )
  else
    D/V * KA/(KA-K) * ( exp(-K*(t-tD))/(1-exp(-K*Tau)) - exp(-KA*(t-tD))/(1-exp(-KA*Tau))    )
}
## Convenience function to assign helper variables
prepare2cpt_ <- function(env=parent.frame()) {
  assign('tD', env$TIME, envir=env)
  assign('D', env$AMT, envir=env)
  assign('Tau', env$II, envir=env)
  assign('TInf', env$AMT / env$RATE, envir=env)
  if( length(env$TInf) > 0 && env$TInf > env$Tau ) warning("Infusion time larger than interdose interval")

  Beta =with(env, 1/2*(K12 + K21 + K - sqrt((K12+K21+K)^2-4*K21*K) ) )
  assign('Beta', Beta, envir=env)

  Alpha=with(env, K21*K/Beta )
  assign('Alpha', Alpha, envir=env)

  A=with(env, 1/V*(Alpha-K21)/(Alpha-Beta) )
  assign('A', A, envir=env)

  B=with(env, 1/V*(Beta-K21)/(Beta-Alpha) )
  assign('B', B, envir=env)
}
#' @rdname pk_
#' @export
pk2cptiv_ <- function(t, TIME, AMT, V, K, K12, K21, SS=0, II=Inf) {
  prepare2cpt_()
  if(SS==0)
    D*(A*exp(-Alpha*(t-tD)) + B*exp(-Beta*(t-tD)) )
  else
    D*(A*exp(-Alpha*(t-tD)) / (1-exp(-Alpha*Tau)) + B*exp(-Beta*(t-tD)) / (1-exp(-Beta*Tau)) )
}
#' @rdname pk_
#' @export
pk2cptinfusion_ <- function(t, TIME, AMT, RATE, V, K, K12, K21, SS=0, II=Inf) {
  prepare2cpt_()
  if(SS==0)
    D/TInf * ifelse(t-tD <= TInf,
                    A/Alpha*(1-exp(-Alpha*(t-tD))) +
                      B/Beta*(1-exp(-Beta*(t-tD)))
                             ,
                    A/Alpha*(1-exp(-Alpha*TInf)) * exp(-Alpha*(t-tD-TInf)) +
                      B/Beta*(1-exp(-Beta*TInf)) * exp(-Beta*(t-tD-TInf))
    )
  else
    D/TInf * ifelse(t-tD <= TInf,
                    A/Alpha*(
                      (1-exp(-Alpha*(t-tD))) +
                        exp(-Alpha*Tau) * (1-exp(-Alpha*TInf))* exp(-Alpha*(t-tD-TInf)) / (1-exp(-Alpha*Tau))
                    ) +
                      B/Beta*(
                        (1-exp(-Beta*(t-tD))) +
                          exp(-Beta*Tau) * (1-exp(-Beta*TInf))* exp(-Beta*(t-tD-TInf)) / (1-exp(-Beta*Tau))
                      )
                    ,
                    A/Alpha*(1-exp(-Alpha*TInf)) * exp(-Alpha*(t-tD-TInf)) / (1-exp(-Alpha*Tau)) +
                      B/Beta*(1-exp(-Beta*TInf)) * exp(-Beta*(t-tD-TInf)) / (1-exp(-Beta*Tau))
    )
}
#' @rdname pk_
#' @export
pk2cptoral_ <- function(t, TIME, AMT, V, K, KA, K12, K21, SS=0, II=Inf) {
  prepare2cpt_()
  A <- KA / (KA-Alpha) * A
  B <- KA / (KA-Beta) * B
  if(SS==0)
    D * ( A*exp(-Alpha*(t-tD))+ B*exp(-Beta*(t-tD)) - (A+B)*exp(-KA*(t-tD)) )
  else
    D * (A*exp(-Alpha*(t-tD))/(1-exp(-Alpha*Tau)) +
           B*exp(-Beta*(t-tD))/(1-exp(-Beta*Tau)) -
           (A+B)*exp(-KA*(t-tD))/(1-exp(-KA*Tau))
    )
}
## Convenience function to assign helper variables
prepare3cpt_ <- function(env=parent.frame()) {
  assign('tD', env$TIME, envir=env)
  assign('D', env$AMT, envir=env)
  assign('Tau', env$II, envir=env)
  assign('TInf', env$AMT / env$RATE, envir=env)
  if( length(env$TInf) > 0 && env$TInf > env$Tau ) warning("Infusion time larger than interdose interval")

  result <- with(env, {
    a0 = K*K21*K31
    a1 = K*K31+K21*K31+K21*K13+K*K21+K31*K12
    a2 = K+K12+K13+K21+K31

    p=a1-(a2)^2/3
    q=2*a2^3/27-a1*a2/3+a0
    r1=sqrt(-(p^3/27))
    r2=2*r1^(1/3)
    phi=acos(-q/(2*r1))/3
    Alpha=-(cos(phi)*r2-a2/3)
    Beta=-(cos(phi+2*pi/3)*r2-a2/3)
    Gamma=-(cos(phi+4*pi/3)*r2-a2/3)
    list(Alpha=Alpha, Beta=Beta, Gamma=Gamma)
  })

  assign('Alpha', result$Alpha, envir=env)
  assign('Beta', result$Beta, envir=env)
  assign('Gamma', result$Gamma, envir=env)

  A=with(env, 1/V*(K21-Alpha)/(Alpha-Beta)*(K31-Alpha)/(Alpha-Gamma) )
  assign('A', A, envir=env)

  B=with(env, 1/V*(K21-Beta)/(Beta-Alpha)*(K31-Beta)/(Beta-Gamma) )
  assign('B', B, envir=env)

  C=with(env, 1/V*(K21-Gamma)/(Gamma-Beta)*(K31-Gamma)/(Gamma-Alpha) )
  assign('C', C, envir=env)
}

#' @rdname pk_
#' @export
pk3cptiv_ <- function(t, TIME, AMT, V, K, K12, K21, K13, K31, SS=0, II=Inf) {
  prepare3cpt_()
  if(SS==0)
    D*(A*exp(-Alpha*(t-tD))+B*exp(-Beta*(t-tD))+C*exp(-Gamma*(t-tD)))
  else
    D*(
      A*exp(-Alpha*(t-tD))/(1-exp(-Alpha*Tau))+
        B*exp(-Beta*(t-tD))/(1-exp(-Beta*Tau))+
        C*exp(-Gamma*(t-tD))/(1-exp(-Gamma*Tau))
    )
}
#' @rdname pk_
#' @export
pk3cptinfusion_ <- function(t, TIME, AMT, RATE, V, K, K12, K21, K13, K31, SS=0, II=Inf) {
  prepare3cpt_()
  if(SS==0)
    D/TInf * ifelse(t-tD <= TInf,
                         A/Alpha*(1-exp(-Alpha*(t-tD))) +
                           B/Beta *(1-exp(-Beta *(t-tD))) +
                           C/Gamma*(1-exp(-Gamma*(t-tD)))  ,
                         A/Alpha*(1-exp(-Alpha*TInf))  *exp(-Alpha*(t-tD-TInf)) +
                           B/Beta *(1-exp(-Beta *TInf))*exp(-Beta*(t-tD-TInf)) +
                           C/Gamma*(1-exp(-Gamma*TInf))*exp(-Gamma*(t-tD-TInf))
    )
  else
    D/TInf * ifelse(t-tD <= TInf,
                    A/Alpha*(
                      (1-exp(-Alpha*(t-tD))) +
                        exp(-Alpha*Tau) * (1-exp(-Alpha*TInf))* exp(-Alpha*(t-tD-TInf)) / (1-exp(-Alpha*Tau))
                    ) +
                      B/Beta*(
                        (1-exp(-Beta*(t-tD))) +
                          exp(-Beta*Tau) * (1-exp(-Beta*TInf))* exp(-Beta*(t-tD-TInf)) / (1-exp(-Beta*Tau))
                      ) +
                      C/Gamma*(
                        (1-exp(-Gamma*(t-tD))) +
                          exp(-Gamma*Tau) * (1-exp(-Gamma*TInf))* exp(-Gamma*(t-tD-TInf)) / (1-exp(-Gamma*Tau))
                      )
                    ,
                    A/Alpha*(1-exp(-Alpha*TInf)) * exp(-Alpha*(t-tD-TInf)) / (1-exp(-Alpha*Tau)) +
                      B/Beta*(1-exp(-Beta*TInf)) * exp(-Beta*(t-tD-TInf)) / (1-exp(-Beta*Tau)) +
                      C/Gamma*(1-exp(-Gamma*TInf)) * exp(-Gamma*(t-tD-TInf)) / (1-exp(-Gamma*Tau))
    )
}
#' @rdname pk_
#' @export
pk3cptoral_ <- function(t, TIME, AMT, V, K, KA, K12, K21, K13, K31, SS=0, II=Inf) {
  prepare3cpt_()

  A <- KA / (KA-Alpha) * A
  B <- KA / (KA-Beta) * B
  C <- KA / (KA-Gamma) * C
  if(SS==0)
    D*(
       A*exp(-Alpha*(t-tD))+
       B*exp(-Beta*(t-tD)) +
       C*exp(-Gamma*(t-tD)) -
       (A+B+C)*exp(-KA*(t-tD))
    )
  else
    D*(
      A*exp(-Alpha*(t-tD))/(1-exp(-Alpha*Tau)) +
        B*exp(-Beta*(t-tD))/(1-exp(-Beta*Tau)) +
        C*exp(-Gamma*(t-tD))/(1-exp(-Gamma*Tau)) -
        (A+B+C)*exp(-KA*(t-tD))/(1-exp(-KA*Tau))
    )
}


# Getting parameters from parent frame ------------------------------------
pkFromParentFrame <- function(fName=as.character(sys.call())[1], env=parent.frame()) {
  # Infer name of function call
  fAlgebraicName <- paste0(fName, "_")

  # Find that function, and its arguments
  myFun <- match.fun(FUN=fAlgebraicName)
  myFormals <- formals(myFun)

  # Check for the arguments
  args <- mget(names(myFormals), envir=env, ifnotfound=myFormals)
  missing <- unlist( lapply(args, is.name) )
  if(any(missing)) {
    stop("Could not execute algebraic function ", fAlgebraicName,", ",
         "argument(s) ", paste(names(myFormals)[missing], collapse=", "), " not found."  )
  }

  do.call(myFun, args)
}

#' @export
pk1cptiv <- pkFromParentFrame
#' @export
pk1cptinfusion <- pkFromParentFrame
#' @export
pk1cptoral <- pkFromParentFrame
#' @export
pk2cptiv <- pkFromParentFrame
#' @export
pk2cptinfusion <- pkFromParentFrame
#' @export
pk2cptoral <- pkFromParentFrame
#' @export
pk3cptiv <- pkFromParentFrame
#' @export
pk3cptinfusion <- pkFromParentFrame
#' @export
pk3cptoral <- pkFromParentFrame


# NONMEM-style definition -------------------------------------------------
#
# TODO: Let's keep this as a TODO for now
# NONMEM has plenty of extra parameters that need to be supported.
# As an example, the simple ADVAN1 has the following parameters:
#S1 - Scale for central compartment (also called SC)
#S2 - Scale for output compartment (also called S0)
#F1 - Bioavailability for central compartment
#R1 - Rate for central compartment
#D1 - Duration for central compartment
#ALAG1 - Absorption lag for central compartment
#F0 - Output fraction (also called F2, FO)
#XSCALE - Time scale
#MTIME(i) - Model event times
#
# Due to the amount of work required, we will keep this on the back-burner.
#
# ADVAN1 One Compartment Linear Model
# ADVAN2 One Compartment Linear Model with First Order Absorption
# ADVAN3 Two Compartment Linear Model
# ADVAN4 Two Compartment Linear Model with First Order Absorption
# ADVAN5 General Linear Model
# ADVAN6 General Nonlinear Model
# ADVAN7 General Linear Model with Real Eigenvalues
# ADVAN8 General Nonlinear Model with Stiff Differential Equations
# ADVAN9 General Nonlinear Model with Equilibrium Compartments
# ADVAN10 One Compartment Model with Michaelis-Menten Elimination
# ADVAN11 Three Compartment Linear Model
# ADVAN12 Three Compartment Linear Model with First Order Absorption
# ADVAN13 General Nonlinear Model (sometimes better than ADVAN6) (NM7)

# TRANS1 Dummy, or null, translator; (K,V)
# TRANS2 Used with ADVAN1 and ADVAN2 (CL,V).
# TRANS3 Used with ADVAN3 and ADVAN4 (CL,V,Q,VSS).
# TRANS4 Used with ADVAN3, ADVAN4(CL,V1,Q,V2)
# TRANS4 Used with ADVAN11, ADVAN12 (CL,V1,Q2,V2,Q3,V3)
# TRANS5 Used with ADVAN3 and ADVAN4 (AOB, ALPHA, BETA)
# TRANS6 Used with ADVAN3, ADVAN4, ADVAN11, ADVAN12 (ALPHA, BETA, K21)

# This function captures all required parameters for the PK
# function evaluation from eval.env, and assigns them to assign.env
# applyTrans <- function(advan, trans, eval.env, assign.env) {
#   if(advan==1) {
#     if(trans==1) {
#       assign('K', eval.env$K, assign.env)
#       assign('V', eval.env$S1, assign.env)
#     } else if (trans ==2) {
#
#     } else {
#       stop("TRANS", trans," cannot be used with ADVAN", advan)
#     }
#   }
#
#   if(trans == 1) {
#     ## Nothing to do, just check
#     K=with(eval.env, {K})
#     V=with(eval.env, {V})
#   }
# }
# ADVAN1 <- function(env=parent.frame(), trans=1) {
#   applyTrans(1, trans, eval.env=env, assign.env=sys.frame())
#   if(!is.null(env$RATE))
#     pk1cptiv()
#   else
#     pk1cptivbolus()
# }
# ADVAN2 <- function(env=parent.frame(), trans=1) {
#   applyTrans(2, trans, eval.env=env, assign.env=sys.frame())
#   pk1cptoral()
# }
# ADVAN3 <- function(env=parent.frame()) {
#   applyTrans(3, trans, eval.env=env, assign.env=sys.frame())
#   if(!is.null(env$RATE))
#     pk2cptiv()
#   else
#     pk2cptivbolus()
# }
# #' Two Compartment Linear Model with First Order Absorption
# ADVAN4 <- function(env=parent.frame()) {
#   applyTrans(4, trans, eval.env=env, assign.env=sys.frame())
#   pk2cptoral()
# }
# #' General Linear Model
# ADVAN5 <- function(env=parent.frame()) {
#   stop("Not supported")
# }
# #' General Nonlinear Model
# ADVAN6 <- function(env=parent.frame()) {
#   stop("Not supported")
# }
# #' General Linear Model with Real Eigenvalues
# ADVAN7 <- function(env=parent.frame()) {
#   stop("Not supported")
# }
# #' General Nonlinear Model with Stiff Differential Equations
# ADVAN8 <- function(env=parent.frame()) {
#   stop("Not supported")
# }
# #' General Nonlinear Model with Equilibrium Compartments
# ADVAN9 <- function(env=parent.frame()) {
#   stop("Not supported")
# }
# #' One Compartment Model with Michaelis-Menten Elimination
# ADVAN10 <- function(env=parent.frame()) {
#   stop("Not supported")
# }
# #' Three Compartment Linear Model
# ADVAN11 <- function(env=parent.frame()) {
#   applyTrans(11, trans, eval.env=env, assign.env=sys.frame())
#   if(!is.null(env$RATE))
#     pk3cptiv()
#   else
#     pk3cptivbolus()
# }
# #' Three Compartment Linear Model with First Order Absorption
# ADVAN12 <- function(env=parent.frame()) {
#   applyTrans(12, trans, eval.env=env, assign.env=sys.frame())
#   pk3cptoral()
# }
# #' General Nonlinear Model (sometimes better than ADVAN6) (NM7)
# ADVAN13 <- function(env=parent.frame()) {
#   stop("Not supported")
# }

autoPick <- function(values, env, params, transform) {
  ## Search what is the appropriate function
  match <- FALSE
  for(i in seq_along(params)) {
    if(all(params[[i]] %in% names(values) )) {
      match <- TRUE
      break
    }
  }

  if(!match)
    stop("No matching model could be found.
         Specify the model using any of the following variables:",
         as.character(params))

  ## Is there some kind of transformation defined for that function?
  transfun <- transform[[i]]
  if(!is.null(transfun)) {
    # do this again with the adapted parameters
    values <- transfun(values)
    return( autoPick(values, env, params, transform) )
  }

  fName <- names(params[i])

  ## Translate names and assign to environment
  values <- values[ unlist( params[[i]] ) ] #take the required parameters for the model

  tNames <- names( params[[i]] ) #names to translate
  keepNames <- nchar(tNames)==0
  tNames[keepNames] <- names(values)[keepNames]

  for(j in seq_along(values)) assign(tNames[j], values[[j]], envir=env) #and assign to environment

  ## And delegate to function
  do.call(fName, list(env=env))
}

# nlmixr-style definition -------------------------------------------------
### from https://github.com/nlmixrdevelopment/nlmixr/blob/21f1165370358dbd31038ed1fd9406d48898fe46/R/saem_fit.R#L599
nlmixrParams = list(
  pk3cptoral=list(K="KE", "V", "K12", "K21", "K13", "K31", "KA"),
  pk3cptinfusion=list(K="KE", "V", "K12", "K21", "K13", "K31", "RATE"),
  pk3cptiv=list(K="KE", "V", "K12", "K21", "K13", "K31"),
  pk2cptoral=list(K="KE", "V", "K12", "K21", "KA"),
  pk2cptinfusion=list(K="KE", "V", "K12", "K21", "RATE"),
  pk2cptiv=list(K="KE", "V", "K12", "K21"),
  pk1cptoral=list(K="KE", "V", "KA"),
  pk1cptinfusion=list(K="KE", "V", "RATE"),
  pk1cptiv=list(K="KE", "V"),

  pk3cptoral=list("CL", "V", "CLD", "VT", "CLD2", "VT2", "KA"),
  pk3cptinfusion=list("CL", "V", "CLD", "VT", "CLD2", "VT2", "RATE"),
  pk3cptiv=list("CL", "V", "CLD", "VT", "CLD2", "VT2"),
  pk2cptoral=list("CL", "V", "CLD", "VT", "KA"),
  pk2cptinfusion=list("CL", "V", "CLD", "VT", "RATE"),
  pk2cptiv=list("CL", "V", "CLD", "VT"),
  pk1cptoral=list("CL", "V", "KA"),
  pk1cptinfusion=list("CL", "V", "RATE"),
  pk1cptiv=list("CL", "V")
)

#' This function guesses which model should be calculated, based on the available variables
#' in the caller environment
#'
#' Please note this function also modifies the caller environment with the
#' appropriate values for the call to the algebraic function.
#' As an example, it will define 'K'.
#'
#' @export
linCmt <- function(env=parent.frame()) {
  values <- as.list(env) ## TODO: does not include parent values...
  transform <- list()
  transform[10:18] <- list(function(values){
    values$KE = values$CL / values$V
    values$K12 = values$CLD / values$V
    values$K21 = values$CLD / values$VT
    values$K13 = values$CLD2 / values$V
    values$K31 = values$CLD2 / values$VT2
    values
  })
  autoPick(values, env, nlmixrParams, transform)
}

# Monolix-style definition ------------------------------------------------
monolixParams = list(
  pk3cptoral=list(K="k", "V", K12="k12", K21="k21", K13="k13", K31="k31", KA="ka"),
  pk3cptinfusion=list(K="k", "V", K12="k12", K21="k21", K13="k13", K31="k31", "RATE"),
  pk3cptiv=list(K="k", "V", K12="k12", K21="k21", K13="k13", "k31"),
  pk2cptoral=list(K="k", "V", K12="k12", K21="k21", KA="ka"),
  pk2cptinfusion=list(K="k", "V", K12="k12", K21="k21", "RATE"),
  pk2cptiv=list(K="k", "V", K12="k12", K21="k21"),
  pk1cptoral=list(K="k", "V", KA="ka"),
  pk1cptinfusion=list(K="k", "V", "RATE"),
  pk1cptiv=list(K="k", "V"),

  pk3cptoral=list("Cl", "V", "Q2", "V2", "Q3", "V3", "ka"),
  pk3cptinfusion=list("Cl", "V", "Q2", "V2", "Q3", "V3", "RATE"),
  pk3cptiv=list("Cl", "V", "Q2", "V2", "Q3", "V3"),
  pk2cptoral=list("Cl", "V", "Q2", "V2", "ka"),
  pk2cptinfusion=list("Cl", "V", "Q2", "V2", "RATE"),
  pk2cptiv=list("Cl", "V", "Q2", "V2"),
  pk1cptoral=list("Cl", "V", "ka"),
  pk1cptinfusion=list("Cl", "V", "RATE"),
  pk1cptiv=list("Cl", "V")
)
## See http://monolix.lixoft.com/data-and-models/pkmodel1/
#' @export
pkmodel <- function(...) {
  env=parent.frame()
  qArgsNamed = rlang::quos(..., .named=TRUE)
  values = list(...)
  names(values) = names(qArgsNamed)

  transform <- list()
  transform[10:18] <- list(function(values){
    #clearance-volume
    values$k = values$Cl / values$V
    values$k12 = values$Q2 / values$V
    values$k21 = values$Q2 / values$V2
    values$k13 = values$Q3 / values$V
    values$k31 = values$Q3 / values$V3
    values
  })

  autoPick(values, env, monolixParams, transform)
}

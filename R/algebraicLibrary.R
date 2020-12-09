## Library of algebraic functions
### See https://www.sciencedirect.com/science/article/pii/S1056871915000362

### Ask R CMD CHECK to ignore the use of locally assigned variables
if(getRversion() >= "2.15.1")  utils::globalVariables(c("D", "tD", "TInf", "Tau", "Alpha", "Beta", "Gamma", "A", "B", "C", "lambda1", "lambda2", "lambda3", "E1", "E2", "E3"))

#'
#' Algebraic equations for PK models
#'
#' @param t sampling times
#' @param TIME dosing time
#' @param AMT dosing amount
#' @param V volume of distribution
#' @param K elimination constant
#' @param SS 0 for no steady-state, 1 for steady-state
#' @param II inter-dose interval, required for steady-state calculations
#' @param RATE infusion rate
#' @param KA absorption constant
#' @param K12 transfer rate from central compartment to peripheral compartment A
#' @param K21 transfer rate from peripheral compartment A to central compartment
#' @param K13 transfer rate from central compartment to peripheral compartment B
#' @param K31 transfer rate from peripheral compartment B to central compartment
#' @param A0 initial state for A0 at TIME
#' @param A1 initial state for A1 at TIME
#' @param A2 initial state for A2 at TIME
#' @param A3 initial state for A3 at TIME
#'
#' @return vector of same size as `t`, with the concentrations
#' @name pk_
NULL

prepare1cpt_ <- function(env=parent.frame()) {
  assign('tD', env$TIME, envir=env)
  assign('D', env$AMT, envir=env)
  assign('Tau', env$II, envir=env)
  if(!is.null(env$RATE)) {
    TInf <- env$AMT / env$RATE
    if(env$RATE == 0) TInf <- 0
    assign('TInf', TInf, envir=env)
  }
  if( length(env$TInf) > 0 && env$TInf > env$Tau ) warning("Infusion time larger than interdose interval")
}

#' @rdname pk_
#' 
pk1cptiv_ <- function(t, A1=0, TIME, AMT, V, K, SS=0, II=Inf) {
  prepare1cpt_()
  if(SS==0) {
    A1 = A1 + D
    A1 = A1 * exp(-K*(t-tD))
  } else {
    A1 = D * exp(-K*(t-tD)) / (1 - exp(-K*Tau) ) #ignore initial state
  }
  list(A1=A1, CONC=A1/V)
}

#' @rdname pk_
#' 
pk1cptinfusion_ <- function(t, A1=0, TIME, AMT, RATE, V, K, SS=0, II=Inf) {
  prepare1cpt_()
  if(SS==0) {
    A1i = A1
    A1 = RATE * 1/(K)* ifelse(t-tD <= TInf,
           (1 - exp(-K*(t-tD)) ),
           (1 - exp(-K*TInf) )*exp(-K*(t-tD-TInf))
    )
    A1 = A1 + A1i * exp(-K*(t-tD) ) # add initial
  } else {
    A1 = RATE * 1/(K) * ifelse(t-tD <= TInf,
                              (1-exp(-K*(t-tD))) +   exp(-K*Tau)*(1-exp(-K*TInf))*exp(-K*(t-tD-TInf))/(1-exp(-K*Tau)),
                              (1-exp(-K*TInf)) * exp(-K*(t-tD-TInf)) / (1-exp(-K*Tau))
    )
  }
  list(A1=A1, CONC=A1/V)
}

#' @rdname pk_
#' 
pk1cptoral_ <- function(t, A0=0, A1=0, TIME, AMT, V, K, KA, SS=0, II=Inf) {
  prepare1cpt_()
  if(SS==0) {
    A0i = A0 + D
    A1i = A1
    A1 = A0i*KA/(KA-K)*( exp(-K*(t-tD)) - exp(-KA*(t-tD))    ) +
       A1i*exp(-K*(t-tD))
    A0 = A0i * exp(-KA*(t-tD))
  } else {
    A0 = D * exp(-KA*(t-tD)) / (1 - exp(-KA*Tau) ) #same formula as for 1cpt IV bolus SS
    A1 = D * KA/(KA-K) * ( exp(-K*(t-tD))/(1-exp(-K*Tau)) - exp(-KA*(t-tD))/(1-exp(-KA*Tau))    )
  }
  list(A0=A0, A1=A1, CONC=A1/V)
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


  k10 <- env$K
  k12 <- env$K12
  k21 <- env$K21
  k20 <- 0

  #calculate hybrid rate constants
  lambda1 <- 0.5*(k12+k21+k10+sqrt((k12+k21+k10)^2-4*k21*k10))
  lambda2 <- 0.5*(k12+k21+k10-sqrt((k12+k21+k10)^2-4*k21*k10))
  assign('lambda1', lambda1, envir=env)
  assign('lambda2', lambda2, envir=env)
}

#' @rdname pk_
#' 
pk2cptiv_ <- function(t, A1=0, A2=0, TIME, AMT, V, K, K12, K21, SS=0, II=Inf) {
  prepare2cpt_()
  if(SS==0) {
    A1i = A1 + D
    A2i = A2
    E1 <- K+K12
    E2 <- K21
    t <- t-tD
    A1 = (((A1i*E2+A2i*K21)-A1i*lambda1)*exp(-t*lambda1)-((A1i*E2+A2i*K21)-A1i*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
    A2 = (((A2i*E1+A1i*K12)-A2i*lambda1)*exp(-t*lambda1)-((A2i*E1+A1i*K12)-A2i*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
    list(A1=A1, A2=A2, CONC=A1/V)
  } else {
    CONC = D*(A*exp(-Alpha*(t-tD)) / (1-exp(-Alpha*Tau)) + B*exp(-Beta*(t-tD)) / (1-exp(-Beta*Tau)) )
    list(CONC=CONC)
  }
}
#' @rdname pk_
#' 
pk2cptinfusion_ <- function(t, A1=0, A2=0, TIME, AMT, RATE, V, K, K12, K21, SS=0, II=Inf) {
  prepare2cpt_()
  if(SS==0) {
    k10 <- K
    k12 <- K12
    k21 <- K21
    k20 <- 0
    E1 <- k10+k12
    E2 <- k21+k20

    #calculate hybrid rate constants
    lambda1 = 0.5*((E1+E2)+sqrt((E1+E2)^2-4*(E1*E2-k12*k21)))
    lambda2 = 0.5*((E1+E2)-sqrt((E1+E2)^2-4*(E1*E2-k12*k21)))

    ## Calculate until end of infusion
    InfEnd <- tD + AMT / RATE
    DuringInfusion <- t < InfEnd
    if(all(DuringInfusion)) stop("Infusion will not be done when calculating next state! ERROR")
    tDuring <- c( t[DuringInfusion], InfEnd )
    tAfter <- c(InfEnd, t[ !DuringInfusion ])
    t <- tDuring - tD

    A1last <- A1
    A2last <- A2
    Doserate <- RATE

    A1term1 = (((A1last*E2+Doserate+A2last*k21)-A1last*lambda1)*exp(-t*lambda1)-((A1last*E2+Doserate+A2last*k21)-A1last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
    A1term2 = Doserate*E2*(1/(lambda1*lambda2)+exp(-t*lambda1)/(lambda1*(lambda1-lambda2))-exp(-t*lambda2)/(lambda2*(lambda1-lambda2)))
    A1During <- A1term1+A1term2    #Amount in the central compartment

    A2term1 = (((A2last*E1+A1last*k12)-A2last*lambda1)*exp(-t*lambda1)-((A2last*E1+A1last*k12)-A2last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
    A2term2 = Doserate*k12*(1/(lambda1*lambda2)+exp(-t*lambda1)/(lambda1*(lambda1-lambda2))-exp(-t*lambda2)/(lambda2*(lambda1-lambda2)))
    A2During <- A2term1+A2term2   #Amount in the peripheral compartment

    outAfter <- pk2cptiv_(t=tAfter, A1=A1During[length(A1During)], A2=A2During[length(A2During)], TIME=InfEnd, AMT=0, V=V, K=K, K12=K12, K21=K21)
    A1 = c(A1During[-length(A1During)], outAfter$A1[-1])
    A2 = c(A2During[-length(A2During)], outAfter$A2[-1])

    list(A1=A1, A2=A2, CONC=A1/V)
  } else {
    CONC = RATE * ifelse(t-tD <= TInf,
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
    list(CONC=CONC)
  }
}
#' @rdname pk_
#' 
pk2cptoral_ <- function(t, A0=0, A1=0, A2=0, TIME, AMT, V, K, KA, K12, K21, SS=0, II=Inf) {
  prepare2cpt_()
  A <- KA / (KA-Alpha) * A
  B <- KA / (KA-Beta) * B
  if(SS==0) {
    k20 <- K
    k23 <- K12
    k32 <- K21
    KA  <- KA
    k30 <- 0
    E2 <- k20+k23
    E3 <- k32+k30

    #calculate hybrid rate constants
    lambda1 = 0.5*((E2+E3)+sqrt((E2+E3)^2-4*(E2*E3-k23*k32)))
    lambda2 = 0.5*((E2+E3)-sqrt((E2+E3)^2-4*(E2*E3-k23*k32)))

    t <- t-tD
    A2last <- A1
    A3last <- A2
    A1last <- A0 + AMT

    A2term1 = (((A2last*E3+A3last*k32)-A2last*lambda1)*exp(-t*lambda1)-((A2last*E3+A3last*k32)-A2last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
    A2term2 = A1last*KA*(exp(-t*KA)*(E3-KA)/((lambda1-KA)*(lambda2-KA))+exp(-t*lambda1)*(E3-lambda1)/((lambda2-lambda1)*(KA-lambda1))+exp(-t*lambda2)*(E3-lambda2)/((lambda1-lambda2)*(KA-lambda2)))
    A1 = A2term1+A2term2  #Amount in the central compartment

    A3term1 = (((A3last*E2+A2last*k23)-A3last*lambda1)*exp(-t*lambda1)-((A3last*E2+A2last*k23)-A3last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
    A3term2 = A1last*KA*k23*(exp(-t*KA)/((lambda1-KA)*(lambda2-KA))+exp(-t*lambda1)/((lambda2-lambda1)*(KA-lambda1))+exp(-t*lambda2)/((lambda1-lambda2)*(KA-lambda2)))
    A2 = A3term1+A3term2  #Amount in the peripheral compartment

    A0 = A1last*exp(-t*KA)

    list(A0=A0, A1=A1, A2=A2, CONC=A1/V)
  } else {
    CONC = D * (A*exp(-Alpha*(t-tD))/(1-exp(-Alpha*Tau)) +
           B*exp(-Beta*(t-tD))/(1-exp(-Beta*Tau)) -
           (A+B)*exp(-KA*(t-tD))/(1-exp(-KA*Tau))
    )
    list(CONC=CONC)
  }
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


  k10 <- env$K
  k12 <- env$K12
  k21 <- env$K21
  k13 <- env$K13
  k31 <- env$K31
  k20 <- 0
  k30 <- 0
  E1 <- k10+k12+k13
  E2 <- k21+k20
  E3 <- k31+k30
  assign('E1', E1, envir=env)
  assign('E2', E2, envir=env)
  assign('E3', E3, envir=env)

  #calculate hybrid rate constants
  a <- E1+E2+E3
  b <- E1*E2+E3*(E1+E2)-k12*k21-k13*k31
  c <- E1*E2*E3-E3*k12*k21-E2*k13*k31

  m <- (3*b - a^2)/3
  n <- (2*a^3 - 9*a*b + 27*c)/27
  Q <- (n^2)/4 + (m^3)/27

  alpha <- sqrt(-1*Q)
  beta <- -1*n/2
  gamma <- sqrt(beta^2+alpha^2)
  theta <- atan2(alpha,beta)

  lambda1 <- a/3 + gamma^(1/3)*(cos(theta/3) + sqrt(3)*sin(theta/3))
  lambda2 <- a/3 + gamma^(1/3)*(cos(theta/3) - sqrt(3)*sin(theta/3))
  lambda3 <- a/3 -(2*gamma^(1/3)*cos(theta/3))
  assign('lambda1', lambda1, envir=env)
  assign('lambda2', lambda2, envir=env)
  assign('lambda3', lambda3, envir=env)
}

#' @rdname pk_
#' 
pk3cptiv_ <- function(t, A1=0, A2=0, A3=0, TIME, AMT, V, K, K12, K21, K13, K31, SS=0, II=Inf) {
  prepare3cpt_()
  if(SS==0) {
    t <- t - tD
    k21 = K21
    k12 = K12
    k13 = K13
    k31 = K31
    A1last <- A1 + AMT
    A2last <- A2
    A3last <- A3

    B = A2last*k21+A3last*k31
    C = E3*A2last*k21+E2*A3last*k31
    I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31
    J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21

    A1term1 = A1last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A1term2 = exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))

    A1 <- (A1term1+A1term2)  #Amount in the central compartment

    A2term1 = A2last*(exp(-t*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A2term2 = exp(-t*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))

    A2 <- A2term1+A2term2             #Amount in the first-peripheral compartment

    A3term1 = A3last*(exp(-t*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A3term2 = exp(-t*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))

    A3 <- A3term1+A3term2            #Amount in the second-peripheral compartment

    list(A1=A1, A2=A2, A3=A3, CONC=A1/V)
  } else {
    CONC <- D*(
      A*exp(-Alpha*(t-tD))/(1-exp(-Alpha*Tau))+
        B*exp(-Beta*(t-tD))/(1-exp(-Beta*Tau))+
        C*exp(-Gamma*(t-tD))/(1-exp(-Gamma*Tau))
    )
    list(CONC=CONC)
  }
}
#' @rdname pk_
#' 
pk3cptinfusion_ <- function(t, A1=0, A2=0, A3=0, TIME, AMT, RATE, V, K, K12, K21, K13, K31, SS=0, II=Inf) {
  prepare3cpt_()
  if(SS==0) {
    k10 <- K
    k12 <- K12
    k21 <- K21
    k13 <- K13
    k31 <- K31
    k20 <- 0
    k30 <- 0
    E1 <- k10+k12+k13
    E2 <- k21+k20
    E3 <- k31+k30

    #calculate hybrid rate constants
    a <- E1+E2+E3
    b <- E1*E2+E3*(E1+E2)-k12*k21-k13*k31
    c <- E1*E2*E3-E3*k12*k21-E2*k13*k31

    m <- (3*b - a^2)/3
    n <- (2*a^3 - 9*a*b + 27*c)/27
    Q <- (n^2)/4 + (m^3)/27

    alpha <- sqrt(-1*Q)
    beta <- -1*n/2
    gamma <- sqrt(beta^2+alpha^2)
    theta <- atan2(alpha,beta)

    lambda1 <- a/3 + gamma^(1/3)*(cos(theta/3) + sqrt(3)*sin(theta/3))
    lambda2 <- a/3 + gamma^(1/3)*(cos(theta/3) - sqrt(3)*sin(theta/3))
    lambda3 <- a/3 -(2*gamma^(1/3)*cos(theta/3))

    A1last <- A1
    A2last <- A2
    A3last <- A3
    Doserate <- RATE

    ## Calculate until end of infusion
    InfEnd <- tD + AMT / RATE
    DuringInfusion <- t < InfEnd
    if(all(DuringInfusion)) stop("Infusion will not be done when calculating next state! ERROR")
    tDuring <- c( t[DuringInfusion], InfEnd )
    tAfter <- c(InfEnd, t[ !DuringInfusion ])
    t <- tDuring - tD

    B = A2last*k21+A3last*k31
    C = E3*A2last*k21+E2*A3last*k31
    I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31
    J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21

    A1term1 = A1last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A1term2 = exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))
    A1term3 = Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))

    A1During <- A1term1+A1term2+A1term3    #Amount in the central compartment

    A2term1 = A2last*(exp(-t*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A2term2 = exp(-t*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))
    A2term3 = Doserate*k12*(E3/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))

    A2During <- A2term1+A2term2+A2term3    #Amount in the first-peripheral compartment

    A3term1 = A3last*(exp(-t*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A3term2 = exp(-t*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))
    A3term3 = Doserate*k13*(E2/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))

    A3During <- A3term1+A3term2+A3term3  #Amount in the second-peripheral compartment


    outAfter <- pk3cptiv_(t=tAfter,
                          A1=A1During[length(A1During)],
                          A2=A2During[length(A2During)],
                          A3=A3During[length(A3During)],
                          TIME=InfEnd,
                          AMT=0,
                          V=V, K=K, K12=K12, K21=K21, K13=K13, K31=K31)
    A1 = c(A1During[-length(A1During)], outAfter$A1[-1])
    A2 = c(A2During[-length(A2During)], outAfter$A2[-1])
    A3 = c(A3During[-length(A3During)], outAfter$A3[-1])

    list(A1=A1, A2=A2, A3=A3, CONC=A1/V)
  } else {
    CONC = RATE * ifelse(t-tD <= TInf,
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
    list(CONC=CONC)
  }
}
#' @rdname pk_
#' 
pk3cptoral_ <- function(t, A0=0, A1=0, A2=0, A3=0, TIME, AMT, V, K, KA, K12, K21, K13, K31, SS=0, II=Inf) {
  prepare3cpt_()

  A <- KA / (KA-Alpha) * A
  B <- KA / (KA-Beta) * B
  C <- KA / (KA-Gamma) * C
  if(SS==0) {
    k20 <- K
    k23 <- K12
    k32 <- K21
    k24 <- K13
    k42 <- K31
    KA  <- KA
    k30 <- 0
    k40 <- 0
    E2 <- k20+k23+k24
    E3 <- k32+k30
    E4 <- k42+k40

    #calculate hybrid rate constants
    a <- E2+E3+E4
    b <- E2*E3+E4*(E2+E3)-k23*k32-k24*k42
    c <- E2*E3*E4-E4*k23*k32-E3*k24*k42

    m <- (3*b - a^2)/3
    n <- (2*a^3 - 9*a*b + 27*c)/27
    Q <- (n^2)/4 + (m^3)/27

    alpha <- sqrt(-1*Q)
    beta <- -1*n/2
    gamma <- sqrt(beta^2+alpha^2)
    theta <- atan2(alpha,beta)

    lambda1 <- a/3 + gamma^(1/3)*(cos(theta/3) + sqrt(3)*sin(theta/3))
    lambda2 <- a/3 + gamma^(1/3)*(cos(theta/3) - sqrt(3)*sin(theta/3))
    lambda3 <- a/3 -(2*gamma^(1/3)*cos(theta/3))

    t <- t-tD
    A1last <- A0 + AMT
    A2last <- A1
    A3last <- A2
    A4last <- A3

    B = A3last*k32+A4last*k42
    C = E4*A3last*k32+E3*A4last*k42
    I = A2last*k23*E4-A3last*k24*k42+A4last*k23*k42
    J = A2last*k24*E3+A3last*k24*k32-A4last*k23*k32

    A2term1 = A2last*(exp(-t*lambda1)*(E3-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E3-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E3-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A2term2 = exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))
    A2term3 = A1last*KA*(exp(-t*lambda1)*(E3-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(KA-lambda1))+exp(-t*lambda2)*(E3-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(KA-lambda2))+exp(-t*lambda3)*(E3-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(KA-lambda3))+exp(-t*KA)*(E3-KA)*(E4-KA)/((lambda1-KA)*(lambda2-KA)*(lambda3-KA)))
    A1 = A2term1+A2term2+A2term3   #Amount in the central compartment

    A3term1 = A3last*(exp(-t*lambda1)*(E2-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A3term2 = exp(-t*lambda1)*(I-A2last*k23*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A2last*k23*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A2last*k23*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))
    A3term3 = A1last*KA*k23*(exp(-t*lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(KA-lambda1))+exp(-t*lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(KA-lambda2))+exp(-t*lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(KA-lambda3))+exp(-t*KA)*(E4-KA)/((lambda1-KA)*(lambda2-KA)*(lambda3-KA)))
    A2 = A3term1+A3term2+A3term3  #Amount in the first-peripheral compartment

    A4term1 = A4last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))
    A4term2 = exp(-t*lambda1)*(J-A2last*k24*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A2last*k24*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A2last*k24*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))
    A4term3 = A1last*KA*k24*(exp(-t*lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(KA-lambda1))+exp(-t*lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(KA-lambda2))+exp(-t*lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(KA-lambda3))+exp(-t*KA)*(E3-KA)/((lambda1-KA)*(lambda2-KA)*(lambda3-KA)))
    A3 = A4term1+A4term2+A4term3  #Amount in the second-peripheral compartment

    A0 = A1last*exp(-t*KA)
    list(A0=A0, A1=A1, A2=A2, A3=A3, CONC=A1/V)
  } else {
    CONC = D*(
      A*exp(-Alpha*(t-tD))/(1-exp(-Alpha*Tau)) +
        B*exp(-Beta*(t-tD))/(1-exp(-Beta*Tau)) +
        C*exp(-Gamma*(t-tD))/(1-exp(-Gamma*Tau)) -
        (A+B+C)*exp(-KA*(t-tD))/(1-exp(-KA*Tau))
    )
    list(CONC=CONC)
  }
}


# Getting parameters from parent frame ------------------------------------
pkFromParentFrame <- function(fun) {
  fAlgebraicName <- quote(fun)

  # Find that function, and its arguments
  myFun <- fun
  myFormals <- formals(myFun)

  # Create a wrapper function
  function(env=parent.frame()) {
    # Check for the arguments
    args <- mget(names(myFormals), envir=env, ifnotfound=myFormals)
    missing <- unlist( lapply(args, is.name) )
    if(any(missing)) {
      stop("Could not execute algebraic function ", fAlgebraicName,", ",
           "argument(s) ", paste(names(myFormals)[missing], collapse=", "), " not found."  )
    }

    do.call(myFun, args)
  }
}

#' Executes the requested PK model, and
#' fetches the arguments from the caller's environment.
#'
#' @param env environment in which to search the function arguments
#' filled in automatically
#'
#' @name pk
NULL

#' 
#' @rdname pk
pk1cptiv <- pkFromParentFrame(pk1cptiv_)
#' 
#' @rdname pk
pk1cptinfusion <- pkFromParentFrame(pk1cptinfusion_)
#' 
#' @rdname pk
pk1cptoral <- pkFromParentFrame(pk1cptoral_)
#' 
#' @rdname pk
pk2cptiv <- pkFromParentFrame(pk2cptiv_)
#' 
#' @rdname pk
pk2cptinfusion <- pkFromParentFrame(pk2cptinfusion_)
#' 
#' @rdname pk
pk2cptoral <- pkFromParentFrame(pk2cptoral_)
#' 
#' @rdname pk
pk3cptiv <- pkFromParentFrame(pk3cptiv_)
#' 
#' @rdname pk
pk3cptinfusion <- pkFromParentFrame(pk3cptinfusion_)
#' 
#' @rdname pk
pk3cptoral <- pkFromParentFrame(pk3cptoral_)


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
#' @param env environment in which to search for parameters
#'
#' 
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

#' Execute the appropriate PK model, based on the same rules as used by Monolix
#'
#' @details
#' See http://monolix.lixoft.com/data-and-models/pkmodel1/ for more information
#'
#' @param ... the parameters for the PK model
#' 
#' @importFrom rlang quos
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

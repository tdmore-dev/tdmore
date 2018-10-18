library(tdmore)
library(RxODE)
library(testthat)

context("Test that the algebraic library works as intented")

omegas=c(V=0.28^2, CL=0.19^2)
mod1 <- pk1cptivbolusVCL(THETA = list(V = 61, CL = 3.7),
                   OMEGA = omegas) %>% algebraic() %>% tdmore(exp = 0.30)

modelCode <- "
Vc = 61 * exp(ETA_V);
CL = 3.7 * exp(ETA_CL);
CONC = centr / Vc;
d/dt(centr) = -CL/Vc*centr;
"
mod2 <- RxODE::RxODE(modelCode) %>% tdmore(exp=0.30, omega=convertVectorToDiag(omegas))

regimen <- data.frame(
  TIME=seq(0, 7)*24,
  AMT=5 #5mg
)
observed <- data.frame(TIME=c(2, 28), CONC=c(0.04, 0.08))

ipred1 <- mod1 %>% estimate(observed, regimen)
p1 <- plot(ipred1, newdata=data.frame(TIME=seq(0, 34, length.out=100), CONC=NA))

ipred2 <- mod2 %>% estimate(observed, regimen)
p2 <- plot(ipred2, newdata=data.frame(TIME=seq(0, 34, length.out=100), CONC=NA))

# Check resuls are same (rxode vs algebraic)
expect_equal(round(ipred1$res, digits=4), round(ipred2$res, digits=4))

gridExtra::grid.arrange(p1 + ggtitle("Algebraic"), p2 + ggtitle("RxODE"))
cat("Algebraic method:\n")
print(summary(ipred1))
cat("RxODE method:\n")
print(summary(ipred2))

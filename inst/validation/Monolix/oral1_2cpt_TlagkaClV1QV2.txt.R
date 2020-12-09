m1 <- RxODE::RxODE("
V = V1;
k = Cl/V;
k12 = Q/V;
k21 = Q/V2;

d/dt(A0) = -ka*A0;
alag(A0) = Tlag;
d/dt(A1) = ka*A0 - k * A1 - k12*A1 + k21*A2;
d/dt(A2) = k12*A1 - k21*A2;
CONC = A1 / V;
")

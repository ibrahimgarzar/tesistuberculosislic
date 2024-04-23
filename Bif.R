setwd("Desktop/Tesis")

library(deSolve)
library(phaseR)

# declarar parametros
k2=1.68164470296576e-09 
k4=1677.71529516028
k5=3377.90445527686
k6=8881221.80265334
k9=4.89887206065021e-10
k10=0.0968962836928783
k13=0.000208989059062582
k14=168.675283030697
k1=0.00217688531461890
k3=0.0116821323710160
k8=2.67341167302623
k12=1.46921454819336
k15=23979612.7030956
k16=19.1664743109897
k17=9.23387740720957

# declarar los valores de parámetros
parms<- c(k2=1.68164470296576e-09, k4=1677.71529516028, k5=3377.90445527686, 
       k6=8881221.80265334, k9=4.89887206065021e-10, k10=0.0968962836928783, 
       k13=0.000208989059062582, k14=168.675283030697, k1=0.00217688531461890, 
       k3=0.0116821323710160, k8=2.67341167302623, k12=1.46921454819336, 
       k15=23979612.7030956, k16=19.1664743109897, k17=9.23387740720957)

# establecer modelo
Tuberculosis <- function(t, y, parms) {
  dM <- Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;  
  dMf <- M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
  dT <- k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
  dTf <- k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;
  list(c(dM,dMf,dT,dTf))
}













# tiempo de integración
tspan <- seq(from = 0, to = 10, by = 0.01)

source('Grind.r') 

model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
    dM <- M*k6-N*(O/k16)*k2*(1+M*k13*k17)-N*k1-N*T_t*k9;  
    dMf <- N*(1/k16)*O*k2*(1+M*k13*k17)-M*k3-M*(O*k10+k14)-M*k5-M*k13;
    dT <- k8*O*(1-O/k15)+O*k13*k16-N*(T_t/k16)*k2*(1+M*k13*k17);
    dTf <- k12*P*(1-P/(1+k15*M))+N*(1/k16)*O*k2*(1+M*k13*k17)-M*P*k4;
    return(list(c(dM,dMf,dT,dTf)))
  })
}

model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
    dM <- Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;  
    dMf <- M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
    dT <- k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
    dTf <- k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;
    return(list(c(dM,dMf,dT,dTf)))
  })
}

# declarar los valores de parámetros
p<- c(k2=1.68164470296576e-09, k4=1677.71529516028, k5=3377.90445527686, 
          k6=8881221.80265334, k9=4.89887206065021e-10, k10=0.0968962836928783, 
          k13=0.000208989059062582, k14=168.675283030697, k1=0.00217688531461890, 
          k3=0.0116821323710160, k8=2.67341167302623, k12=1.46921454819336, 
          k15=23979612.7030956, k16=19.1664743109897, k17=9.23387740720957)

# condición inicial
s <- c(x=0,y=0)

# graficar plano de fase
plane(xmax=4)

mid <- newton(s)
low <- newton(c(x=1,y=0))
hig <- newton(c(x=0,y=1))

continue(state=hig, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1) # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
# 15-64세 기준
# 성접촉으로 인힌한 감염만 고려 



# HIV mathematical model

Expansion <- function(year, state_values, parameters) {
  
  with(as.list(c(state_values, parameters)),{
    
    # Total population
    dN_M = (N_M + N_G) * kappa - zeta * S_M + tau_3 * (S_G + UD_G + ID_G + A_G) - d_1 * ( S_M + ID_M + UD_M + A_M) -  d_2 * A_M
    
    dN_W = N_W * kappa - d_1 * (S_W + ID_W + UD_W + A_W) -  d_2 * A_W
    
    dN_G = zeta * S_M  - tau_3 * (S_G + UD_G + ID_G + A_G) - d_1 * ( S_G + ID_G + UD_G + A_G) -  d_2 * A_G
    
    # Suceptible Population
    dS_M = (N_M + N_G) * kappa - zeta * S_M + tau_3 * S_G - (1-tau_1) * (1-gamma_1) * (1-omega_1) * eta_1 * phi * UD_W * S_M /N_W -
      (1-tau_1) * (1-gamma_1) * (1-omega_1) * eta_1 * (1-phi) * UD_W * S_M /N_W - d_1 * S_M
    
    dS_W = N_W * kappa - (1-tau_1) * (1-gamma_1) * (1-omega_2) * eta_1 * phi * UD_M * S_W /N_M -
      (1-tau_1) * (1-gamma_1) * (1-omega_2) * eta_1 * (1-phi) * UD_M * S_W /N_M - d_1 * S_W
    
    dS_G = zeta * S_M - tau_3 * S_G - (1-tau_1) * (1-tau_2) * (1-gamma_2) * (1-omega_3) * eta_2 * phi * UD_G * S_G /N_G -
      (1-tau_1) * (1-tau_2) * (1-gamma_2) * (1-omega_3) * eta_2 * (1-phi) * UD_G * S_G /N_G - d_1 * S_G
    
    
    #Infected Populatoin but Diagnosis
    
    dID_M = (1-tau_1) * (1-gamma_1) * (1-omega_1) * eta_1 * phi * UD_W * S_M /N_W + phi * UD_M + theta * A_M - beta * ID_M + tau_3 * ID_G - d_1 * ID_M
    
    dID_W = (1-tau_1) * (1-gamma_1) * (1-omega_2) * eta_1 * phi * UD_M * S_W /N_M + phi * UD_W + theta * A_W - beta * ID_W - d_1 * ID_M
    
    dID_G = (1-tau_1) * (1-tau_2) * (1-gamma_2) * (1-omega_3) * eta_2 * phi * UD_G * S_G /N_G + phi * UD_G + theta * A_G - beta * ID_M - tau_3 * ID_G - d_1 * ID_G
  
    
    #Infected Populatoin but Undiagnosis
    
    
    dUD_M = (1-tau_1) * (1-gamma_1) * (1-omega_1) * eta_1 * (1-phi) * UD_W * S_M /N_W - phi * UD_M - beta * ID_M + tau_3 * UD_G - d_1 * UD_M
    
    dUD_W = (1-tau_1) * (1-gamma_1) * (1-omega_2) * eta_1 * (1-phi) * UD_M * S_W /N_M - phi * UD_W  - beta * ID_W - d_1 * UD_W
    
    dUD_G = (1-tau_1) * (1-tau_2) * (1-gamma_2) * (1-omega_3) * eta_2 * (1-phi) * UD_G * S_G /N_G - phi * UD_G  - beta * UD_G - tau_3 * UD_G - d_1 * UD_G
    
    #AIDS Progress Populatoin 
    
    dA_M = beta * (UD_M + ID_M) - theta * A_M - (d_1+ d_2) * A_M + tau_3 * A_G 
    
    dA_W = beta * (UD_W + ID_W) - theta * A_W - (d_1+ d_2) * A_W
    
    dA_G = beta * (UD_G + ID_G) - theta * A_G - (d_1+ d_2) * A_G - tau_3 * A_G 
    
    
    # dD = d_1 *(S_M + S_W + S_G + ID_M + ID_W + ID_G + UD_M + UD_W + UD_G + A_M + A_W + A_G) + d_2 * (A_M + A_W + A_G)

    
    return(list(c(dN_M,dN_W,dN_G,dS_M,dS_W,dS_G,dID_M,dID_W,dID_G,dUD_M,dUD_W,dUD_G,dA_M,dA_W,dA_G)))
  })
}

# Load Data

# library(readr)
# data <- read_csv("Documents/GitHub/Mathematical-modeling/sample.csv")
# data<-as.data.frame(data)

# Load Data



param<-c(tau_1=0,tau_2=0,tau_3=0,theta = 0.2,beta=0.7,eta_1 = 0.1,eta_2 = 0.2, zeta = 0.005, kappa = 0.0065, phi=0.1,
         gamma_1=0,gamma_2=0,omega_1=0.1,omega_2=0.2,omega_3=0.3,
         d_1=0.005825,d_2=0.1642)


N_M = 22346173
N_W = 22196339
N_G = 11179

S_M= 22346129
S_W =22196358
S_G =11135

ID_M = 44
ID_W = 19
ID_G = 44

UD_M = 0
UD_W = 0
UD_G = 0

A_M = 0
A_W = 0
A_G = 0

D=0



init<-c(N_M,N_W,N_G,S_M,S_W,S_G,ID_M,ID_W,ID_G,UD_M,UD_W,UD_G,A_M,A_W,A_G)

times <-seq(1995,2018,1)

library(deSolve)

out <-ode(y=init, parms=param, times = times, func = Expansion)





HIV.model = function(data, parameters){ 
  
  G0 <- as.vector(0.005 * data$total_M)
  Ym0 <- as.vector(data$total_M) - G0
  Yw0 <- as.vector(data$total_W)
  IM0 <- as.vector(data$Men)
  IW0 <- as.vector(data$Women)
  IG0 <- as.vector(data$Gay)
  
  param<-c(tau_1=0,tau_2=0,tau_3=0,theta = 0,beta=0,eta_1 = 0,eta_2 = 0, zeta = 0, kappa = 0, phi=0,
           gamma1=0,gamma2=0,omega_1=parameters[2],omega_2=parameters[3],omega_3=parameters[3],
           d_1=0.007165862,d_2=parameters[5])
  
  Ym= 23842278
  Yw =23770470
  G =119810.4
  Im = 58
  Iw = 25
  Ig = 136
  
  
  
  init<-c(Ym= Ym,Yw =Yw, G =G,
          Im = Im,Iw = Iw,Ig = Ig)
  
  time<-seq(2000,2018,1)
  
  Ym0_predicted <- as.vector(ode(init, time, HIV,param)[,2])
  Yw0_predicted <- as.vector(ode(init, time, HIV,param)[,3])
  G0_predicted <- as.vector(ode(init, time, HIV,param)[,4])
  IM0_predicted <- as.vector(ode(init, time, HIV,param)[,5])
  IW0_predicted <- as.vector(ode(init, time, HIV,param)[,6])
  IG0_predicted <- as.vector(ode(init, time, HIV,param)[,7])
  
  
  #error1 <- Ym0 - YM_predicted 
  #error2<- Yw0 - Ym0_predicted 
  #error3<- G0 - G0_predicted 
  error4<- IM0 - IM0_predicted 
  error6<- IW0 - IW0_predicted 
  error8<- IG0 - IG0_predicted 
  #error10<- D0 - D0_predicted 
  
  
  #SSE1 = sum(error1^2)
  #SSE2 = sum(error2^2)
  #SSE3 = sum(error3^2)
  SSE4 = sum(error4^2)
  SSE6 = sum(error6^2)
  SSE8 = sum(error8^2)
  #SSE10 = sum(error10^2)
  
  return(SSE4)
}



library(deSolve)




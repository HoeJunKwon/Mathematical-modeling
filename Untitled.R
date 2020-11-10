library(deSolve)
library(GenSA)

# sex oriented HIV Transmission mathematical model

sex_oriented <- function(year, state_values, parameters) {
  with(as.list(c(state_values, parameters)),{
    
    #Total population
    # N_M = (S_M + I_M + T_M)
    # N_W = (S_W + I_W + T_W)
    # N_MSM= (S_MSM + I_MSM+T_MSM)
    
    # transmission_rate
    
    # lambda_M = beta_1 * gamma_1 * ( 1 - omega_1 ) * ( 1- kappa_1 ) * (I_W + theta_1 * T_W) / (S_M + I_M + T_M)
    # lambda_W = beta_2 * gamma_2 * ( 1 - omega_1 ) * ( 1- kappa_1 ) * (I_M + theta_2 * T_M)  / (S_W + I_W + T_W)
    # lambda_MSM = beta_3 * gamma_3 * ( 1 - omega_1 ) * ( 1 - omega_2) * ( 1- kappa_2 ) * (I_MSM + theta_3 * T_MSM)  / (S_MSM + I_MSM + T_MSM)
    
    # Suceptible Population
    
    dS_M = Pi * (1 - tau) * psi1 + rho * S_MSM - beta_1 * gamma_1 * ( 1 - omega_1 ) * ( 1- kappa_1 ) * (I_W + theta_1 * T_W) / (S_M + I_M + T_M) * S_M - ( mu_d * zeta ) * S_M
    
    dS_W = Pi * tau * psi2 - beta_2 * gamma_2 * ( 1 - omega_1 ) * ( 1- kappa_1 ) * (I_M + theta_2 * T_M)  / (S_W + I_W + T_W) * S_W - mu_d * S_W
    
    dS_MSM = zeta * S_M -  beta_3 * gamma_3 * ( 1 - omega_1 ) * ( 1 - omega_2) * ( 1- kappa_2 ) * (I_MSM + theta_3 * T_MSM)  / (S_MSM + I_MSM + T_MSM) * S_MSM - (mu_d * rho) * S_MSM
    
    
    #Infected Population 
    
    dI_M = beta_1 * gamma_1 * ( 1 - omega_1 ) * ( 1- kappa_1 ) * (I_W + theta_1 * T_W) / (S_M + I_M + T_M)* S_M - (mu_d + delta_a) * I_M
    
    dI_W = beta_2 * gamma_2 * ( 1 - omega_1 ) * ( 1- kappa_1 ) * (I_M + theta_2 * T_M)  / (S_W + I_W + T_W) * S_W - (mu_d + delta_a) * I_W
    
    dI_MSM =  beta_3 * gamma_3 * ( 1 - omega_1 ) * ( 1 - omega_2) * ( 1- kappa_2 ) * (I_MSM + theta_3 * T_MSM)  / (S_MSM + I_MSM + T_MSM) * S_MSM - (mu_d + delta_a) * I_MSM
    
    
    #Treated Populatoio but 
    
    
    dT_M = delta_a * I_M - (mu_d) * T_M
    
    dT_W = delta_a * I_W - (mu_d) * T_W
    
    dT_MSM = delta_a * I_MSM - (mu_d) * T_MSM
    
    return(list(c(dS_M,dS_W,dS_MSM,dI_M,dI_W,dI_MSM,dT_M,dT_W,dT_MSM)))
  })
}

Total_Men <-c(19084999,19142958,19192036,19263033,19305761,19230937,19183985)#,19303494) # 2012~2019 15-64
Total_Women <-c(18258526,18314400,18360874,18429691,18478656,18405536,18363056)#,18286058) # 2012~2019 15-64

Survival_Men <-c(6691,7436,8248,9010,9788,10489,11236)#,11938) # 2012~2019 생존감염자  15-64
Survival_Women <-c(552,593,624,658,690,716,742)#,764) # 2012~2019 생존감염자  15-64

Men_Report <-c(776,912,974,937,957,913,896)#,922) # 2012~2019 Men Report  15-64
Women_Report <-c(56,59,56,42,53,45,42)#,45) # 2012~2019 생존감염자  15-64




N_M0 = 18935443
S_M0 = 18934025
I_M0 =  326
T_M0 = 1092

N_W0 = 18258526
S_W0 = 18257974
I_W0 = 127
T_W0 = 425


N_MSM0 = 149556
S_MSM0 = 144283
I_MSM0 =  1213
T_MSM0 = 4060




# for (i in seq(1:length(Survival_Men)-1)) {
#   
#   AIDS_Death_Men[i]   <- Survival_Men[i] + Men_Report[i]  - Survival_Men[i+1]
#   AIDS_Death_Women[i] <- Survival_Women[i] + Women_Report[i] - Survival_Women[i+1]
#   
# }
# 
# AIDS_Death_Men   <- AIDS_Death_Men[-8]   ; AIDS_Death_Men # 2012~2018 AIDS 사망자
# AIDS_Death_Women <- AIDS_Death_Women[-8] ; AIDS_Death_Women # 2012~2018 AIDS 사망자


derivs <-function(parameters){
  
  
  Years <- seq(2012,2018,1)
  
  init <- c(S_M=S_M0,S_W=S_W0,S_MSM=S_MSM0,I_M=I_M0,I_W=I_W0,I_MSM=I_MSM0,T_M=T_M0,T_W=T_W0,T_MSM=T_MSM0 )
  
  
  parameter <- c(beta_1=0.004,beta_2=0.0247,beta_3=0.03674,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
                 mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.77, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
  
  I_M0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,5])
  I_W0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,6])
  I_MSM0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,7])
  T_M0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,8])
  T_W0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,9])
  T_MSM0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,10])
  
  
  MSE1<- sum(Survival_Men - ((I_M0_predicted + T_M0_predicted)+(I_MSM0_predicted + T_MSM0_predicted)))^2/length(Years)
  MSE2<- sum(Survival_Women - (I_W0_predicted + T_W0_predicted))^2/length(Years)


  
  MSE <- (MSE1 + MSE2) / 2
  
  return(MSE)
}

## GenSA

dimension <- 3

lower <- rep(0, dimension)
upper <- rep(1, dimension)

out <- GenSA(lower = lower, upper = upper, fn = derivs, control=list(max.time = 60,verbose=TRUE))



out[c("value","par","counts")]
ou1<- as.data.frame(out[c("trace.mat")])

## Genetic Algorithsm




Years <- seq(2012,2018,1)

init <- c(S_M=S_M0,S_W=S_W0,S_MSM=S_MSM0,I_M=I_M0,I_W=I_W0,I_MSM=I_MSM0,T_M=T_M0,T_W=T_W0,T_MSM=T_MSM0 )


parameter <- c(beta_1=0.004,beta_2=0.0247,beta_3=0.03674,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.77, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)

ode(init, Years, sex_oriented,parameter)

Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]

Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

Men
Women

####

library(readr)
Data <-as.data.frame(read_csv("~/GitHub/Data.csv"))
Data1 <- Data[seq(1:14),]

Data2 <- Data[c(15:28),]

plot<-ggplot(data = Data2,aes(x=Year,y=population))+geom_point(aes(color=Name))
plot
 


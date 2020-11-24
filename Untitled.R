library(deSolve)
library(GenSA)
library(ggplot2)
library(readr)
# library(grid)


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




# # for (i in seq(1:length(Survival_Men)-1)) {
# #   
# #   AIDS_Death_Men[i]   <- Survival_Men[i] + Men_Report[i]  - Survival_Men[i+1]
# #   AIDS_Death_Women[i] <- Survival_Women[i] + Women_Report[i] - Survival_Women[i+1]
# #   
# # }
# # 
# # AIDS_Death_Men   <- AIDS_Death_Men[-8]   ; AIDS_Death_Men # 2012~2018 AIDS 사망자
# # AIDS_Death_Women <- AIDS_Death_Women[-8] ; AIDS_Death_Women # 2012~2018 AIDS 사망자
# 
# 
# derivs <-function(parameters){
#   
#   
#   Years <- seq(2012,2018,1)
#   
#   init <- c(S_M=S_M0,S_W=S_W0,S_MSM=S_MSM0,I_M=I_M0,I_W=I_W0,I_MSM=I_MSM0,T_M=T_M0,T_W=T_W0,T_MSM=T_MSM0 )
#   
#   
#   parameter <- c(beta_1=0.004,beta_2=0.0247,beta_3=0.03674,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
#                  mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.77, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
#   
#   I_M0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,5])
#   I_W0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,6])
#   I_MSM0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,7])
#   T_M0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,8])
#   T_W0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,9])
#   T_MSM0_predicted <- as.vector(ode(init, Years, sex_oriented,parameter)[,10])
#   
#   
#   MSE1<- sum(Survival_Men - ((I_M0_predicted + T_M0_predicted)+(I_MSM0_predicted + T_MSM0_predicted)))^2/length(Years)
#   MSE2<- sum(Survival_Women - (I_W0_predicted + T_W0_predicted))^2/length(Years)
# 
# 
#   
#   MSE <- (MSE1 + MSE2) / 2
#   
#   return(MSE)
# }
# 
# ## GenSA
# 
# dimension <- 3
# 
# lower <- rep(0, dimension)
# upper <- rep(1, dimension)
# 
# out <- GenSA(lower = lower, upper = upper, fn = derivs, control=list(max.time = 60,verbose=TRUE))
# 
# 
# 
# out[c("value","par","counts")]
# ou1<- as.data.frame(out[c("trace.mat")])





Years <- seq(2012,2018,1)


init <- c(S_M=S_M0,S_W=S_W0,S_MSM=S_MSM0,I_M=I_M0,I_W=I_W0,I_MSM=I_MSM0,T_M=T_M0,T_W=T_W0,T_MSM=T_MSM0 )

parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)

ode(init, Years, sex_oriented,parameter)

Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]

Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

Men
Women

(mean((Men/Survival_Men)[1:6]) + (Survival_Men/Men)[7])/2
(mean((Women/Survival_Women)[1:6]) + (Survival_Women/Women)[7])/2


sample <-as.data.frame(read_csv("C:/Users/Monokuma/Desktop/sample.csv"))

ggplot(data=sample, aes(x=Sex,y= Ratio,fill = Sex))+
  geom_boxplot(fill='slategrey',color='darkslategrey',width=0.3)+
theme_classic(base_size = 20)





####

library(readr)
Data <- read_csv("Documents/GitHub/Mathematical-modeling/Data.csv")
Data1 <- Data[seq(1:14),]

Data2 <- Data[c(15,28),]


S_Men <-ode(init, Years, sex_oriented,parameter)[,2] +ode(init, Years, sex_oriented,parameter)[,4]
S_Women <-ode(init, Years, sex_oriented,parameter)[,3] 
# 
# write.csv(S_Men,"/Users/victoria/Downloads/S_Men.csv")
# write.csv(S_Women,"/Users/victoria/Downloads/S_Women.csv")
# 



Data3 <- Data[seq(15,21),]
Data4 <- Data[seq(22,28),]
Data5 <- Data[seq(29,35),]
Data6 <- Data[seq(36,42),]

plot_base1 <-ggplot(data = Data3,aes(x=Year,y=population))+geom_point(aes(color=Name),size = 0.5) + 
  geom_line(data = Data3,aes(x=Year,y=population))+geom_point(aes(color=Name), size=0.5) 

plot_base <- plot_base1 + geom_point(data = Data4,aes(x=Year,y=population,color=Name),size = 0.5) + theme_classic()  +
  theme(legend.position = c(0.15, 0.8),legend.background=element_rect(fill="white", color="black"),legend.text = element_text(size=13),
        axis.title=element_text(size=15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  guides(color=guide_legend(title = NULL)) + geom_line(data = Data4,aes(x=Year,y=population,color=Name),size = 0.5)


q <-plot_base  + geom_point(data = Data5,aes(x=Year,y=population,color=Name)) + 
  geom_line(data = Data5,aes(x=Year,y=population,color=Name),size = 0.5)  +
  geom_point(data = Data6,aes(x=Year,y=population,color=Name)) + 
  geom_line(data = Data6,aes(x=Year,y=population,color=Name),size = 0.5)








## NO prevention effect 
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")



Data <-as.data.frame(read_csv("/Users/victoria/Desktop/Infected.csv"))

# fitting data
Data1 <- Data[c(1:14),]
Data2 <- Data[c(15:28),]

plot<-ggplot(data = Data,aes(x=Year,y=population))+geom_point(aes(color=Name))
plot1 <- plot + geom_line(data = Data2,aes(x=Year,y=population,color=Name))
plot1


## Current Situation

Data1 <- Data[c(1:14),]

plot_base1 <-ggplot(data = Data1[c(1:7),],aes(x=Year,y=population))+geom_point(aes(color=Name),size = 0.5) + 
  geom_line(size=0.5,color = "red") 

plot_base <- plot_base1 + geom_point(data = Data1[c(8:14),],aes(x=Year,y=population,color=Name),size = 0.5) + theme_classic()  +
  theme(legend.position = c(0.15, 0.8),legend.background=element_rect(fill="white", color="black"),legend.text = element_text(size=13),
        axis.title=element_text(size=15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  guides(color=guide_legend(title = NULL)) + geom_line(data = Data1[c(8:14),],aes(x=Year,y=population,color=Name),size = 0.5)


plot_base  + geom_point(data = Data2[c(1:7),],aes(x=Year,y=population,color=Name)) + 
  geom_line(data = Data2[c(1:7),],aes(x=Year,y=population,color=Name),size = 0.5)  +
  geom_point(data = Data2[c(8:14),],aes(x=Year,y=population,color=Name)) + 
  geom_line(data = Data2[c(8:14),],aes(x=Year,y=population,color=Name),size = 0.5)




# No Prevention

Data3 <- Data[c(29:49),]
Data4 <- Data[c(50:70),]

plot <-plot_base +geom_point(data = Data3,aes(x=Year,y=population,color=Name),size = 0.5) + 
  geom_vline(xintercept = as.numeric(Data3$Year[1]), color = "red", linetype = 2) + 
  geom_line(data = Data3,aes(x=Year,y=population,color=Name),size=0.5)

plot_prev <- plot + geom_point(data = Data4,aes(x=Year,y=population,color=Name),size = 0.5) + 
  geom_line(data = Data4,aes(x=Year,y=population,color=Name),size=0.5) + 
  theme(legend.position = c(0.125, 0.8)) +
  scale_x_continuous(breaks=seq(2012,2038,2))

  
plot_prev


####################
####################

## Heterosexual --> 0.758% 고정

## Only PrEP 50%(MSM) --> 43.137%
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0.00758,kappa_2=0.43137,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

# write.csv(Men,"/Users/victoria/Downloads/Men.csv")
# write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data5 <- Data[c(71:91),]
Data6 <- Data[c(92:112),]

plot_PrEP1 <- plot_prev + geom_point(data = Data5,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data5,aes(x=Year,y=population,color=Name),size=0.5) 

plot_PrEP_C1 <- plot_PrEP1 + geom_point(data = Data6,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data6,aes(x=Year,y=population,color=Name),size=0.5) 

plot_PrEP_C1



C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.5* 5011450

C1 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C1) <-c("Men","Women","MSM")
write.csv(C1,"/Users/victoria/Downloads/C1.csv")
## ## Only PrEP 60%(MSM) --> 51.765%

Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0.00758,kappa_2=0.51765,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

# write.csv(Men,"/Users/victoria/Downloads/Men.csv")
# write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data7 <- Data[c(113:133),]
Data8 <- Data[c(134:154),]

plot_PrEP2 <- plot_PrEP_C1 + geom_point(data = Data7,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data7,aes(x=Year,y=population,color=Name),size=0.5) 

plot_PrEP_C2 <- plot_PrEP2 + geom_point(data = Data8,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data8,aes(x=Year,y=population,color=Name),size=0.5) 

plot_PrEP_C2


C2A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C2B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C2C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.6* 5011450

C2 <-as.data.frame(cbind(C2A,C2B,C2C))
colnames(C2) <-c("Men","Women","MSM")
write.csv(C2,"/Users/victoria/Downloads/C2.csv")

## ## Only PrEP 70%(MSM) --> 60.392%
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0.00758,kappa_2=0.60392,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

# write.csv(Men,"/Users/victoria/Downloads/Men.csv")
# write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data9 <- Data[c(155:175),]
Data10 <- Data[c(176:196),]

plot_PrEP3 <- plot_PrEP_C2 + geom_point(data = Data9,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data9,aes(x=Year,y=population,color=Name),size=0.5) 

plot_PrEP_C3 <- plot_PrEP3 + geom_point(data = Data10,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data10,aes(x=Year,y=population,color=Name),size=0.5) 

plot_PrEP_C3


C3A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C3B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C3C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.7* 5011450

C3 <-as.data.frame(cbind(C3A,C3B,C3C))
colnames(C3) <-c("Men","Women","MSM")
write.csv(C3,"/Users/victoria/Downloads/C3.csv")


## ## Only PrEP 80%(MSM) --> 69.02%
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0,omega_2=0,kappa_1=0.00758,kappa_2=0.6902,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

# write.csv(Men,"/Users/victoria/Downloads/Men.csv")
# write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data11 <- Data[c(197:217),]
Data12 <- Data[c(218:238),]

plot_PrEP4 <- plot_PrEP_C3 + geom_point(data = Data11,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data11,aes(x=Year,y=population,color=Name),size=0.5) 

plot_PrEP_C4 <- plot_PrEP4 + geom_point(data = Data12,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data12,aes(x=Year,y=population,color=Name),size=0.5) + theme(legend.position = c(0.1, 0.8))

plot_PrEP_C4


C4A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C4B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C4C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.8* 5011450

C4 <-as.data.frame(cbind(C4A,C4B,C4C))
colnames(C4) <-c("Men","Women","MSM")
write.csv(C4,"/Users/victoria/Downloads/C4.csv")

###########################################
###########################################
###########################################


# 100% 인 경우
# rho = 0.02
# omega1 =  0.78
# omega2 =  0.33


## 50% Education omega1 =  0.39, omega2 =  0.165
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.39,omega_2=0.165,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0.02,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]
write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data13 <- Data[c(239:259),]
Data14 <- Data[c(260:280),]

plot_Edu1 <- plot_prev + geom_point(data = Data13,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data13,aes(x=Year,y=population,color=Name),size=0.5) 
plot_Edu_C1 <- plot_Edu1 + geom_point(data = Data14,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data14,aes(x=Year,y=population,color=Name),size=0.5) 
plot_Edu_C1

## 60% Education omega1 =  0.468, omega2 =  0.198
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.468,omega_2=0.198,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0.02,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]
write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data15 <- Data[c(281:301),]
Data16 <- Data[c(302:322),]

plot_Edu2 <- plot_Edu_C1 + geom_point(data = Data15,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data15,aes(x=Year,y=population,color=Name),size=0.5) 
plot_Edu_C2 <- plot_Edu2 + geom_point(data = Data16,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data16,aes(x=Year,y=population,color=Name),size=0.5)  
plot_Edu_C2

## 70% Education omega1 =  0.546, omega2 =  0.231
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.546,omega_2=0.231,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0.02,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]
write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data17 <- Data[c(323:343),]
Data18 <- Data[c(344:364),]

plot_Edu3<- plot_Edu_C2 + geom_point(data = Data17,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data17,aes(x=Year,y=population,color=Name),size=0.5)   
plot_Edu_C3 <- plot_Edu3 + geom_point(data = Data18,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data18,aes(x=Year,y=population,color=Name),size=0.5)  
plot_Edu_C3

## 80% Education omega1 =  0.624, omega2 =  0.264

Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.624,omega_2=0.264,kappa_1=0,kappa_2=0,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0.02,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]
write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data19 <- Data[c(365:385),]
Data20 <- Data[c(386:406),]

plot_Edu4 <- plot_Edu_C3 + geom_point(data = Data19,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data19,aes(x=Year,y=population,color=Name),size=0.5)  
plot_Edu_C4 <- plot_Edu4 + geom_point(data = Data20,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data20,aes(x=Year,y=population,color=Name),size=0.5) + theme(legend.position = c(0.1, 0.8))
plot_Edu_C4


###########################################
###########################################
###########################################

## 80% Education omega1 =  0.624, omega2 =  0.264
## ## Only PrEP 80%(MSM) --> 69.02%

Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.624,omega_2=0.264,kappa_1=0.00758,kappa_2=0.6902,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0.02,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]
write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

Data21 <- Data[c(407:427),]
Data22 <- Data[c(428:448),]

plot_Conv <- plot_prev + geom_point(data = Data21,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data21,aes(x=Year,y=population,color=Name),size=0.5)  
plot_Conv_C <- plot_Conv + geom_point(data = Data22,aes(x=Year,y=population,color=Name),size = 0.5) +
  geom_line(data = Data22,aes(x=Year,y=population,color=Name),size=0.5) + theme(legend.position = c(0.1, 0.8))
plot_Conv_C


#2013
# 9966415 + 5011450 * 80000
# 
# x = 9966415
# k<-c()
# for (i in 1:25) {
#   x = x+x * 0.0441294
#   k[i+1]<- x
# }
# k
# cost <-k[7:26]
# Data1 <-as.data.frame(read_csv("/Users/victoria/Desktop/infected2.csv"))
# An_cost <-Data1[3]*0.8 *cost
# write.csv(An_cost,"/Users/victoria/Downloads/cost.csv")
# 
# 
# 
# t<-20000000000
# ed<-c()
# for (i in 1:25) {
#   t = t+t * 0.0441294
#   ed[i]<- t
# }
# ed
# write.csv(ed,"/Users/victoria/Downloads/ed.csv")








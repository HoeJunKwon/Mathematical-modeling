#PrEP <-c(50,60,70,80)
#Edu <-c(50,60,70,80)

## Only PrEP 50%(MSM) --> 0.43137
## Only PrEP 60%(MSM) --> 0.51765
## Only PrEP 70%(MSM) --> 0.60392
## Only PrEP 80%(MSM) --> 0.6902
## 50% Education omega_1 =  0.39, omega_2 =  0.165
## 60% Education omega_1 =  0.468, omega_2 =  0.198
## 70% Education omega_1 =  0.546, omega_2 =  0.231
## 80% Education omega_1 =  0.624, omega_2 =  0.264


# PrEP = 50 , Edu = 50
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.39,omega_2=0.165,kappa_1=0.00758,kappa_2=0.43137,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")



C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.5* 5011450


C1 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C1) <-c("Men","Women","MSM")
write.csv(C1,"/Users/victoria/Downloads/C/C1.csv")


# PrEP = 50 , Edu = 60
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.468,omega_2=0.198,kappa_1=0.00758,kappa_2=0.43137,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.5* 5011450

C2 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C2) <-c("Men","Women","MSM")
write.csv(C2,"/Users/victoria/Downloads/C/C2.csv")

# PrEP = 50 , Edu = 70
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.546,omega_2=0.231,kappa_1=0.00758,kappa_2=0.43137,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.5* 5011450

C3 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C3) <-c("Men","Women","MSM")
write.csv(C3,"/Users/victoria/Downloads/C/C3.csv")

# PrEP = 50 , Edu = 80
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1=0.624,omega_2=0.264,kappa_1=0.00758,kappa_2=0.43137,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.5* 5011450

C4 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C4) <-c("Men","Women","MSM")
write.csv(C4,"/Users/victoria/Downloads/C/C4.csv")


# PrEP = 60 , Edu = 50
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.39, omega_2 =  0.165,kappa_1=0.00758,kappa_2=0.51765,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.6* 5011450

C5 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C5) <-c("Men","Women","MSM")
write.csv(C5,"/Users/victoria/Downloads/C/C5.csv")

# PrEP = 60 , Edu = 60
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.468, omega_2 =  0.198,kappa_1=0.00758,kappa_2=0.51765,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.6* 5011450

C6 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C6) <-c("Men","Women","MSM")
write.csv(C6,"/Users/victoria/Downloads/C/C6.csv")

# PrEP = 60 , Edu = 70
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.546, omega_2 =  0.231,kappa_1=0.00758,kappa_2=0.51765,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")


C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.6* 5011450

C7 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C7) <-c("Men","Women","MSM")
write.csv(C7,"/Users/victoria/Downloads/C/C7.csv")

# PrEP = 60 , Edu = 80
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.624, omega_2 =  0.264,kappa_1=0.00758,kappa_2=0.51765,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.6* 5011450

C8 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C8) <-c("Men","Women","MSM")
write.csv(C8,"/Users/victoria/Downloads/C/C8.csv")

# PrEP = 70 , Edu = 50
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.39, omega_2 =  0.165,kappa_1=0.00758,kappa_2=0.60392,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.7* 5011450

C9 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C9) <-c("Men","Women","MSM")
write.csv(C9,"/Users/victoria/Downloads/C/C9.csv")

# PrEP = 70 , Edu = 60
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.468, omega_2 =  0.198,kappa_1=0.00758,kappa_2=0.60392,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 * 0.01 * 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 * 0.01 * 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 * 0.7 * 5011450

C10 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C10) <-c("Men","Women","MSM")
write.csv(C10,"/Users/victoria/Downloads/C/C10.csv")

# PrEP = 70 , Edu = 70
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.546, omega_2 =  0.231,kappa_1=0.00758,kappa_2=0.60392,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.7* 5011450

C11 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C11) <-c("Men","Women","MSM")
write.csv(C11,"/Users/victoria/Downloads/C/C11.csv")

# PrEP = 70 , Edu = 80
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.624, omega_2 =  0.264,kappa_1=0.00758,kappa_2=0.60392,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.7* 5011450

C12 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C12) <-c("Men","Women","MSM")
write.csv(C12,"/Users/victoria/Downloads/C/C12.csv")

# PrEP = 80 , Edu = 50
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.39, omega_2 =  0.165,kappa_1=0.00758,kappa_2=0.6902,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.8* 5011450

C13 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C13) <-c("Men","Women","MSM")
write.csv(C13,"/Users/victoria/Downloads/C/C13.csv")

# PrEP = 80 , Edu = 60
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.468, omega_2 =  0.198,kappa_1=0.00758,kappa_2=0.6902,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.8* 5011450

C14 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C14) <-c("Men","Women","MSM")
write.csv(C14,"/Users/victoria/Downloads/C/C14.csv")

# PrEP = 80 , Edu = 70
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.546, omega_2 =  0.231,kappa_1=0.00758,kappa_2=0.6902,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.8* 5011450

C15 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C15) <-c("Men","Women","MSM")
write.csv(C15,"/Users/victoria/Downloads/C/C15.csv")

# PrEP = 80 , Edu = 80
Years <- seq(2018,2038,1)
init <- c(S_M=19079868,S_W=18361009,S_MSM=148919.7,I_M=180.2634,I_W=203.6377,I_MSM=4095.125,T_M=1253.061,T_W=538.3682,T_MSM=5707.558 )
parameter <- c(beta_1=0.004,beta_2=0.01765,beta_3=0.033786,gamma_1=3,gamma_2=3,gamma_3=6,omega_1 =  0.624, omega_2 =  0.264, kappa_1=0.00758,kappa_2=0.6902,Pi = 19084999 + 18258526,tau=0.499246,
               mu_d=0.00005825,zeta = 0.0000804, delta_a = 0.11, rho=0,psi1 = 0.0013,psi2 = 0.00098,theta_1=0.1,theta_2=0.3,theta_3 =0.3)
ode(init, Years, sex_oriented,parameter)
Men <-ode(init, Years, sex_oriented,parameter)[,5] +ode(init, Years, sex_oriented,parameter)[,7]+
  ode(init, Years, sex_oriented,parameter)[,8]+ode(init, Years, sex_oriented,parameter)[,10]
Women <-ode(init, Years, sex_oriented,parameter)[,6] +ode(init, Years, sex_oriented,parameter)[,9]

write.csv(Men,"/Users/victoria/Downloads/Men.csv")
write.csv(Women,"/Users/victoria/Downloads/Women.csv")

C1A <-(ode(init, Years, sex_oriented,parameter)[,5] + ode(init, Years, sex_oriented,parameter)[,8])[2:21] * 3 *0.01* 5011450
C1B <-(ode(init, Years, sex_oriented,parameter)[,6] + ode(init, Years, sex_oriented,parameter)[,9])[2:21]  * 3 *0.01* 5011450
C1C <-(ode(init, Years, sex_oriented,parameter)[,7] + ode(init, Years, sex_oriented,parameter)[,10])[2:21]  * 6 *0.8* 5011450

C16 <-as.data.frame(cbind(C1A,C1B,C1C))
colnames(C16) <-c("Men","Women","MSM")
write.csv(C16,"/Users/victoria/Downloads/C/C16.csv")



















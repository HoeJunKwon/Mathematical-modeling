HIV <- function(year, state_values, parameters) {

  with(as.list(c(state_values, parameters)),{
    
    dYm<-  M-(1-tau3)*zeta*Ym+tau3*G-omega1*(1-gamma1)*(1-tau1)*((Iw)/(Iw+Yw))*Ym-d1*Ym
    
    dYw<-  W-omega2*(1-gamma1)*(1-tau1)*((Im)/(Im+Ym))*Yw-d1*Yw
    
    dG<-  (1-tau3)*zeta*Ym-tau3*G-(1-gamma2)*(1-tau2)*(omega4)*(Ig/(Ig+G))*G-d1*G
    
    dIm<- omega1*(1-gamma1)*((Iw)/(Iw+Yw))*Ym-d2*Im
    
    dIw<- omega2*(1-gamma1)*(1-tau1)*((Im)/(Im+Ym))*Yw-d2*Iw                                                                 
    
    dIg<- (1-gamma2)*(1-tau2)*(omega4)*(Ig/(Ig+G))*G-d2*Ig               
    
    dD <- d1*Ym+d1*Yw+d1*G+d2*Im+d2*Iw+d2*Ig
    
    return(list(c(dYm, dYw, dG, dIm, dIw, dIg, dD)))
  })
}


HIV.model = function(data, parameters){ 
  
  Ym0 <- as.vector(data$Ym)
  Yw0 <- as.vector(data$Yw)
  G0 <- as.vector(data$G)
  IM0 <- as.vector(data$IM)
  IW0 <- as.vector(data$IW)
  IG0 <- as.vector(data$IG)
  D0<-as.vector(data$D)
  
  param<-c(tau1=0,tau2=0,tau3=0,tau4=0,zeta=parameters[1],
           gamma1=0,gamma2=0,omega1=parameters[2],omega2=parameters[3],
           omega4=parameters[4],d1=0.007165862,d2=parameters[5],M=246538,W=257026)
  
  
  init<-c(Ym= 16380526,Yw =16576010, G =131045,
          Im = 1951,Iw = 7,Ig = 5854,D = 236162)
  
  time<-seq(1992,2018,1)
  
  Ym0_predicted <- as.vector(ode(init, time, HIV,param)[,2])
  Yw0_predicted <- as.vector(ode(init, time, HIV,param)[,3])
  G0_predicted <- as.vector(ode(init, time, HIV,param)[,4])
  IM0_predicted <- as.vector(ode(init, time, HIV,param)[,5])
  IW0_predicted <- as.vector(ode(init, time, HIV,param)[,6])
  IG0_predicted <- as.vector(ode(init, time, HIV,param)[,7])
  D0_predicted<-as.vector(ode(init, time, HIV,param)[,8])
  
  
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


setwd("~/Desktop")
library(readr)
data <- read_csv("fit.csv")
data<-as.data.frame(data)
data


new.seed = as.integer(runif(1)*2e9) 
cat("Random seed: ", new.seed, "\n")
set.seed(new.seed) 

par = runif(5)


min.global = 0      # Expected global minimum
tol = 1e-13         # Tolerance

# 6. Set the search range of each parameter ###########################

dimension=5
# 
lower <- rep(0,dimension) 
upper <- rep(1,dimension) 

lower <- rep(0.03,0.000001,0.000001,0.000001,0.000005) 
upper <- rep(0.07,0.01,0.01,0.1,0.000007) 




# 7. Run SA to estimate (optimize) the unknown parameters #############
library(GenSA)
library(deSolve)
#################################################################################

par<-c(0.05,0.00144,0.00288,0.04986,5.5047E-05)

example<-matrix(0, nrow = 1, ncol = 5)
example<-as.data.frame(example)
colnames(example)<-c("par1","par2","par3","par4","par5")


out <- GenSA(par = par, lower = lower, upper = upper, fn = HIV.model, data = data,
             control=list(threshold.stop=min.global+tol,
                          verbose=FALSE))


####SSE4
K<-as.vector(out$par)
example1<-rbind(example,K)
example1<-example1[2,]
example1



HIV <- function(year, state_values, parameters) {

  with(as.list(c(state_values, parameters)),{
    
    dYm<-  M-(1-tau3)*zeta*Ym+tau3*G-omega1*(1-gamma1)*(1-tau1)*((Iw)/(Iw+Yw))*Ym-d1*Ym
    
    dYw<-  W-omega2*(1-gamma1)*(1-tau1)*((Im)/(Im+Ym))*Yw-d1*Yw
    
    dG<-  (1-tau3)*zeta*Ym-tau3*G-(1-gamma2)*(1-tau2)*(omega4)*(Ig/(Ig+G))*G-d1*G
    
    dIm<- omega1*(1-gamma1)*((Iw)/(Iw+Yw))*Ym-d2*Im
    
    dIw<- omega2*(1-gamma1)*(1-tau1)*((Im)/(Im+Ym))*Yw-d2*Iw                                                                 
    
    dIg<- (1-gamma2)*(1-tau2)*(omega4)*(Ig/(Ig+G))*G-d2*Ig               
    
    return(list(c(dYm, dYw, dG, dIm, dIw, dIg)))
  })
}

library(readr)
data <- read_csv("sample.csv")
data<-as.data.frame(data)
data



HIV.model = function(data, parameters){ 
  
  Ym0 <- as.vector(data$total_M - 0.005 * data$total_M )
  Yw0 <- as.vector(data$total_W)
  G0 <- as.vector(0.005 * data$total_M )
  IM0 <- as.vector(data$Men)
  IW0 <- as.vector(data$Women)
  IG0 <- as.vector(data$Gay)
  
  param<-c(tau1=0,tau2=0,tau3=0,tau4=0,zeta=parameters[1],
           gamma1=0,gamma2=0,omega1=parameters[2],omega2=parameters[3],
           omega4=parameters[4],d1=0.007165862,d2=parameters[5],M=246538,W=257026)
  
  Ym= 23842278
  Yw =23770470
  G =119810.4
  Im = 58
  Iw = 25
  Ig = 136
  
  N_M = Ym + Im
  N_W = Yw + Im
  N_G = G +  Ig
  
  init<-c(Ym= Ym/N_M,Yw =Yw/N_W, G =G/N_G,
          Im = Im/N_M,Iw = Iw/N_W,Ig = Ig/N_G)
  
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





new.seed = as.integer(runif(1)*2e9) 
cat("Random seed: ", new.seed, "\n")
set.seed(new.seed) 

par = runif(5)


# min.global = 0      # Expected global minimum
# tol = 1e-13         # Tolerance

# 6. Set the search range of each parameter ###########################

dimension=5
# 
lower <- rep(0,dimension) 
upper <- rep(1,dimension) 

# lower <- c(0.03,0.000001,0.000001,0.000001,0.000005) 
# upper <- c(0.07,0.01,0.01,0.1,0.000007) 







ode(init, time, HIV,param)
param<-c(tau1=0,tau2=0,tau3=0,tau4=0,zeta= 0 ,
         gamma1=0,gamma2=0,omega1=0.8332876,omega2=0.9091849,
         omega4=0,d1=0.007165862,d2=0.1422957,M=246538,W=257026)








# 7. Run SA to estimate (optimize) the unknown parameters #############
library(GenSA)
library(deSolve)
#################################################################################



example<-matrix(0, nrow = 1, ncol = 5)
example<-as.data.frame(example)
colnames(example)<-c("par1","par2","par3","par4","par5") #  par1 = 0.0000000 par2 =  0.8332876 par3 = 0.9091849 par4 = 0.0000000 par5 = 0.1422957


out <- GenSA(lower = lower, upper = upper, fn = HIV.model, data = data,
             control=list(max.time = 300))


####SSE4

optimization_view<-function(seed,dim,max.time,lower,upper,func){
  #================================================#
  ## Setting Pakcage state ##
  
  list.of.packages <- c("GenSA","deSolve","ggplot2","reshape2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(new.packages)){
    print("Package does not exist. Install the package.")
    install.packages(new.packages)
    lapply(list.of.packages, require ,  
           character.only = TRUE)
  } else{
    print("Activate package")
    lapply(list.of.packages, require ,  
           character.only = TRUE)
  }
  #================================================#
  ## Setting initial state ##
  set.seed(seed)
  dimension  <-  dim
  # global.min <- global_min
  # tol <- tolearnce
  
  lower_bound <- rep(lower, dimension)
  upper_bound <- rep(upper, dimension)
  out <- GenSA(lower = lower_bound, upper = upper_bound, fn = func,
               control=list(max.time = 300,verbose=TRUE))
  out1<-out[c("value","par","counts")]
  ou1<- as.data.frame(out[c("trace.mat")])
  
  #================================================#
  ## visualizae optimization function state ##
  
  Basic_plot1<-ggplot(data=ou1,aes(x=seq(1:nrow(ou1)),y=trace.mat.temperature))+geom_line()
  plot1<-Basic_plot1+ggtitle("\n Observed temperature in time \n")+ 
    labs(x="Count", y="Temperature")+
    theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 20, color = "Black"))+
    theme(axis.title = element_text(face = "bold", size = 13, color = "Black"))
  
  
  Basic_plot2<-ggplot(data=ou1,aes(x=seq(1:nrow(ou1)),y=trace.mat.current.minimum))+geom_line()
  plot2<-Basic_plot2+ggtitle("\n Observed minimum in time \n")+ 
    labs(x="Count", y="Current minimum")+
    theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 20, color = "Black"))+
    theme(axis.title = element_text(face = "bold", size = 13, color = "Black"))
  
  #================================================#
  print(plot1)
  print(plot2)
  print(out1)
}


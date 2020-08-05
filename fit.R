HIV <- function(year, state_values, parameters) {
  
  
  Ym = state_values [1]        # over 15 years-old men
  Yw = state_values [2]        # over 15 years-old women
  G = state_values [3]      # over 15 years-old MSM
  IUm = state_values [4]       # over 15 years-old infected men
  IDm = state_values [5] 
  IUw = state_values [6]       # over 15 years-old infected women
  IDw = state_values [7] 
  IUG = state_values [8]        # over 15 years-old MSM
  IDG = state_values [9]
  D = state_values [10]        # cumlavtive Death
  
  M = parameters[1]
  W = parameters[2]
  tau1 = parameters[3]
  tau2 = parameters[4]
  tau3 = parameters[5]
  tau4 = parameters[6]
  gamma1 = parameters[7]
  gamma2 = parameters[8]
  omega1 = parameters[9]
  omega2 = parameters[10]
  #omega3 = parameters[11]
  omega4 = parameters[12]
  eta1 = parameters[13]
  eta2 = parameters[14]
  psi= parameters[15]
  d1= parameters[16]
  d2= parameters[17]

  
  with(as.list(c(state_values, parameters)),{
    
  dYm<-  M-(1-tau3)*zeta*Ym+tau3*G-eta1*omega1*(1-gamma1)*(1-tau1)*((IUw+IDw)/(IUw+IDw+Yw))*Ym-d1*Ym
    
  dYw<-  W-eta1*omega2*(1-gamma1)*(1-tau1)*((IUm+IDm)/(IUm+IDm+Ym))*Yw-d1*Yw
    
  dG<-  (1-tau3)*zeta*Ym-tau3*G-eta2*(1-gamma2)*(1-tau2)*(omega4)*((IUG+IDG)/(IUG+IDG+G))*G-d1*G
    
  dIUm<- eta1*omega1*(1-gamma1)*(1-tau1)*((IUw+IDw)/(IUw+IDw+Yw))*Ym*(1-(1-tau4)*psi)-d2*IUm
    
  dIDm<- eta1*omega1*(1-gamma1)*(1-tau1)*(1-tau4)*psi*((IUw+IDw)/(IUw+IDw+Yw))*Ym-d2*IDm                    
    
  dIUw<- eta1*omega2*(1-gamma1)*(1-tau1)*(1-(1-tau4)*psi)*((IUm+IDm)/(IUm+IDm+Ym))*Yw-d2*IUw                                                                 
    
  dIDw<- eta1*omega2*(1-gamma1)*(1-tau1)*(1-tau4)*psi*((IUm+IDm)/(IUm+IDm+Ym))*Yw-d2*IDw                                       
    
  dIUG<- eta2*(1-gamma2)*(1-tau2)*(omega4)*(1-(1-tau4)*psi)*((IUG+IDG)/(IUG+IDG+G))*G-d2*IUG                
    
  dIDG<- eta2*(1-gamma2)*(1-tau2)*(omega4)*(1-tau4)*psi*((IUG+IDG)/(IUG+IDG+G))*G-d2*IDG
    
  dD <- d1*Ym+d1*Yw+d1*G+d2*IUm+d2*IDm+d2*IUw +d2*IDw+d2*IUG+d2*IDG
    
  return(list(c(dYm, dYw, dG, dIUm, dIDm, dIUw, dIDw, dIUG, dIDG, dD)))
  })
  }


setwd("~/Desktop")
library(readr)
pop <- read_csv("model.csv")
pop<-as.data.frame(pop)

pop$cum_reporter <-cumsum(pop$reporter)
pop$cum_death <-cumsum(pop$death)

library(ggplot2)

ggplot(data = pop,aes(x = year, y = reporter)) +geom_point(alpha=0.3,color="red")

ggplot(data = pop,aes(x = year, y = cum_reporter)) +geom_point(alpha=0.3,color="red")

ggplot(data = pop,aes(x = year, y = death)) +geom_point(alpha=0.3,color="red")

ggplot(data = pop,aes(x = year, y = cum_death)) +geom_point(alpha=0.3,color="red")

library(reshape2) 

melt_data <- melt(pop[,c(1,5,6)], id.vars = c("year"))

g <- ggplot(melt_data) + 
  geom_line(aes(x = year, y = value, colour = variable), cex = 0.8, show.legend = T)
g

pop1 <- read_csv("pop.csv")
pop1<-as.data.frame(pop1)


ggplot(data = pop1,aes(x = year, y = Ym)) +geom_point(alpha=0.3,color="red")





HIV.model = function(data, parameters){ 
  
  Ym0 <- as.vector(data$Ym)
  Yw0 <- as.vector(data$Yw)
  G0 <- as.vector(data$G)
  IUM0 <- as.vector(data$IUM)
  IDM0 <- as.vector(data$IDM)
  IUW0 <- as.vector(data$IUW)
  IDW0 <- as.vector(data$IDW)
  IUG0 <- as.vector(data$IUG)
  IDG0 <- as.vector(data$IDG)
  D0<-as.vector(data$D)
  
  param<-c(tau1=0,tau2=0,tau3=0,tau4=0,zeta=parameters[1],
               gamma1=0,gamma2=0,omega1=0.00144,omega2=0.00288,
               omega4=0.04986,eta1=parameters[2],
               eta2=parameters[3],psi=0.1667,d1=0.007165862,d2=5.5047E-05,M=246538,W=257026)
  
  
  init<-c(Ym= 16380526,Yw =16576010, G =131045,IUm = 7804,
          IDm = 1951,IUw = 28 ,IDw = 7,IUG = 7207,IDG = 5854, D = 236162)
  
  time<-seq(1992,2018,1)
  
  Ym0_predicted <- as.vector(ode(init, time, HIV,param)[,2])
  Yw0_predicted <- as.vector(ode(init, time, HIV,param)[,3])
  G0_predicted <- as.vector(ode(init, time, HIV,param)[,4])
  IUM0_predicted <- as.vector(ode(init, time, HIV,param)[,5])
  IDM0_predicted <- as.vector(ode(init, time, HIV,param)[,6])
  IUW0_predicted <- as.vector(ode(init, time, HIV,param)[,7])
  IDW0_predicted <- as.vector(ode(init, time, HIV,param)[,8])
  IUG0_predicted <- as.vector(ode(init, time, HIV,param)[,9])
  IDG0_predicted <- as.vector(ode(init, time, HIV,param)[,10])
  D0_predicted<-as.vector(ode(init, time, HIV,param)[,11])
  
  
  
  #error1 <- Ym0 - YM_predicted 
  #error2<- Yw0 - Ym0_predicted 
  #error3<- G0 - G0_predicted 
  error4<- IUM0 - IUM0_predicted 
  error5<- IDM0 - IDM0_predicted 
  error6<- IUW0 - IUW0_predicted 
  error7<- IDW0 - IDW0_predicted 
  error8<- IUG0 - IUG0_predicted 
  error9<- IDG0 - IDG0_predicted 
  #error10<- D0 - D0_predicted 
  
  
  #SSE1 = sum(error1^2)
  #SSE2 = sum(error2^2)
  #SSE3 = sum(error3^2)
  SSE4 = sum(error4^2)
  SSE5 = sum(error5^2)
  SSE6 = sum(error6^2)
  SSE7 = sum(error7^2)
  SSE8 = sum(error8^2)
  SSE9 = sum(error9^2)
  #SSE10 = sum(error10^2)
 
  return(SSE5)
}


#################################################################################

new.seed = as.integer(runif(1)*2e9) 
cat("Random seed: ", new.seed, "\n")
set.seed(new.seed) 

par = runif(3)


min.global = 0      # Expected global minimum
tol = 1e-13         # Tolerance

# 6. Set the search range of each parameter ###########################

dimension=3

lower <- rep(0,dimension) 
upper <- rep(1,dimension) 


# 7. Run SA to estimate (optimize) the unknown parameters #############
library(GenSA)

#################################################################################


example<-matrix(0, nrow = 1, ncol = 4)
example<-as.data.frame(example)
colnames(example)<-c("par1","par2","par3")


out <- GenSA(par = par, lower = lower, upper = upper, fn = HIV.model, data = pop,
                 control=list(threshold.stop=min.global+tol,
                              verbose=FALSE))


####SSE4
K<-as.vector(out$par)
example1<-rbind(example,K)
example1<-example1[2,]
example1
# ####SSE5
# K<-as.vector(out$par)
# example2<-rbind(example,K)
# example2<-example2[2,]
# example2
# ####SSE6
# K<-as.vector(out$par)
# example3<-rbind(example,K)
# example3<-example3[2,]
# example3
# ####SSE7
# K<-as.vector(out$par)
# example4<-rbind(example,K)
# example4<-example4[2,]
# example4
# ####SSE8
# K<-as.vector(out$par)
# example5<-rbind(example,K)
# example5<-example5[2,]
# example5
# ####SSE9
# K<-as.vector(out$par)
# example6<-rbind(example,K)
# example6<-example6[2,]
# example6



#################################################################################

fit <- mle2(HIV.model, 
            start=list(zeta=0.09, 
                       eta1=0.1, 
                       eta2=0.1, 0.1),  
            method="Nelder-Mead",
            control=list(maxit=1E5,trace=0),
            trace=FALSE)






param<-c(tau1=0,tau2=0,tau3=0,tau4=0,zeta=0.009,
         gamma1=0,gamma2=0,omega1=0.00144,omega2=0.00288,
         omega4=0.04986,eta1=,
         eta2=,psi=,d1=0.006911133,d2=5.5047E-05,M=245482,W=282050)


init<-c(Ym= 21565748,Yw =22107461, G =194092,IUm = 3126,
        IDm = 2120,IUw = 799 ,IDw = 615,IUG = 7294,IDG = 4946, D = 264797)

time<-seq(2012,2018,1)


out<-ode(init, time, HIV,param)

pop$D<-cumsum(pop$D)
pop




pop

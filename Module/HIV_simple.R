#Simple Interaction Model
#================================================#
library(deSolve)

HIV_R <- function(paramters,N_M0 =19084999, N_W0 = 18258526,M_0=19078308, W_0=18257974, I_M0=6691, I_W0=552 ){
  
  derivs <-function(time,initial,paramters){
    
    with(as.list(c(paramters,initial)), {
      dN_M <- M + I_M - mu_1 * M - mu_2 * I_M
      dN_W <- W + I_W - mu_1 * W - mu_2 * I_W
      dM   <- N_M * lambda - beta  * M * I_W  - mu_1 * M
      dW   <- N_W * lambda - gamma * W * I_M  - mu_1 * W
      dI_M <-            beta  * M * I_W  - mu_1 * I_M - mu_2 * I_M
      dI_W <-            gamma * W * I_M  - mu_1 * I_W - mu_2 * I_W
      
      return(list(c(dN_M,dN_W,dM,dW,dI_M,dI_W)))
    })
  }
  
  y <- c(N_M = N_M0,N_W = N_W0 ,M = M_0, W = W_0, I_M = I_M0, I_W = I_W0)
  
  times <-seq(2012,2018,1)
  
  out <-ode(y=y, parms=paramters, times = times, func = derivs)
  
  as.data.frame(out)
  
}

parameters<-c(lambda = 0.005, beta = 0.000000024, gamma = 0.0000000004, mu_1 = 0.005825, mu_2 = 0.1642)
out<- HIV_R(parameters)
out




=======
library(GenSA)

#http://27.101.213.4/ageStatMonth.do
#http://kosis.kr/statHtml/statHtml.do?orgId=101&tblId=DT_2KAA206_OECD
#================================================#



M   <-c(19078308,19135522,19183788,19254023,19295973,19220448,19172749)
W   <-c(18257974,18313807,18360250,18429033,18477966,18404820,18362314)

I_M <-c(6691,7436,8248,9010,9788,10489,11236)
I_W <-c(552,593,624,658,690,716,742)

N_M <-M+I_M
N_W <-W+I_W  

times <-seq(2012,2018,1)

HIV_data<-cbind(times,N_M,N_W,M,W,I_M,I_W)
HIV_data<-data.frame(HIV_data)
#================================================#
#================================================#

# library(FME)
# HIVcost <- function (pars) {
#   out <- HIV_R(pars)
#   cost <- modCost(model = out, obs = HIV_data, err = "sd")
#   return(modCost(model = out, obs = HIV_data, err = "sd", cost = cost))
# }

#================================================#

derivs <-function(time,initial,paramters){
  
  with(as.list(c(paramters,initial)), {
<<<<<<< HEAD
    dN_M <- M + I_M - mu_1 * M - mu_2 * I_M
    dN_W <- W + I_W - mu_1 * W - mu_2 * I_W
    dM   <- N_M * lambda - beta  * M * I_W  - mu_1 * M
    dW   <- N_W * lambda - gamma * W * I_M  - mu_1 * W
    dI_M <-            beta  * M * I_W  - mu_1 * I_M - mu_2 * I_M
    dI_W <-            gamma * W * I_M  - mu_1 * I_W - mu_2 * I_W
=======
    dM   <- growth_Rate * (M + I_M) - beta  * M * I_W  - mu_1 * M
    dW   <- growth_Rate * (M + I_W) - gamma * W * I_M  - mu_1 * W
    dI_M <-            beta  * M * I_W  - (mu_1+mu_2) * I_M 
    dI_W <-            gamma * W * I_M  - (mu_1+mu_2) * I_W
>>>>>>> ec04cf1bf57d2c52209eea1a44fbc00ea70d4c9a
    
    return(list(c(dN_M,dN_W,dM,dW,dI_M,dI_W)))
  })
}


HIV_R <- function(data,parameters){
  
  M_0 <- as.vector(data[1,4])
  W_0 <- as.vector(data[1,5])
  I_M0 <- as.vector(data[1,6])
  I_W0 <- as.vector(data[1,7])
  
  y <- c(N_M = data[1,2], N_W = data[1,3],M = data[1,4], W = data[1,5], I_M = data[1,6], I_W = data[1,7])
  
  times <-seq(2012,2018,1)
  
<<<<<<< HEAD
  parameter<-c(lambda = 0.005, beta = parameters[1], gamma = parameters[2], mu_1 = 0.005825, mu_2 = 0.1642)
  
  M_0_predicted <- as.vector(ode(y, times, derivs,parameter)[,4])
  W_0_predicted <- as.vector(ode(y, times, derivs,parameter)[,5])
  I_M0_predicted <- as.vector(ode(y, times, derivs,parameter)[,6])
  I_W0_predicted <- as.vector(ode(y, times, derivs,parameter)[,7])
=======
  parameter<-c(growth_Rate = 0.0065, beta = parameters[1], gamma = parameters[2], mu_1 = 0.005825, mu_2 = 0.2)
  
  M_0_predicted <- as.vector(ode(y, times, derivs,parameter)[,2])
  W_0_predicted <- as.vector(ode(y, times, derivs,parameter)[,3])
  I_M0_predicted <- as.vector(ode(y, times, derivs,parameter)[,4])
  I_W0_predicted <- as.vector(ode(y, times, derivs,parameter)[,5])
  
  RMSE1<-sqrt((sum(I_M0 - I_M0_predicted)^2 )/(ncol(data)))
  RMSE2<-sqrt((sum(I_W0 - I_W0_predicted)^2 )/(ncol(data)))
>>>>>>> ec04cf1bf57d2c52209eea1a44fbc00ea70d4c9a
  
  RMSE <- (RMSE1 + RMSE2) /2
  
  return(RMSE)
}

library(GenSA)

set.seed(1234) 
dimension <- 2

<<<<<<< HEAD
lower <- rep(0, dimension)
upper <- rep(1, dimension)
=======
>>>>>>> ec04cf1bf57d2c52209eea1a44fbc00ea70d4c9a

lower <- c(0, dimension)
upper <- rep(1, dimension)

out <- optim(par = c(0.00000024,0.0000000008))


out <- GenSA(lower = lower, upper = upper, fn = HIV_R,
             control=list(max.time = 10,verbose=TRUE),data = HIV_data )
out
out[c("value","par","counts")]
ou1<- as.data.frame(out[c("trace.mat")])









#================================================#
## visualizae optimization function state ##
library(ggplot2)

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


#================================================#


<<<<<<< HEAD
HIV_R <- function(paramters, paramters,N_M0 =19084999, N_W0 = 18258526,M_0=19078308, W_0=18257974, I_M0=6691, I_W0=552 ){
  
=======


>>>>>>> ec04cf1bf57d2c52209eea1a44fbc00ea70d4c9a
  derivs <-function(time,initial,paramters){
    
    with(as.list(c(paramters,initial)), {
      dM   <- growth_Rate * (M + I_M) - beta  * M * I_W  - mu_1 * M
      dW   <- growth_Rate * (M + I_W) - gamma * W * I_M  - mu_1 * W
      dI_M <-            beta  * M * I_W  - (mu_1+mu_2) * I_M 
      dI_W <-            gamma * W * I_M  - (mu_1+mu_2) * I_W
      
      return(list(c(dM,dW,dI_M,dI_W)))
    })
  }
  
  parameters<-c(growth_Rate = 0.0065, beta = 0.00000026, gamma = 0.0000000008, mu_1 = 0.005825, mu_2 = 0.2)
  
  y <- c(M = 19078308, W = 18257974, I_M = 6691, I_W = 552)
  
  times <-seq(2012,2018,1)
  
  out <-ode(y=y, parms=parameters, times = times, func = derivs)
  
  as.data.frame(out)
<<<<<<< HEAD
  
}

# parameters<-c(lambda_1 = 248958,lambda_2 =235592, beta = 7.580485e-06, gamma = 0.000000e+00, mu_1 = 0.051, mu_2 = 0.12)
parameters<-c(lambda = 0.005, beta = 0.000000024, gamma = 0.0000000004, mu_1 = 0.005825, mu_2 = 0.1642)

out<- HIV_R(parameters)
cumsum(out$I_M)
cumsum(out$I_W)

=======
>>>>>>> ec04cf1bf57d2c52209eea1a44fbc00ea70d4c9a

  
  library(readr)

  data1 =  read.csv("C:/Users/Monokuma/Desktop/same.csv", header = F, sep=",",encoding = "UTF-8") 
  
  data1[1,1]<- c("2012")
  
  names(data1) <- c("year","category","sex","report")
  
  tspan <- sort(unique(data1[["year"]]))
  
  library(dplyr)
  ggplot(data1,aes(x=tspan,y=report,group_by(Category)))



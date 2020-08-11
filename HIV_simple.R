#Simple Interaction Model
#================================================#
library(deSolve)
HIV_R <- function(paramters, M_0=19078308, W_0=18257974, I_M0=6691, I_W0=552 ){
  
  derivs <-function(time,initial,paramters){
    
    with(as.list(c(paramters,initial)), {
      dM   <- lambda_1 - beta  * M * I_W  - mu_1 * M
      dW   <- lambda_2 - gamma * W * I_M  - mu_1 * W
      dI_M <-            beta  * M * I_W  - mu_2 * I_M
      dI_W <-            gamma * W * I_M  - mu_2 * I_W
      
      return(list(c(dM,dW,dI_M,dI_W)))
    })
  }
  
  y <- c(M = M_0, W = W_0, I_M = I_M0, I_W = I_W0)
  
  times <-seq(2012,2018,1)
  
  out <-ode(y=y, parms=paramters, times = times, func = derivs)
  
  as.data.frame(out)
  
}

parameters<-c(lambda_1 = 248958,lambda_2 = 235592, beta = 0.00000024, gamma = 0.00000000056, mu_1 = 0.051, mu_2 = 0.12)
out<- HIV_R(parameters)
out





#http://27.101.213.4/ageStatMonth.do
#http://kosis.kr/statHtml/statHtml.do?orgId=101&tblId=DT_2KAA206_OECD
#================================================#

M   <-c(19078308,19135522,19183788,19254023,19295973,19220448,19172749)
W   <-c(18257974,18313807,18360250,18429033,18477966,18404820,18362314)

I_M <-c(6691,7436,8248,9010,9788,10489,11236)
I_W <-c(552,593,624,658,690,716,742)

times <-seq(2012,2018,1)

HIV_data<-cbind(times,M,W,I_M,I_W)
HIV_data<-as.data.frame(HIV_data)
#================================================#
#================================================#

# library(FME)
# HIVcost <- function (pars) {
#   out <- HIV_R(pars)
#   cost <- modCost(model = out, obs = HIV_data, err = "sd")
#   return(modCost(model = out, obs = HIV_data, err = "sd", cost = cost))
# }

#================================================#
library(GenSA)

derivs <-function(time,initial,paramters){
  
  with(as.list(c(paramters,initial)), {
    dM   <- lambda_1 - beta  * M * I_W  - mu_1 * M
    dW   <- lambda_2 - gamma * W * I_M  - mu_1 * W
    dI_M <-            beta  * M * I_W  - mu_2 * I_M
    dI_W <-            gamma * W * I_M  - mu_2 * I_W
    
    return(list(c(dM,dW,dI_M,dI_W)))
  })
}


HIV_R <- function(data,parameters){
  
  M_0 <- as.vector(data[1,2])
  W_0 <- as.vector(data[1,3])
  I_M0 <- as.vector(data[1,4])
  I_W0 <- as.vector(data[1,5])
  
  y <- c(M = data[1,2], W = data[1,3], I_M = data[1,4], I_W = data[1,5])
  
  times <-seq(2012,2018,1)
  
  parameters<-c(lambda_1 = 248958,lambda_2 =235592, beta = parameters[1], gamma = parameters[2], mu_1 = 0.051, mu_2 = 0.12)
  
  M_0_predicted <- as.vector(ode(y, times, derivs,parameters)[,2])
  W_0_predicted <- as.vector(ode(y, times, derivs,parameters)[,3])
  I_M0_predicted <- as.vector(ode(y, times, derivs,parameters)[,4])
  I_W0_predicted <- as.vector(ode(y, times, derivs,parameters)[,5])
  
  # RMSE<-sqrt((sum(I_M0 - I_M0_predicted)^2 )/(ncol(data)))
  RMSE<-sqrt((sum(I_W0 - I_W0_predicted)^2 )/(ncol(data)))
  
  return(RMSE)
}


set.seed(1234) 
dimension <- 2
# global.min <- 5908.02
# tol <- 1e-13
# lower <- rep(0, dimension)
# upper <- rep(1, dimension)

lower <- rep(0, dimension)
upper <- rep(0.0001, dimension)

out <- GenSA(par<-c(0.00000002,0.00000002), lower = lower, upper = upper, fn = HIV_R,
             control=list(max.call=10^3,verbose=TRUE),data = HIV_data )
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
## paratmer 2
#  max.call  10^3,     RMSE = 6460.595 // 2.000000e-08 1.934009e-08
#  max.call  10^4,     RMSE = 2524.749 // 2.556424e-08 1.000000e-08
#  max.call 4 * 10^4,  RMSE = 2524.749 // 2.556424e-08 1.000000e-08
#  max.call 4 * 10^4,  RMSE = 66800.46 // 2.556424e-08 1.000000e-08

#================================================#
#  max.call 4 * 10^5,  RMSE = 3822.926 //  7.580485e-06 0.000000e+00


HIV_R <- function(paramters, M_0=19078308, W_0=18257974, I_M0=6691, I_W0=552 ){
  
  derivs <-function(time,initial,paramters){
    
    with(as.list(c(paramters,initial)), {
      dM   <- lambda_1 - beta  * M * I_W  - mu_1 * M
      dW   <- lambda_2 - gamma * W * I_M  - mu_1 * W
      dI_M <-            beta  * M * I_W  - mu_2 * I_M
      dI_W <-            gamma * W * I_M  - mu_2 * I_W
      
      return(list(c(dM,dW,dI_M,dI_W)))
    })
  }
  
  y <- c(M = M_0, W = W_0, I_M = I_M0, I_W = I_W0)
  
  times <-seq(2012,2018,1)
  
  out <-ode(y=y, parms=paramters, times = times, func = derivs)
  
  as.data.frame(out)
  
}

# parameters<-c(lambda_1 = 248958,lambda_2 =235592, beta = 7.580485e-06, gamma = 0.000000e+00, mu_1 = 0.051, mu_2 = 0.12)
parameters<-c(lambda_1 = 248958,lambda_2 =235592, beta = 7.747783e-09, gamma = 0.000000e+00, mu_1 = 0.051, mu_2 = 0.12)

out<- HIV_R(parameters)
cumsum(out$I_M)
cumsum(out$I_W)






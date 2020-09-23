library(deSolve)
library(GenSA)


M   <-c(19078308,19135522,19183788,19254023,19295973,19220448,19172749)
W   <-c(18257974,18313807,18360250,18429033,18477966,18404820,18362314)

I_M <-c(6691,7436,8248,9010,9788,10489,11236)
I_W <-c(552,593,624,658,690,716,742)

times <-seq(2012,2018,1)

HIV_data<-cbind(times,M,W,I_M,I_W)
HIV_data<-data.frame(HIV_data)




derivs <-function(time,initial,paramters){
  
  with(as.list(c(paramters,initial)), {
    
    dM   <- growth_Rate * (M + I_M) - beta  * M * I_W / (W + I_W)  - mu_1 * M
    dW   <- growth_Rate * (W + I_W) - gamma * W * I_M / (M + I_M)  - mu_1 * W
    dI_M <-               beta  * M * I_W / (W + I_W)  - (mu_1 + mu_2) * I_M 
    dI_W <-               gamma * W * I_M / (M + I_M)  - (mu_1 + mu_2) * I_W
    
    return(list(c(dM,dW,dI_M,dI_W)))
  })
}





HIV_R <- function(data,parameters){
  
  M_0 <- as.vector(HIV_data[1,2])
  W_0 <- as.vector(HIV_data[1,3])
  I_M0 <- as.vector(HIV_data[1,4])
  I_W0 <- as.vector(HIV_data[1,5])
  
  y <- c(M = M_0, W = W_0, I_M = I_M0, I_W = I_W0)
 
  parameter<-c(growth_Rate = 0.0065, beta = parameters[1], gamma = parameters[2], mu_1 = 0.005825, mu_2 = 0.2)
  
  times <-seq(2013,2018,1)
  
  M_0_predicted <- as.vector(ode(y, times, derivs,parameter)[-1,2] )
  W_0_predicted <- as.vector(ode(y, times, derivs,parameter)[-1,3] )
  I_M0_predicted <- as.vector(ode(y, times, derivs,parameter)[-1,4] )
  I_W0_predicted <- as.vector(ode(y, times, derivs,parameter)[-1,5] )
  
  
  RMSE1<-sqrt((sum(I_M0 - I_M0_predicted)^2 )/(ncol(data)))
  RMSE2<-sqrt((sum(I_W0 - I_W0_predicted)^2 )/(ncol(data)))
  RMSE <- (RMSE1 + RMSE2) /2
  
  return(RMSE)
}


set.seed(1234) 
dimension <- 2

lower <- rep(0, dimension)
upper <- rep(0.1, dimension)


out <- GenSA(lower = lower, upper = upper, fn = HIV_R,
             control=list(max.call = 10^5, verbose=TRUE) , data = HIV_data )
out
out[c("value","par","counts")]
ou1<- as.data.frame(out[c("trace.mat")])




M_0 <- as.vector(HIV_data[1,2])
W_0 <- as.vector(HIV_data[1,3])
I_M0 <- as.vector(HIV_data[1,4])
I_W0 <- as.vector(HIV_data[1,5])

y <- c(M = M_0, W = W_0, I_M = I_M0, I_W = I_W0)
times <-seq(2013,2018,1)
parameter<-c(growth_Rate = 0.0065, beta = 1.0000000, gamma = 0.1014039, mu_1 = 0.005825, mu_2 = 0.2)

ode(y, times, derivs,parameter)

out



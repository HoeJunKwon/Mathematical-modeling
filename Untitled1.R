library(deSolve)
library(GenSA)


M   <-c(19078308,19135522,19183788,19254023,19295973,19220448,19172749)
W   <-c(18257974,18313807,18360250,18429033,18477966,18404820,18362314)

I_M <-c(6691,7436,8248,9010,9788,10489,11236)
I_W <-c(552,593,624,658,690,716,742)

times <-seq(2012,2018,1)

HIV_data<-cbind(times,M,W,I_M,I_W)
HIV_data<-data.frame(HIV_data)


Report <-c(107,104,125,129,186,219,327,397,533,610,680,749,740,797,768,773,888,868,1013,1081,1018,1060,1008,989)
Years <- seq(1995,2018,1)
# N <- 50000000 # population of the South Korea

N_M <- 22357352
N_W <- 22196358

N <- N_M + N_W

M_report <-c(88,93,107,111,160,194,292,363,502,557,640,687,698,743,710,723,827,808,946,1016,974,1000,958,945)
W_report <-c(19,11,18,18,26,25,35,34,31,53,40,62,42,54,58,50,61,60,67,65,44,60,50,44)






# old <- par(mfrow = c(1, 2))
# plot(Years, Report, type ="b")
# plot(Years, Report, log = "y")
# abline(lm(log10(Report) ~ Years))
# title("Confirmed infections HIV/AIDS in the South Korea", outer = TRUE, line = -2)


SIR <-function(time,initial,paramters){
  
  with(as.list(c(paramters,initial)), {

    dM   <- - alpha * beta  * M * I_W / N   - 0.005825 * M + zeta * I_M
    dW   <- - delta * gamma * W * I_M / N  - 0.005825 * W + yota * I_W 
    dI_M <-                      alpha * beta  *  M  * I_W / N  - 0.005825 * I_M - 0.1642 * I_M -  zeta * I_M
    dI_W <-                      delta * gamma *  W  * I_M / N     - 0.005825 * I_W - 0.1642 * I_W -  yota * I_W
    
    return(list(c(dM,dW,dI_M,dI_W)))
  })
}



library(deSolve)

init <- c(M = N_M - M_report[1], I_M = M_report[1], W = N_W - W_report[1], I_W = W_report[1])

RSS <- function(parameters) {
  parameters<-c(0.5, 0.5, 0.5, 0.5,0.5,0.5)
  names(parameters) <- c("alpha","beta","delta", "gamma","zeta","yota")
  
  out <- ode(y = init, times = Years, func = SIR, parms = parameters)
  
  fit1 <- out[ , 3]
  
  fit2 <- out[ , 5]
  
  A <- sum((M_report - fit1) ^ 2)
  B <- sum((W_report - fit2) ^ 2)
  A

}


library(GenSA)

dimension <- 6
lower <- rep(0, dimension)
upper <- rep(1, dimension)

out <- GenSA(lower = lower, upper = upper, fn = RSS,control=list(max.time = 10, verbose=TRUE))



Opt <- optim(c(0.5, 0.5, 0.5, 0.5,0.5,0.5), RSS, method = "L-BFGS-B", lower = c(0,0,0,0,0,0), upper = c(1,1,1,1,1,1)) # optimize with some sensible conditions
Opt$message


Opt_par <- setNames(Opt$par,  c("alpha","beta","delta", "gamma","zeta","yota"))
print(Opt_par)

t <- 1995:2020 # time in days
fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
col <- 1:3 # colour
print(fit) 
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col)
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col, log = "y")

points(Years, Report)
legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 2, col = col, inset = 0.05)
title("Predicted Cases 2019-nCoV UK (worst case)", outer = TRUE, line = -2)








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



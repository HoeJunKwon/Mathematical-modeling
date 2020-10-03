library(deSolve)
library(GenSA)


Years <- seq(1995,2018,1)
# N <- 50000000 # population of the South Korea

N_M <- 22357352
N_W <- 22196358

N <- N_M + N_W

M_report <-c(88,93,107,111,160,194,292,363,502,557,640,687,698,743,710,723,827,808,946,1016,974,1000,958,945)
W_report <-c(19,11,18,18,26,25,35,34,31,53,40,62,42,54,58,50,61,60,67,65,44,60,50,44)



Report <- M_report + W_report

derivs <-function(time,initial,paramters){
  
  with(as.list(c(paramters,initial)), {
    dM   <- growth_Rate * (M + I_M) - beta  * M * I_W  - mu_1 * M + tau * I_M
    dW   <- growth_Rate * (W + I_W) - gamma * W * I_M   - mu_1 * W + tau * I_W
    dI_M <-            beta  * M * I_W    - (mu_1 + mu_2) * I_M - tau * I_M
    dI_W <-            gamma * W * I_M    - (mu_1 + mu_2) * I_W - tau * I_W
    
    return(list(c(dM,dW,dI_M,dI_W)))
  })
}

parameters<-c(growth_Rate = 0.0065, beta = 0.00000026, gamma = 0.0000000008, mu_1 = 0.005825, mu_2 = 0.2)

init <- c(M = N_M - M_report[1], I_M = M_report[1], W = N_W - W_report[1], I_W = W_report[1])

out <-ode(y=init, parms=parameters, times = Years, func = derivs)

as.data.frame(out)


HIV_R <- function(parameters){
  
  I_M0 <- M_report
  I_W0 <- W_report
  
  init <- c(M = N_M - M_report[1], I_M = M_report[1], W = N_W - W_report[1], I_W = W_report[1])
  
  Years <-seq(1995,2018,1)
  
  parameter<-c(growth_Rate = 0.0065, beta = parameters[1], gamma = parameters[2], tau =parameters[3], mu_1 = 0.005825, mu_2 = 0.2)

  I_M0_predicted <- as.vector(ode(init, Years, derivs,parameter)[,3])
  I_W0_predicted <- as.vector(ode(init, Years, derivs,parameter)[,5])
  
  RMSE1<-sqrt(sum(I_M0 - I_M0_predicted)^2/length(Years))
  RMSE2<-sqrt(sum(I_W0 - I_W0_predicted)^2/length(Years))
  RMSE <- (RMSE1 + RMSE2) /2
  
  return(RMSE)
}

dimension <- 3

lower <- rep(0, dimension)
upper <- rep(1, dimension)

out <- GenSA(par=c(0.01,0.01,0.01),lower = lower, upper = upper, fn = HIV_R, control=list(max.call=10^5,verbose=TRUE))



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






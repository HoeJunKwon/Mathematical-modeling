Men <-c(6691,7436,8248,9010,9788,10489,11236) # 2012~2018 생존감염자  15-64
Women <-c(552,593,624,658,690,716,742) # 2012~2018 생존감염자  15-64
 
Men_Report <-c(776,912,974,937,957,913,896) # 2012~2018 Men Report  15-64
Women_Report <-c(56,59,56,42,53,45,42) # 2012~2018 생존감염자  15-64


Death_rate_Men <-c(0.00162,0.00253,0.00322,0.00426,0.00598,0.00944, 0.01599,0.02464,0.03518,0.05140,0.07844) # 2012    15-64       
Death_rate_Women <-c(0.00090,0.00130,0.00183,0.00257,0.00301,0.00446,0.00621,0.00862,0.01176,0.01807,0.03101) # 2012   15-64


Death_rate_Men <-c(0.00150,0.00237,0.00307,0.00419,0.00577,0.00938,0.01523,0.02350,0.03399,0.04918,0.07376) # 2013     15-64
Death_rate_Women <-c(0.00079,0.00122,0.00176,0.00244,0.00315,0.00433,0.00616,0.00861,0.01155,0.01798,0.02965) # 2013    15-64


Death_rate_Men <-c(0.00158,0.00218,0.00294,0.00399,0.00560,0.00872,0.01440,0.02287,0.03246,0.04687,0.07037) # 2014
Death_rate_Women <-c(0.00093,0.00111,0.00159,0.00234,0.00315,0.00411,0.00591,0.00825,0.01138,0.01675,0.02785) # 2014

Death_rate_Men <-c(0.00130,0.00214,0.00298,0.00365,0.00504,0.00826,0.01384,0.02156,0.03097,0.04578,0.06852)  # 2015
Death_rate_Women <-c(0.00075,0.00112,0.00168,0.00241,0.00300,0.00413,0.00580,0.00766,0.01078,0.01627,0.02677) # 2015

Death_rate_Men <-c(0.00139,0.00204,0.00284,0.00358,0.00488,0.00808,0.01298,0.02046,0.02996,0.04402,0.06647) # 2016
Death_rate_Women <-c(0.00076,0.00117,0.00149,0.00220,0.00297,0.00398,0.00550,0.00780,0.01097,0.01611,0.02605) # 2016

Death_rate_Men <-c(0.00131,0.00203,0.00280,0.00340,0.00479,0.00741,0.01236,0.01959,0.02921,0.04177,0.06355) # 2017
Death_rate_Women <-c(0.00071,0.00102,0.00141,0.00206,0.00280,0.00386,0.00542,0.00726,0.00975,0.01460,0.02458) # 2017

Death_rate_Men <-c(0.00125,0.00199,0.00264,0.00355,0.00499,0.00743,0.01223,0.01934,0.02852,0.04138,0.06352) # 2018
Death_rate_Women <-c(0.00084,0.00129,0.00147,0.00216,0.00292,0.00387,0.00538,0.00755,0.01009,0.01457,0.02400) # 2018


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








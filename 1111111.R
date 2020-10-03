
N_M <- 22357352
N_W <- 22196358

N <- N_M + N_W

M_report <-c(88,93,107,111,160,194,292,363,502,557,640,687,698,743,710,723,827,808,946,1016,974,1000,958,945)
W_report <-c(19,11,18,18,26,25,35,34,31,53,40,62,42,54,58,50,61,60,67,65,44,60,50,44)



Report <- M_report + W_report




Infected <- Report

Day <- 1:(length(Infected))
# N <- 66000000 # pupulation of the UK

old <- par(mfrow = c(1, 2))
plot(Day, Infected, type ="b")
plot(Day, Infected, log = "y")
abline(lm(log10(Infected) ~ Day))
title("Annual Reported infections HIV/AIDS in the South Korea", outer = TRUE, line = -2)

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {

    dS <- -beta/N * I * S + gamma * I
    dI <- beta/N * I * S - gamma * I
    list(c(dS, dI))
  })
}



# SIR <- function(time, state, parameters) {
#   par <- as.list(c(state, parameters))
#   with(par, {
#     
#     dM   <- -beta1/N_W * I_W * M - mu1*M
#     dI_M <-  beta1/N_W * I_W * M  - mu1*I_M - mu2*I_M
#     
#     dW   <- -beta2/N_M * I_M * W - mu1*W
#     dI_W <-  beta2/N_M * I_M * W - mu1*I_W - mu2*I_W
#     
#     list(c(dM,dI_M,dW,dI_W))
#   })
# }








library(deSolve)
# init <- c(M = N_M-M_report[1], I_M = M_report[1],W = N_W-W_report[1], I_W = W_report[1])
init <- c(S = N-Infected[1], I = Infected[1])

RSS <- function(parameters) {
  # names(parameters) <- c("beta1", "beta2","mu1","mu2")
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  # fit2 <- out[ , 5]
  
  # sum((sum((M_report - fit1)^2) + sum((W_report - fit2)^2))/2)
  sum((Infected - fit)^2)
}

library(GenSA)

dimension <- 2
lower <- rep(0, dimension)
upper <- rep(1, dimension)

out <- GenSA(par = c(0.5,0.5), lower = lower, upper = upper, fn = RSS,control=list(max.time = 10, verbose=TRUE))




Opt <- optim(rep(0.5,2), RSS, method = "L-BFGS-B", lower = rep(0,2), upper = rep(1,2)) # optimize with some sensible conditions
Opt$message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

# Opt_par <- setNames(Opt$par, c("beta1", "beta2","mu1","mu2"))
Opt_par <- setNames(Opt$par, c("beta", "gamma"))
print(Opt_par)

t <- 1995:2018 # time in days
fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
col <- 1:3 # colour
fit$cumsum_M <-cumsum(fit$I_M)
fit$cumsum_W <-cumsum(fit$I_W)
print(fit) 


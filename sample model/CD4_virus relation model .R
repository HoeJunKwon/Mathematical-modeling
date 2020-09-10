HIV_R <- function(pars, V_0 = 50000, dV_0 = -200750, T_0=100){
  
  derivs <-function(time,y,pars){
    with(as.list(c(pars,y)), {
      dT <- lambda - rho*T-beta*T*V
      dI <- beta*T*V- delta*I
      dV <- n*delta*I-c*V - beta*T*V
      return(list(c(dT,dI,dV),logV=log(V)))
    })
  }
  
  I_0 <- with(as.list(pars), (dV_0 +c*V_0) / (n*delta))
  y <- c(T=T_0, I = I_0, V=V_0)
  
  times <-c(seq(0,0.8,0.1),seq(2,60,2))
  out <-ode(y=y, parms=pars, times = times, func = derivs)
  as.data.frame(out)
}

HIV <- function (pars, V_0 = 50000, dV_0 = -200750, T_0 = 100) {
  I_0 <- with(as.list(pars), (dV_0 + c * V_0) / (n * delta))
  y <- c(T = T_0, I = I_0, V = V_0)
  times <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, seq(2, 60, by = 2))
 out <- ode(y = y, parms = pars, times = times, func = "derivshiv",initfunc = "inithiv", nout = 1, 
            outnames = "logV", dllname = "FME")
 as.data.frame(out) 
 }


pars <- c(beta = 0.00002, rho = 0.15, delta = 0.55, c = 5.5,lambda = 80, n = 900)

out<- HIV(pars)
# par(mfrow = c(1, 2))
# plot(out$time, out$logV, main = "Viral load", ylab = "log(V)",xlab = "time", type = "b")
# plot(out$time, out$T, main = "CD4+ T", ylab = "-", xlab = "time", type = "b")
# par(mfrow = c(1, 1))

DataLogV <- cbind(time = out$time,logV = out$logV + rnorm(sd = 0.45, n = length(out$logV)), sd = 0.45)

ii <- which (out$time %in% seq(0, 56, by = 4))
DataT <- cbind(time = out$time[ii],T = out$T[ii] + rnorm(sd = 4.5, n = length(ii)), sd = 4.5)

HIVcost <- function (pars) {
  out <- HIV(pars)
  cost <- modCost(model = out, obs = DataLogV, err = "sd")
  return(modCost(model = out, obs = DataT, err = "sd", cost = cost))
}

HIVcost(pars)$model
HIVcost(pars)$residuals

ggplot(HIVcost(pars)$residuals,aes(x=x,y=res,group=name))+geom_point(alpha = 1, aes(color=name))

Sfun <- sensFun(HIVcost, pars)
summary(Sfun)
plot(Sfun, which = c("logV", "T"), lwd = 2)
pairs(Sfun, which = c("logV", "T"), col = c("blue", "green"))






HIVcost2 <- function(lpars) HIVcost(c(exp(lpars), n = 900))
Pars <- pars[1:5] * 2
Fit <- modFit(f = HIVcost2, p = log(Pars))
exp(coef(Fit))
deviance(Fit)

ini <- HIV(pars = c(Pars, n = 900))
final <- HIV(pars = c(exp(coef(Fit)), n = 900))


par(mfrow = c(1, 2))
plot(DataLogV, xlab = "time", ylab = "logV", ylim = c(7, 11))
lines(ini$time, ini$logV, lty = 2)
lines(final$time, final$logV)
legend("topright", c("data", "initial", "fitted"),lty = c(NA, 2, 1), pch = c(1, NA, NA))
plot(DataT, xlab = "time", ylab = "T")
lines(ini$time, ini$T, lty = 2)
lines(final$time, final$T)
par(mfrow = c(1, 1))




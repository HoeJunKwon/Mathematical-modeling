library(deSolve)
library(GenSA)
set.seed(1234) # The user can use any seed.



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
  D = state_values [7]        # Death
  
  with(as.list(c(state_values, parameters)), {
    
    dYm<-M-(1-tau3)*zeta*Ym+tau3*G-eta1*omega1*(1-gamma1)*(1-tau1)*(IUw+IDw)-d1*Ym
    
    dYw<- W-eta1*omega2*(1-gamma1)*(1-tau1)*(IUm+IDm)-d1*Yw
    
    dG<-(1-tau3)*zeta*Ym-tau3*G-eta2*(1-gamma2)*(1-tau2)*(omega3+omega4)*(IUG+IDG)-d1*G
    
    dIUm<- eta1*omega1*(1-gamma1)*(1-tau1)*(IUw+IDw)*(1-tau4*psi)-d2*IUm
    
    dIDm<- eta1*omega1*(1-gamma1)*(1-tau1)*tau4*psi*(IUw+IDw)*-d2*IDm                    
    
    dIUw<- eta1*omega2*(1-gamma1)*(1-tau1)*(IUm+IDm)*(1-tau4*psi)-d2*IUw                                                                 
    
    dIDw<- eta1*omega2*(1-gamma1)*(1-tau1)*tau4*psi*(IUm+IDm)-d2*IDw                                       
    
    dIUG<- eta2*(1-gamma2)*(1-tau2)*(omega3+omega4)*(IUG+IDG)*(1-tau4*psi)-d2*IUG                
    
    dIDG<- eta2*(1-gamma2)*(1-tau2)*(omega3+omega4)*tau4*psi*(IUG+IDG)-d2*IDG
    
    dD <- d1*Ym+d1*Yw+d1*G+d2*IUm+d2*IDm+d2*IUw +d2*IDw+d2*IUG+d2*IDG
    
    return(list(c(dYm, dYw, dG, dIUm, dIDm, dIUw, dIDw, dIUG, dIDG, dD)))
  })
}


parameters_values <- c(
  tau1   =0.1,
  tau2   =0.1,
  tau3   =0.1,
  tau4   =0.1,
  zeta   = 0.009,
  gamma1 = 0.1,
  gamma2 = 0.44,
  omega1 = 0.1,
  omega2 = 0.1,
  omega3 = 0.1,
  omega4 = 0.1,
  eta1   = 0.1,
  eta2   = 0.1,
  psi    = 0.1,
  d1     = 0.1,
  d2     = 0.1,
  M      = 0.1, 
  W      =0.1 
)

initial_values <- c(
  Ym  =  16380526 ,  
  Yw  =  16576010,
  G   =  147425,
  IUm =  140 ,
  IDm =  14,
  IUw =  40 ,    
  IDw =  4,
  IUG =  600,       
  IDG =  60,
  D   =  147425
)


year <- seq(1992,2018,1)

HIV_fit <- ode(
  y = initial_values,
  times = year,
  func = HIV,
  parms = parameters_values 
)

HIV_fit



global.min <- 0
tol <- 1e-13
lower <- c(tau1=0.1,tau2=0.1,tau3=0.1,tau4=0.1,zeta=0.1,
             gamma1=0.1,gamma2=0.1,omega1=0.1,omega2=0.1,omega3=0.1,omega4=0.1,
             eta1=0.1,eta2=0.1,psi=0.1,d1=0.1,d2=0.1,M=162640,W=184854)

upper <- c(tau1=0.9,tau2=0.9,tau3=0.9,tau4=0.9,zeta=0.9,
             gamma1=0.9,gamma2=0.9,omega1=0.9,omega2=0.9,omega3=0.9,omega4=0.9,
             eta1=0.9,eta2=0.9,psi=0.9,d1=0.9,d2=0.9,M=372379,W=564107)

out <- GenSA(lower = lower, upper = upper, fn = HIV,control=list(threshold.stop=global.min+tol,verbose=TRUE))
out[c("value","par")]


set.seed(1234) # The user can use any seed.
dimension <- 18
global.min <- 0
tol <- 1e-13
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
out <- GenSA(lower = lower, upper = upper, fn = HIV,
             control=list(threshold.stop=global.min+tol,verbose=TRUE))
out[c("value","par","counts")]

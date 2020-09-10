
library(deSolve)

HIV.model<-function(year, state_values, parameters){
  with(as.list(c(state_values, parameters)),{
    
    dS   = rho_A * N - rho_B * (1-omega_C) * S - beta_A * gamma_A * (1-omega_A) * (1-omega_D) * S * eta_A - beta_B * gamma_A * (1-omega_A) * (1-omega_D) * S * eta_A - mu_A * S
    
    dG   = rho_B * (1-omega_C) * S - beta_C * gamma_B * (1-omega_A) * (1-omega_D) * (1-omega_B) * G*beta_B*gamma_A * (1-omega_A) * (1-omega_D) * S* eta_B - beta_D * gamma_B * (1-omega_A) * (1-omega_D) * (1-omega_B) * G * eta_B - mu_A * G
    
    dU_S = beta_A * gamma_A * (1-omega_A) * (1-omega_D) * S*eta_A - tau_A*U_S -sigma_B * U_S - mu_A * U_S
    
    dD_S = tau_A * U_S + beta_B*gamma_A * (1-omega_A) * (1-omega_D)*S*eta_A - sigma_A * D_S - mu_A * D_S
    
    dU_G = beta_C * gamma_B * (1-omega_A) * (1-omega_D) * (1-omega_B)*G*eta_B - tau_B * U_G - sigma_D * U_G - mu_A * U_G
    
    dD_G = beta_D * gamma_B * (1-omega_A) * (1-omega_D) * (1-omega_B)*G*eta_B - tau_B * D_G - sigma_C * D_G - mu_A * D_G
    
    dA_S = sigma_A*D_S + sigma_B*U_S - mu_B*A_S
    
    dA_G = sigma_C*D_S + sigma_D*U_S - mu_B*A_G
    
    dD   = mu_A*(S+G+U_S+D_S+U_G+D_G)+mu_B*(A_S+A_G)

    return(list(c(dS, dG, dU_S, dD_S, dU_G, dD_G, dA_S,dA_G,dD)))
  })
}

paramaters_list<-c(rho_A = 0.01 ,rho_B =0.005 , omega_A =0.3 ,omega_B = 0.15 ,omega_C = 0.1 ,
                  omega_D = 0.1 ,gamma_A = 0.0003,gamma_B = 0.0003 , beta_A = 0.0003, beta_B = 0.0003, beta_C = 0.0005,beta_D = 0.0007,mu_A = 0.003,mu_B = 0.01,
                  tau_A = 0.07 ,tau_B = 0.12,eta_A = 0.098 ,eta_B = 0.19,sigma_A = 0.0005,sigma_B= 0.0005,sigma_C= 0.0005,sigma_D= 0.0005)

# the initial values of variables:
initial_values<- c( S = 16500000, U_S = 1500 , D_S = 1000, G = 143066, U_G = 9000 ,D_G = 7000,A_S = 100 ,A_G = 200, D=100)
times <-seq(2010,2015,1)
# solving

out <- ode(initial_values, times, HIV.model, paramaters_list)
out

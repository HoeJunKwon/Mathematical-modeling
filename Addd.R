Hlibrary(readr)
data <- read_csv("sample.csv")
data<-as.data.frame(data)
data





HIV.model<-function(year, state_values, parameters){
  with(as.list(c(state_values, parameters)),{
    
    dM   = Sexual_Identity_Education * G - Annual_change_rate * M   - ((Diagnosis_rate * Heterosexual_Men_Sex_Partner_Shift * Insertive_Penile_vaginal_intercourse * D_M * (D_W + U_W)) + 
      ((1 - Diagnosis_rate) * Heterosexual_Men_Sex_Partner_Shift * Insertive_Penile_vaginal_intercourse * U_M * (D_W + U_W)))  - Death_rate * M
      
    dD_M = Additive_Diagnosis_rate * U_M + Treatment_rate * A_M  + (Diagnosis_rate * Heterosexual_Men_Sex_Partner_Shift * Insertive_Penile_vaginal_intercourse * D_M * (D_W + U_W)) - 
      AIDS_Progress_rate * D_M - Death_rate* D_M
    
    dU_M = (1 - Diagnosis_rate) * Heterosexual_Men_Sex_Partner_Shift * Insertive_Penile_vaginal_intercourse * U_M * (D_W + U_W) - Additive_Diagnosis_rate * U_M - AIDS_Progress_rate * U_M -
    Death_rate * U_M
      
    dA_M = Additive_Diagnosis_rate * (D_M + U_M) - Treatment_rate * A_M - (Death_rate + AIDS_Death_rate)* A_M
      
    dW   =  -((Diagnosis_rate * Heterosexual_Women_Sex_Partner_Shift * Receptive_Penile_vaginal_intercourse * D_W * (D_M + U_M)) + 
        ((1 - Diagnosis_rate) * Heterosexual_Women_Sex_Partner_Shift * Receptive_Penile_vaginal_intercourse * U_W * (D_M + U_M)))  - Death_rate * W
            
    dD_W = Additive_Diagnosis_rate * U_W + Treatment_rate * A_W + (Diagnosis_rate * Heterosexual_Women_Sex_Partner_Shift * Receptive_Penile_vaginal_intercourse * D_W * (D_M + U_M)) - 
      AIDS_Progress_rate * D_W - Death_rate * D_W
      
    dU_W = (1 - Diagnosis_rate) * Heterosexual_Women_Sex_Partner_Shift * Receptive_Penile_vaginal_intercourse * U_W * (D_M + U_M) - 
           AIDS_Progress_rate * U_W - Death_rate * U_W
      
    dA_W =  Additive_Diagnosis_rate * (D_W + U_W) - Treatment_rate * A_W - (Death_rate + AIDS_Death_rate)* A_W  
    
    dG   = Annual_change_rate * M - Sexual_Identity_Education * G - ((Diagnosis_rate * Homosexual_Sex_Partner_Shift * Receptive_Anal_intercourse * D_G * (D_G + U_G))  + 
          ((1 - Diagnosis_rate) * Homosexual_Sex_Partner_Shift * Receptive_Anal_intercourse * U_G * (D_G + U_G)))  - Death_rate * G
      
    dD_G = Additive_Diagnosis_rate * U_G + Treatment_rate * A_G  + (Diagnosis_rate * Homosexual_Sex_Partner_Shift * Receptive_Anal_intercourse * D_G * (D_G + U_G)) - 
      AIDS_Progress_rate * D_G - Death_rate * D_G
      
    dU_G = (1 - Diagnosis_rate) * Homosexual_Sex_Partner_Shift * Receptive_Anal_intercourse * U_G * (D_G + U_G) - Additive_Diagnosis_rate * U_G - AIDS_Progress_rate * U_G -
  Death_rate * U_G
      
    dA_G = Additive_Diagnosis_rate * (D_G + U_G) - Treatment_rate * A_G - (Death_rate + AIDS_Death_rate)* A_G
      
    dD   = Death_rate * (M + D_M + U_M + W + D_W + U_W + G + D_G + U_G) +  (Death_rate + AIDS_Death_rate) * (A_M + A_W + A_G)
    
    return(list(c(dM,dD_M,dU_M,dA_M,dW,dD_W,dU_W,dA_W,dG,dD_G,dU_G,dA_G,D)))
  })
}




HIV.model = function(data, parameters){ 
  
  Ym0 <- as.vector(data$total_M - 0.005 * data$total_M )
  Yw0 <- as.vector(data$total_W)
  G0 <- as.vector(0.005 * data$total_M )
  IM0 <- as.vector(data$Men)
  IW0 <- as.vector(data$Women)
  IG0 <- as.vector(data$Gay)
  
  param<-c(tau1=0,tau2=0,tau3=0,tau4=0,zeta=parameters[1],
           gamma1=0,gamma2=0,omega1=parameters[2],omega2=parameters[3],
           omega4=parameters[4],d1=0.007165862,d2=parameters[5],M=246538,W=257026)
  
  Ym= 23842278
  Yw =23770470
  G =119810.4
  Im = 58
  Iw = 25
  Ig = 136
  
  N_M = Ym + Im
  N_W = Yw + Im
  N_G = G +  Ig
  
  init<-c(Ym= Ym/N_M,Yw =Yw/N_W, G =G/N_G,
          Im = Im/N_M,Iw = Iw/N_W,Ig = Ig/N_G)
  
  time<-seq(1992,2018,1)
  
  Ym0_predicted <- as.vector(ode(init, time, HIV,param)[,2])
  Yw0_predicted <- as.vector(ode(init, time, HIV,param)[,3])
  G0_predicted <- as.vector(ode(init, time, HIV,param)[,4])
  IM0_predicted <- as.vector(ode(init, time, HIV,param)[,5])
  IW0_predicted <- as.vector(ode(init, time, HIV,param)[,6])
  IG0_predicted <- as.vector(ode(init, time, HIV,param)[,7])
  
  
  #error1 <- Ym0 - YM_predicted 
  #error2<- Yw0 - Ym0_predicted 
  #error3<- G0 - G0_predicted 
  error4<- IM0 - IM0_predicted 
  error6<- IW0 - IW0_predicted 
  error8<- IG0 - IG0_predicted 
  #error10<- D0 - D0_predicted 
  
  
  #SSE1 = sum(error1^2)
  #SSE2 = sum(error2^2)
  #SSE3 = sum(error3^2)
  SSE4 = sum(error4^2)
  SSE6 = sum(error6^2)
  SSE8 = sum(error8^2)
  #SSE10 = sum(error10^2)
  
  return(SSE4)
}







parameters_list <- c(Growth_Rate, Insertive_Penile_vaginal_intercourse,
                    Receptive_Penile_vaginal_intercourse,Receptive_Anal_Intercourse,
                    Annual_change_rate,Condom_Usage_Education,UAIC_Prevention_Education,
                    Sexual_Identity_Education, Heterosexual_PrEP,
                    Homosexual_PrEP, Heterosexual_Men_Sex_Partner_Shift,Heterosexual_Women_Sex_Partner_Shift,
                    Homosexual_Sex_Partner_Shift, Death_rate, AIDS_Death_rate,
                    Diagnosis_rate,Additive_Diagnosis_rate, Treatment_rate,AIDS_Progress_rate)

initial_values <-c(M,D_M,A_M,W,D_W,U_W,G,D_G,U_G,D)


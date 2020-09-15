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






parameters_list <- c(Growth_Rate, Insertive_Penile_vaginal_intercourse,
                    Receptive_Penile_vaginal_intercourse,Receptive_Anal_Intercourse,
                    Annual_change_rate,Condom_Usage_Education,UAIC_Prevention_Education,
                    Sexual_Identity_Education, Heterosexual_PrEP,
                    Homosexual_PrEP, Heterosexual_Men_Sex_Partner_Shift,Heterosexual_Women_Sex_Partner_Shift
                    Homosexual_Sex_Partner_Shift, Death_rate, AIDS_Death_rate,
                    Diagnosis_rate,Additive_Diagnosis_rate, Treatment_rate,AIDS_Progress_rate)

initial_values <-c(M,D_M,A_M,W,D_W,U_W,G,D_G,U_G,D)


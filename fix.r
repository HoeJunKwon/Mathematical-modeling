# Poulation --> MSM

N_T = N_M + N_W + N_G 

N_M = S_M + UI_M + DI_M + A_M
N_W = S_W + UI_W + DI_W + A_W

N_G = S_G + UI_G + DI_G + A_G


# Model without control


sirc_model = function (current_timepoint, state_values, parameters){
    
    # create state variables (local variables)
    S_M = state_values [1] # susceptibles Man
    S_W = state_values [2] # susceptibles Woman
    S_MSM = state_values [3] # susceptibles Gay

    I_M = state_values [4] # susceptibles Man
    I_W = state_values [5] # susceptibles Woman
    I_MSM = state_values [6] # susceptibles Gay

    
    UI_M = state_values [7] # Undiagnosised Men
    UI_W = state_values [8] # Undiagnosised Women
    UI_MSM = state_values [9] # Undiagnosised Gay
    
    DI_M = state_values [10] # Diagnosised Men
    DI_W = state_values [11] # Diagnosised Women
    DI_MSM = state_values [12] # Diagnosised Gay
    
    A_M = state_values [13] # AIDS Men
    A_W = state_values [14] # AIDS Women
    A_MSM = state_values [15] # AIDS Gay
    
    with (as.list (parameters), # variable names within parameters can be used
    {
    # compute derivatives
    dS_M = B1 - eta * S_M - c1 * beta1 * I_W / N_W * S_M - mu * S_M
    dS_W = B2 - c2 * beta2 * I_M / N_M * S_W - mu * S_W
    dS_MSM = eta * S_M - c3 * beta3 * I_MSM / N_MSM * S_MSM - mu * S_MSM

    dI_M = c1 * beta1 * I_W / N_W * S_M
    dI_W = c2 * beta2 * I_M / N_M * S_W
    dI_MSM = c3 * beta3 * I_MSM / N_MSM * S_MSM


    dUI_M = (1-v) * c1 * beta1 * I_W / N_W * S_M - p * UI_M - mu * UI_M
    dUI_W = (1-v) * c2 * beta2 * I_M / N_M * S_W - p * UI_W - mu * UI_W
    dUI_MSM = (1-v) * c3 * beta3 * I_G / N_MSM * S_MSM - p * UI_MSM - mu * UI_MSM
    
    dDI_M = v * c1 * beta1 * I_W / N_W * S_M - p * DI_M - mu * DI_M   
    dDI_W = v * c2 * beta2 * I_M / N_M * S_W - p * DI_W - mu * DI_W
    dDI_MSM = v * c3 * beta3 * I_MSM / N_MSM * S_MSM - p * DI_MSM - mu * DI_MSM

    dA_M = p * (UI_M + DI_M) - mu * A_M - d * A_M
    dA_W = p * (UI_M + DI_M) - mu * A_W - d * A_w
    dA_MSM = p * (UI_M + DI_M) - mu * A_MSM - d * A_MSM

    # combine results
    results = c (dS_M,dS_W,dS_MSM,dI_M,dI_W,dI_MSM,UI_M,UI_W,UI_MSM,DI_M,DI_W,DI_MSM,A_M,A_W,A_MSM)
    list (results)
    })
}

B1 # annual flow population for Men     r*N(1-N/K)   r == 0.01 N --> total Men, K -> Density   515  
B2 # annual flow population for Women     r*N(1-N/K)   r == 0.01 N --> total Men, K -> Density   515  

beta1 0.00114
beta2 0.00228
beta3 0.04986
eta 0.005 

N_M + N_MSM -> 19078308

N_W -> 18257974

DI_M + DI_MSM + A_M + A_MSM-> 6691

DI_W+ A_W -> 552



parameters <- c(beta1 = 0.00114, beta2 = 0.00228, beta3 = 0.04986, eta = 0.005,c1 = 2,c2 = 2,c3 = 6)
initial <- c(S_M, S_W, S_MSM, I_M, I_W, I_MSM, UI_M, UI_W, UI_MSM, DI_M,DI_W,DI_MSM,A_M,A_W,A_MSM)
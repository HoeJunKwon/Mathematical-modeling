# Mathematical model for HIV spreads control program with ART treatment
HIV_spreads_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  X1 = state_values [1]        # susceptibles human compartment
  X2 = state_values [2]        # infected human of acute level
  X3 = state_values [3]        # infected human of chronic level
  X4 = state_values [4]        # infected human who receive ART treatment intervention
  X5 = state_values [5]        # infected human who receive ART treatment intervention but with a low level of awareness
  X6 = state_values [6]        # infected human who have failed in the ART tretment precedure
  X7 = state_values [7]        # infected human with HIV/AIDS complication

  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      dX1 = A - (Beta_a * X1 * X2 / N_X) - (Beta_c * X1 * (Theta_h * X4 + Theta_l * X5 + X3 + X6) / N_X) - (Beta_s * X1 * X7 / N_X) - Mu * X1
      dX2 = (Beta_a * X1 * X2 / N_X) + (Beta_c * X1 * (Theta_h * X4 + Theta_l * X5 + X3 + X6) / N_X) + (Beta_s * X1 * X7 / N_X) - Gamma_a * X2 - Mu * X2
      dX3 = Gamma_a * X2 + (1 - R_c) * Gamma_b * X6 - U * X3 - Gamma_c * X3 - Mu * X3 
      dX4 = -Mu * X4 + U * X3 - X4 * Zeta4 + X5 *Zeta5
      dX5 = -Mu * X5 - Row_c * X5 + X4 * Zeta4 - X5 *Zeta5
      dX6 = Row_c * X5 - (1 - R_c) * Gamma_b * X6 - R_c * Gamma_b * X6 - Mu * X6
      dX7 = R_c * Gamma_b * X6 - Delta * X7 - Mu * X7 + Gamma_c * X3
      
      # combine results
      results = c (dX1,dX2,dX3,dX4,dX5,dX6,dX7)
      list (results)
    }
  )
}

parameter_list = c (A = 1000/(65*365), Beta_a = 0.000025, Beta_c = 0.008 * 0.000025, Beta_s = 1.5 * 0.000025, Gamma_a = 1/(5*365), Gamma_c = 1/(3*365), Delta = 0, Mu = 1/(65*365),
                    N_X = 1000, Zeta4 = 0.75 / 365, Zeta5 = 0.25/364, Theta_h = 1, Theta_l = 2, R_c = 0, Row_c = 0.1/365, Gamma_b = 0.5/365)











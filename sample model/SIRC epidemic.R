remove (list = objects() )
library (deSolve)

sirc_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1] # susceptibles
  I = state_values [2] # infectious
  R = state_values [3] # recovered
  C = state_values [4] # carriers
  
  with (
    as.list (parameters), # variable names within parameters can be used
    {
      # compute derivatives
      dS = (-beta * S * I) - (epsilon * beta * S * C)
      dI = (beta * S * I) + (epsilon * beta * S * C) - (gamma * I)
      dR = (1 - rho) * (gamma * I) + (tau * C)
      dC = (rho * gamma * I) - (tau * C)
      
      # combine results
      results = c (dS, dI, dR, dC)
      list (results)
    }
  )
}

contact_rate = 6  # number of contacts per day
transmission_probability = 0.07 # transmission probability
infectious_period = 12  # infectious period
reduced_transmission_rate = 0.09  # chronic carriers compared to acute infections
acute_infections_proportion = 0.67  # acute infections that become carriers
length_of_time_in_carrier_state = 5 # length of time in carrier state

beta_value = contact_rate * transmission_probability
gamma_value = 1 / infectious_period
epsilon_value = 0.09
tau_value = 1 / length_of_time_in_carrier_state
rho_value = 0.67

Ro = beta_value / gamma_value

parameter_list = c (beta = beta_value, gamma = gamma_value, epsilon = epsilon_value, tau = tau_value, rho = rho_value)

V = 4 # carrier hosts
X = 17280 # susceptible hosts
Y = 21 # infectious hosts
Z = 770  # recovered hosts

N = V + X + Y + Z

initial_values = c (S = X/N, I = Y/N, C = V/N, R = Z/N)

timepoints = seq (0, 50, by=1)

output = ode(initial_values, timepoints, sirc_model, parameter_list)


plot (S ~ time, data = output, type='b', col = 'blue')

plot (I ~ time, data = output, type='b', col = 'red')

plot (C ~ time, data = output, type='b', col = 'purple')

plot (R ~ time, data = output, type='b', col = 'green')

# susceptible hosts over time
plot (S ~ time, data = output, type='b', ylim = c(0,1), col = 'blue', ylab = 'S, I, R, C', main = 'SIR/C epidemic')

# remain on same frame
par (new = TRUE)

# infectious hosts over time
plot (I ~ time, data = output, type='b', ylim = c(0,1), col = 'red', ylab = '', axes = FALSE)

# remain on same frame
par (new = TRUE)

# carrier hosts over time
plot (C ~ time, data = output, type='b', ylim = c(0,1), col = 'purple', ylab = '', axes = FALSE)

# remain on same frame
par (new = TRUE)

# recovered hosts over time
plot (R ~ time, data = output, type='b', ylim = c(0,1), col = 'green', ylab = '', axes = FALSE)

library (deSolve)


tiv_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  T = state_values [1]        # target cell
  I = state_values [2]        # infected cell
  V = state_values [3]        # virus 
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      dT = lambda - (d * T) - (k *V * T)
      dI = (k *V * T) - (delta * I) 
      dV = (p * I) - (c * V)
      
      # combine results
      results = c (dT, dI, dV)
      list (results)
    }
  )
}




lambda = 10 #reproduction rate of uninfected cells per cubic millimeter per day

d = 0.01 #death rate of unifected cells per day

k = 0.0000457 #HIV infection rate per cubic millimeter per copies per day

delta = 0.4 #death rate of infected cells per day

p = 38 #virus production rate per copies per cell per day

c = 2.4 # virus clearence rate per day


parameter_list = c (lambda, d, k, delta, p, c)


N = X + Y + Z


initial_values = c (T = X, I = Y, V = Z)



timepoints = seq (0, 900, by=1)


output = lsoda (initial_values, timepoints, tiv_model, parameter_list)


# target cells over time
plot (T ~ time, data = output, type='b', ylim = c(0,2000), col = 'blue', ylab = 'T, I, V', main = 'Viral and Immune System Dynamics of HIV Using TIV Model') 
text("T", col = 'blue')



# remain on same frame
par (new = TRUE)    

# infected cells over time
plot (I ~ time, data = output, type='b', ylim = c(0,2000), col = 'red', ylab = '', axes = FALSE) 
text("I", col = 'red')


# remain on same frame
par (new = TRUE)  

# infectious virus over time
plot (V ~ time, data = output, type='b', ylim = c(0,2000), col = 'purple', ylab = '', axes = FALSE)
text("V", col = 'purple')

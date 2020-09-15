


# Display each parameter distribution
library(ggplot2)

# Aids Report Men
ggplot(data = data, mapping = aes(x = year, y = M_re)) + geom_point()+geom_hline(yintercept=1000)+scale_x_continuous(breaks = seq(1992,2018,1))

# cumulative Aids Report Men
ggplot(data = data, mapping = aes(x = year, y = cumsum(M_re))) + geom_point()+geom_hline(yintercept=10000)+scale_x_continuous(breaks = seq(1992,2018,1))

# Aids Report Women
ggplot(data = data, mapping = aes(x = year, y = W_re)) + geom_point()+scale_x_continuous(breaks = seq(1992,2018,1))

# cumulative Aids Report Women
ggplot(data = data, mapping = aes(x = year, y = cumsum(W_re))) + geom_point()+geom_hline(yintercept=1000)+scale_x_continuous(breaks = seq(1992,2018,1))

# Solving differential equations in R
remove (list = objects() ) 
library (deSolve)

#Step 1: writing the differential equations in R
si_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  I = state_values [2]        # infectious
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      dS = (-beta * S * I)
      dI = ( beta * S * I) - (gamma * I)
      
      # combine results
      results = c (dS, dI)
      list (results)
    }
  )
}

#Step 2: defining some values for the parameters

#prevention(Condom usage,UAIC prevention,Sexual identity,AIDS Recognition and Diagnosis,PreP usage)

# tau1 :Condom usage 
# tau2 :UAIC prevention 
# tau3 :Sexual identity 
# tau4 :AIDS Recognition and Diagnosis -> need more data & fix #
# tau5 :PreP usage 

# zeta1 :AIDS prevention probability by condom

# gamma1 :Heterosexual PrEP
# gamma2 :Homosexual PrEP

# omega1 :Insertive Penile-vaginal intercourse
# omega2 :Receptive penile - vaginal intercourse
# omega3 :Insertive Anal intercouse
# omega4 :Receptive Anal Intercourse

# tau1*zeta1 :actual Condom prevention rate
# tau5*gamma1 :actual Heterosexual PReP prevention rate
# tau5*gamma2 :actual Homosexual PReP prevention rate

# (1-tau2)*omega3 :actual actual infection rate due to UAIC by top
# (1-tau2)*omega4 :actual actual infection rate due to UAIC by bottom

#d1 natural death rate
#d2 HIV/AIDS death rate

contact_rate = 2                     # number of contacts per month
contact_rate = 6                     # number of contacts per month
transmission_probability = 0.04       # transmission probability
transmission_probability = 0.04       # transmission probability

infectious_period = 800                 # infectious period


beta_value = contact_rate * transmission_probability
gamma_value = 1 / infectious_period


Ro = beta_value / gamma_value


parameter_list = c (beta = beta_value, gamma = gamma_value)


X = 130000        # susceptible hosts
Y = 8000          # infectious hosts

N = X + Y


# susceptible hosts over time
plot (S ~ time, data = output, type='b', ylim = c(0,1), col = 'blue', ylab = 'S, I', main = 'SI epidemic') 

# remain on same frame
par (new = TRUE)    

# infectious hosts over time
plot (I ~ time, data = output, type='b', ylim = c(0,1), col = 'red', ylab = '', axes = FALSE) 

par (new = TRUE)  


#Step 3: defining initial values for the variables
initial_values = c (S = X/N, I = Y/N)


#Step 4: the points in time where to calculate variables values

timepoints = seq (1996, 2018, by=1)


output = ode (initial_values, timepoints, si_model, parameter_list)











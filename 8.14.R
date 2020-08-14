setwd("~/Desktop")
library(readr)
pop <- read_csv("data.csv")
data<-as.data.frame(pop)
data$name<-as.factor(data$name)



parameters <- c(beta1 = parameters[1], beta2 = parameters[2], mu = 0.0053, d = 0.11, B1 = 190412 ,B2 = 182225)

initial <- c(S_M = 19071617, S_W = 18257422,I_M = 6691 , I_W = 552 )

time<-seq(2012,2018,1)

ode(initial,time,base_model,parameters)






base_model = function (t, y, parameters){

          dN_M = parameters["growth_rate"] * y ["N_M"] * (1-1/515) - 0.053 * y ["S_M"] - 0.11 * y ["I_M"] - 0.053 * y ["I_M"]
          dN_W = parameters["growth_rate"] * y ["N_W"] * (1-1/515) - 0.053 * y ["S_W"] - 0.11 * y ["I_W"] - 0.053 * y ["I_W"]
            
          dS_M = parameters["growth_rate"] * y ["N_M"] * (1-1/515) - parameters["contact_rate1"] * 0.00114 * y ["I_W"] / (y ["I_W"] + y ["S_W"]) * y ["S_M"] - 0.053 * y ["S_M"]
          dS_W = parameters["growth_rate"] * y ["N_W"] * (1-1/515) - parameters["contact_rate2"] * 0.00228 * y ["I_M"] / (y ["I_M"] + y ["S_M"]) * y ["S_W"] - 0.053 * y ["S_W"]
          
          dI_M = parameters["contact_rate1"] * 0.00114 * y ["I_W"] / (y ["I_W"] + y ["S_W"]) * y ["S_M"] - 0.11 * y ["I_M"] - 0.053 * y ["I_M"] 
          dI_W = parameters["contact_rate2"] * 0.00228 * y ["I_M"] / (y ["I_M"] + y ["S_M"]) * y ["S_W"] - 0.11 * y ["I_W"] - 0.053 * y ["I_W"]
          
          # combine results
          results = c (dN_M ,dN_W, dS_M, dS_W, dI_M, dI_W)
          list (results)
}

fn<- function(x,names_x = c("growth_rate","contact_rate1","contact_rate2"),
              names_y = c("N_M", "N_W", "S_M", "S_W", "I_M", "I_W"),
              tspan = if (!is.null(dat.fitting)) sort(unique(dat.fitting [["year"]])) else NULL,
              base_model, dat.fitting = NULL) {
  
  names(x) <- names_x
  
  y0 <- c(19078308, 19072169, 19071617, 18257422, 6691, 552 )
  
  names(y0) <- names_y
  
  parameters<- x [c("growth_rate","contact_rate1","contact_rate2")]
  
  stopifnot(!is.null(tspan))
  out<- ode(y0, tspan, base_model, parameters)
  
  if (is.null(dat.fitting)){
    rss<- out
  }else {
    MSE<- sum(sapply(c("S_M", "S_W", "I_M", "I_W"), function(nm) {
      o.time<- as.character(dat.fitting [dat.fitting$name == nm,"year"])
      sum((out [, nm ][match(o.time, out [, "time"])] - dat.fitting [dat.fitting$name == nm, "population"])^2)
    }))/nrow(dat.fitting)
  }
  return(MSE)
}

res<- GenSA(fn = fn, lower = c(0,1,1),
            upper = c(1,6,6),
            control = list(max.time = 300), base_model = base_model,
            dat.fitting = data)

################################################################################
################################################################################


base_model = function (t, y, parameters){
  
  dN_M = parameters[1] * y [1] * (1-1/515) - 0.053 * y [3] - 0.11 * y [5] - 0.053 * y [5]
  dN_W = parameters[1] * y [2] * (1-1/515) - 0.053 * y [4] - 0.11 * y [6] - 0.053 * y [6]
  
  dS_M = parameters[1] * y [1] * (1-1/515) - parameters[2] * parameters[4] * y [6] / (y [6] + y [4]) * y [3] - 0.053 * y [3]
  dS_W = parameters[1] * y [2] * (1-1/515) - parameters[3] * parameters[5] * y [6] / (y [5] + y [3]) * y [4] - 0.053 * y [4]
  
  dI_M = parameters[2] * parameters[4] * y [6] / (y [6] + y [4]) * y [3] - 0.11 * y [5] - 0.053 * y [5] 
  dI_W = parameters[3] * parameters[5] * y [5] / (y [5] + y [3]) * y [4] - 0.11 * y [6] - 0.053 * y [6]
  
  # combine results
  results = c (dN_M ,dN_W, dS_M, dS_W, dI_M, dI_W)
  list (results)
}

tspan = if (!is.null(data)) sort(unique(data [["year"]])) else NULL
y0 <- c(19078308, 19072169, 19071617, 18257422, 6691, 552 )
parameters <- c( 0.0543583869, 1.3389311305, 3.9340771551, 0.0560642089,0.7106316541)
out<- ode(y0, tspan, base_model, parameters)
out
data

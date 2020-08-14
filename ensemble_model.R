#set work space

setwd("~/Desktop")
library(readr)
pop <- read_csv("data.csv")
data<-as.data.frame(pop)


HIV_model <-function(){
  
}



ensemble_model<-function(seed,dim,global_min,tolearnce,lower,upper,func){
  #================================================#
  ## Setting Pakcage state ##
  
  list.of.packages <- c("GenSA","deSolve","ggplot2","reshape2","gridExtra")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(new.packages)){
    print("Package does not exist. Install the package.")
    install.packages(new.packages)
    lapply(list.of.packages, require ,  
           character.only = TRUE)
  } else{
    print("Activate package")
    lapply(list.of.packages, require ,  
           character.only = TRUE)
  }
  
  funct<- function(x,names_x = c("", "", "","","","","","",""),
                names_y = c("Ym", "Yw","G","I_M","I_W","I_G"),initial = c(), 
                names_parameters = c("", "",""),tspan = if (!is.null(dat.fitting)) sort(unique(dat.fitting [["year"]])) else NULL,
                func = HIV_model, dat.fitting = NULL) {
    
    names(x) <- names_x
    y0 <- initial   # initial condition
    names(y0) <- names_y
    parameters<- x [c("", "", "","","","","","","")] # parameter condition
    
    stopifnot(!is.null(tspan))
    
    out<- ode(y0, tspan, func, parameters)
    
    if (is.null(dat.fitting)) { # dat.fitting  = data
      MSE<- out
    }else {
      rss<- sum(sapply(c("Wirte each variable name "), function(nm) { # Add when required variables are present
        
        o.time<- as.character(dat.fitting [dat.fitting$year == nm,"year"]) #  Data required must be year For time series data that does not include year, reference date when data is measured is required
        
        sum((out [, nm ][match(o.time, out [, "year"])] - dat.fitting [dat.fitting$year == nm, "population"])^2, na.rm = TRUE) # When the data shape is changed to melt, the value of value must be expressed.
      }))/(2*nrow(dat.fitting))
    }
    return(MSE)
  }
  
  #================================================#
  ## Setting initial state ##
  set.seed(seed)
  dimension  <-  dim
  global.min <- global_min
  tol <- tolearnce
  
  lower_bound <- rep(lower, dimension)
  upper_bound <- rep(upper, dimension)
  out <- GenSA(lower = lower_bound, upper = upper_bound, fn = funct,
               control=list(threshold.stop=global.min+tol,verbose=TRUE))     # or Maxit call =c(10^3,10^4,10^5,...)
  out1<-out[c("value","par","counts")]
  ou1<- as.data.frame(out[c("trace.mat")]) 
  
  #================================================#
  ## visualizae optimization function state ##
  
  Basic_plot1<-ggplot(data=ou1,aes(x=seq(1:nrow(ou1)),y=trace.mat.temperature))+geom_line()
  plot1<-Basic_plot1+ggtitle("\n Observed temperature in time \n")+ 
    labs(x="Count", y="Temperature")+
    theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 20, color = "Black"))+
    theme(axis.title = element_text(face = "bold", size = 13, color = "Black"))
  
  
  Basic_plot2<-ggplot(data=ou1,aes(x=seq(1:nrow(ou1)),y=trace.mat.current.minimum))+geom_line()
  plot2<-Basic_plot2+ggtitle("\n Observed minimum in time \n")+ 
    labs(x="Count", y="Current minimum")+
    theme(plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 20, color = "Black"))+
    theme(axis.title = element_text(face = "bold", size = 13, color = "Black"))
  
  #================================================#
  print(plot1)
  print(plot2)
  print(out1)
  
}


ensemble_model(seed=1234,dim=30,global_min=0,tol=1e-13,lower=,upper=,func=HIV_model) # get optimization parameter


#================================================================================================#
# inset optimization parameter

HIV_fixed_model <-function(){
  
}


tspan = if (!is.null(dat.fitting)) sort(unique(dat.fitting [["year"]])) else NULL
y0<-c()
parameters<-c("" = , "" = , "" = ,"" = ,"" = ,"" = ,"" = ,"" = ,"" = ,)

out <- ode(y0, tspan, HIV_fixed_model, parameters)







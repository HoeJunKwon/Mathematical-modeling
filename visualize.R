# NOT RUN {

# Try Rastrgin function (The objective function value for global minimum
# is 0 with all components of par are 0.)
Rastrigin <- function(x) {
  sum(x^2 - 10 * cos(2 * pi  * x)) + 10 * length(x)
}
#
set.seed(1234) 
dimension <- 30
global.min <- 0
tol <- 1e-13
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
out <- GenSA(lower = lower, upper = upper, fn = Rastrigin,
             control=list(threshold.stop=global.min+tol,verbose=TRUE))
out[c("value","par","counts")]



optimization_view<-function(seed,dim,global_min,tolearnce,lower,upper,func){
  #================================================#
  ## Setting initial state ##
  
  list.of.packages <- c("GenSA","globalOptTests","deSolve","ggolot2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(new.packages)){
    install.packages(new.packages)
    lapply(list.of.packages, require ,  
           character.only = TRUE)
  } else{
    print("Package does not exist. Install the package.")
    lapply(list.of.packages, require ,  
           character.only = TRUE)
  }
  
  #================================================#
  ## Setting initial state ##
  set.seed(seed)
  dimension  <-  dim
  global.min <- global_min
  tol <- tolearnce
  
  lower_bound <- rep(lower, dimension)
  upper_bound <- rep(upper, dimension)
  out <- GenSA(lower = lower_bound, upper = upper_bound, fn = func,
               control=list(threshold.stop=global.min+tol,verbose=TRUE))
  out[c("value","par","counts")]
  ou1<- as.data.frame(out[c("trace.mat")])
  
  #================================================#
  ## visualizae optimization function state ##
  library(ggplot2)
  
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
}


optimization_view(seed=1234,dim=5,global_min=0,tol=1e-13,lower=-5.12,upper=5.12,func=Rastrigin)



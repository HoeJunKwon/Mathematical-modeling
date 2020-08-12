
library(reshape2)
#================================================#
## cost function state ##
#================================================#
# MSE #

if (is.null(dat.fitting)) { # dat.fitting  = data
  MSE<- out
}else {
  rss<- sum(sapply(c("Wirte each variable name "), function(nm) { # Add when required variables are present
    
    o.time<- as.character(dat.fitting [dat.fitting$year == nm,"year"]) #  Data required must be year For time series data that does not include year, reference date when data is measured is required
    
    sum((out [, nm ][match(o.time, out [, "year"])] - dat.fitting [dat.fitting$year == nm, "population"])^2, na.rm = TRUE) # When the data shape is changed to melt, the value of value must be expressed.
  }))/(2*nrow(dat.fitting))
}
return(MSE)


#================================================#
#================================================#


require(mkin)
require(GenSA)
require(deSolve)
options(digits = 9)
set.seed(1234)
str(FOCUS_2006_D)

df<- function(t, y, parameters) {
  d_parent<- -parameters["k_parent_sink"]*y["parent"]-parameters["k_parent_m1"]*y ["parent"]
  d_m1 <- parameters["k_parent_m1"]*y["parent"]-parameters ["k_m1_sink"]* y ["m1"]
  return(list(c(d_parent, d_m1)))
}


fn<- function(x,names_x = c("parent_0", "k_parent_sink", "k_parent_m1","k_m1_sink"),
              names_y = c("parent", "m1"),
              names_parameters = c("k_parent_sink", "k_parent_m1","k_m1_sink"),tspan = if (!is.null(dat.fitting)) sort(unique(dat.fitting [["time"]])) else NULL,
              df, dat.fitting = NULL) {
  m1_0 = 0
  names(x) <- names_x
  y0 <- c(x ["parent_0"], m1_0)
  names(y0) <- names_y
  parameters<- x [c("k_parent_sink", "k_parent_m1", "k_m1_sink")]
  stopifnot(!is.null(tspan))
  out<- ode(y0, tspan, df, parameters)
  
  if (is.null(dat.fitting)) {
    rss<- out
  }else {
    rss<- sum(sapply(c("parent", "m1"), function(nm) {
      o.time<- as.character(dat.fitting [dat.fitting$name == nm,"time"])
      sum((out [, nm ][match(o.time, out [, "time"])] - dat.fitting [dat.fitting$name == nm, "value"])^2,na.rm = TRUE)
    }))
  }
  return(rss)
}

res<- GenSA(fn = fn, lower = c(90, rep(0.001, 3)),
            upper = c(110, rep(0.1, 3)),
            control = list(max.time = 5), df = df,
            dat.fitting = FOCUS_2006_D)

names(res$par) <- c("parent_0", "k_parent_sink", "k_parent_m1","k_m1_sink")

print(round(res$par, digits = 6))

#================================================#
df1<- function(t, y, parameters) {

  d_parent<- 99.5984880305-parameters[1]*y["parent"]-0.050778*y ["parent"]
  d_m1 <- parameters[3]*y["parent"]-parameters [2]* y ["m1"]
  return(list(c(d_parent, d_m1)))
}

#================================================#


# x<-c(99.5984880305,0.0479201577,0.0507776401,0.0052606502)
# m1_0 = 0
# names_x = c("parent_0", "k_parent_sink", "k_parent_m1","k_m1_sink")
# names(x) <- names_x
# y0 <- c(x ["parent_0"], m1_0)
# names_y = c("parent0", "m1")
# names(y0) <- names_y
# parameters<- x [c("k_parent_sink", "k_parent_m1", "k_m1_sink")]
# out<- ode(y0, seq(0,120,1), df1, parameters)

#================================================#
data<-as.data.frame(FOCUS_2006_D)
data1<-data[seq(1,22),]
# data2<-data[seq(23,44),]
# data3<-cbind(data[seq(1,11),1],as.data.frame(melt(out[,2])))
# colnames(data3)<-c("name","value")
# data4<-cbind(as.data.frame(out[,1]),data3)
# colnames(data4)<-c("time","name","value")
# 
# #================================================#

plot<-ggplot(data = data,aes(x=time,y=value))+geom_point(aes(color=name))
plot
# plot3<-ggplot(data = data1,aes(x=time,y=value))+geom_point()
# plot4<-ggplot(data = data2,aes(x=time,y=value))+geom_point()
# 
# plot5<-ggplot(data = as.data.frame(out),aes(x=time,y=parent0))+geom_line()
# plot5
# 
# plot(out)



df<- function(t, y, parameters) {
  d_parent<- -parameters[1]*y[1]-parameters[2] * y [1]
  d_m1 <- parameters[2]*y[1] - parameters [3]* y [2]
  return(list(c(d_parent, d_m1)))
}
tspan = if (!is.null(FOCUS_2006_D)) sort(unique(FOCUS_2006_D [["time"]])) else NULL
y0<-c(99.598491,0)
parameters<-c("k_parent_sink" = 0.047920, "k_parent_m1" = 0.050778, "k_m1_sink" = 0.005261)
out <- ode(y0, tspan, df, parameters)
out[,1]
colnames(out) <- c("time","fit_patent","fit_m1")
out1<-melt(out,id.vars = c("time"))[seq(12,33),seq(2,3,1)]
colnames(out1)<-c("name","value")
out2<-cbind(out1,out[,1])
colnames(out2)<-c("name","value","time")
out2
data
data1<-rbind(out2,data)
data1<-as.data.frame(data1)


library(gridExtra)

cast()


plot<-ggplot(data = data1[seq(23,66),],aes(x=time,y=value))+geom_point(aes(color=name))
plot
plot1<-plot+geom_line(data = data1[seq(1,22),],aes(x=time,y=value,color=name))+geom_point(aes(color=name))
plot1


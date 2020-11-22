library(MASS)
library(plot3D)
book <-read.csv("/Users/victoria/Desktop/Book7.csv")

q1 <-as.numeric( as.character( book$Year[1:20]))
q2 <-as.numeric( as.character( book$Infection[1:20]))
q3 <-as.numeric( as.character( book$Cost[1:20]))

den3d <- kde2d(q2,q3)


# the new part:
library(plotly)
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()


par(mfrow = c(1, 1))
panelfirst <- function(pmat) {
  zmin <- min(-quakes$depth)
  XY <- trans3D(quakes$long, quakes$lat,
                z = rep(zmin, nrow(quakes)), pmat = pmat)
  scatter2D(XY$x, XY$y, colvar = quakes$mag, pch = ".",
            cex = 2, add = TRUE, colkey = FALSE)
  xmin <- min(quakes$long)
  XY <- trans3D(x = rep(xmin, nrow(quakes)), y = quakes$lat,
                z = -quakes$depth, pmat = pmat)
  scatter2D(XY$x, XY$y, colvar = quakes$mag, pch = ".",
            cex = 2, add = TRUE, colkey = FALSE)
}
with(quakes, scatter3D(x = long, y = lat, z = -depth, colvar = mag,
                       pch = 16, cex = 1.5, xlab = "longitude", ylab = "latitude",
                       zlab = "depth, km", clab = c("Richter","Magnitude"),
                       main = "Earthquakes off Fiji", ticktype = "detailed",
                       panel.first = panelfirst, theta = 10, d = 2,
                       colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75))
)



persp3D(z = volcano, zlim = c(-60, 200), phi = 20,
        colkey = list(length = 0.2, width = 0.4, shift = 0.15,
                      cex.axis = 0.8, cex.clab = 0.85), lighting = TRUE, lphi = 90,
        clab = c("","height","m"), bty = "f", plot = TRUE)
# create gradient in x-direction
Vx <- volcano[-1, ] - volcano[-nrow(volcano), ]
# add as image with own color key, at bottom
image3D(z = -60, colvar = Vx/10, add = TRUE,
        colkey = list(length = 0.2, width = 0.4, shift = -0.15,
                      cex.axis = 0.8, cex.clab = 0.85),
        clab = c("","gradient","m/m"), plot = TRUE)
# add contour

contour3D(z = -60+0.01, colvar = Vx/10, add = TRUE,
          col = "black", plot = TRUE)


# library
library(rgl)

# This is to output a rgl plot in a rmarkdown document. Note that you must add webgl=TRUE, results='hide' in the chunck header
#library(knitr)
#knit_hooks$set(webgl = hook_webgl)
iris
# Data: the iris data is provided by R
book <-read.csv("/Users/victoria/Desktop/Book1.csv")

# Add a new column with color
mycolors <- c("aliceblue","aquamarine","azure3","bisque","blue","blueviolet","brown",
              "cadetblue","chartreuse","chartreuse4","chocolate4","darksalmon","gold",
              "darkslateblue","darkslategrey","paleturquoise4")

book$color <- mycolors[as.numeric(as.factor(book$Case))]

book$Case <- as.factor(book$Case)

# Plot
par(mar=c(0,0,0,0))
plot3d( 
  x=book$PrEP, y=book$Edu, z=book$Infected, 
  col = book$color, 
  type = 's', 
  radius = 30,
  xlab="PrEP", ylab="Edu", zlab="Infected")

plot3d( 
  x=book$PrEP, y=book$Edu, z=book$Cost, 
  col = book$color, 
  type = 's', 
  radius = 1000000000,
  xlab="PrEP", ylab="Edu", zlab="Cost")


with(book, plot3d(PrEP,Edu,Infected,type='l',col = as.integer(Case)))


Data1 <-as.matrix(read.csv("/Users/victoria/Desktop/data/Infected.csv",header=FALSE))
Data2 <-as.matrix(read.csv("/Users/victoria/Desktop/data/Cost.csv",header=FALSE))

persp(Data1,theta=-120,phi=20,ltheta = -120,shade = 0.75,xlab = 'Case', ylab='Year',zlab = "Infeted People")
persp(Data2,theta=-120,phi=20,ltheta = -120,shade = 0.75,xlab = 'Case', ylab='Year',zlab = "Total Cost")

persp(Data2,theta=0)


library(readr)
A <- as.matrix(read_csv("Downloads/Case/A.csv",col_names = FALSE))
B <- as.matrix(read_csv("Downloads/Case/B.csv",col_names = FALSE))
C <- as.matrix(read_csv("Downloads/Case/C.csv",col_names = FALSE))
D <- as.matrix(read_csv("Downloads/Case/D.csv",col_names = FALSE))
E <- as.matrix(read_csv("Downloads/Case/E.csv",col_names = FALSE))
F <- as.matrix(read_csv("Downloads/Case/F.csv",col_names = FALSE))
G <- as.matrix(read_csv("Downloads/Case/G.csv",col_names = FALSE))
H <- as.matrix(read_csv("Downloads/Case/H.csv",col_names = FALSE))
I <- as.matrix(read_csv("Downloads/Case/I.csv",col_names = FALSE))
J <- as.matrix(read_csv("Downloads/Case/J.csv",col_names = FALSE))
K <- as.matrix(read_csv("Downloads/Case/K.csv",col_names = FALSE))
L <- as.matrix(read_csv("Downloads/Case/L.csv",col_names = FALSE))
M <- as.matrix(read_csv("Downloads/Case/M.csv",col_names = FALSE))
N <- as.matrix(read_csv("Downloads/Case/N.csv",col_names = FALSE))
O <- as.matrix(read_csv("Downloads/Case/O.csv",col_names = FALSE))
P <- as.matrix(read_csv("Downloads/Case/P.csv",col_names = FALSE))
Q <- as.matrix(read_csv("Downloads/Case/Q.csv",col_names = FALSE))
R <- as.matrix(read_csv("Downloads/Case/R.csv",col_names = FALSE))
S <- as.matrix(read_csv("Downloads/Case/S.csv",col_names = FALSE))
T <- as.matrix(read_csv("Downloads/Case/T.csv",col_names = FALSE))


Alpha <-list(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T)

library(animation)
saveGIF({
    z <- T
    persp(T, theta = 45, phi = 35, expand = 0.4, col = "orange")
  }, movie.name = "20.gif")




# saveGIF({
#   for(i in 1:20){
#     
#     z <- Alpha[i]
#     persp(z, theta = 45, phi = 35, expand = 0.4, col = "orange")
#   }
# }, interval = 0.1, ani.width = 550, ani.height = 550)
# 

















)

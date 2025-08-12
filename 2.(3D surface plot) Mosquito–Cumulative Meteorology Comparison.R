
######################################################## 3D surface plot

# load library
library(plot3D)
library(mgcv)

# User-defined function
fig.3D.f <- function(DB, x.var, y.var){
  
  ##### Preporcessing
  {
  DB$mosq.log <- log(DB$Mosq_abund+1, 10)
  Temp.DB <- DB[c(x.var, y.var, "mosq.log")]
  
  intv.f <- function(T.var){
  if(length(grep("Prec_cum", T.var))>0 ){
    intv <- 10  # 10 mm
  } else {
    intv <- 1 # 1 celsius or 1 day
  }
  return(intv)
  }
  
  x.int <- intv.f(x.var)
  y.int <- intv.f(y.var)
  
  Temp.DB2 <- Temp.DB
  Temp.DB2$var.cut1 <-cut(Temp.DB2[,1], seq(0,(max(Temp.DB2[,1])+x.int), by=x.int), right = F)
  Temp.DB2$var.cut2 <-cut(Temp.DB2[,2], seq(0,(max(Temp.DB2[,2])+y.int), by=y.int), right = F)
  
  Temp.DB3 <- aggregate(mosq.log ~ var.cut1+var.cut2, Temp.DB2, mean )
  x <- (as.numeric(Temp.DB3$var.cut1)-1)*x.int
  y <- (as.numeric(Temp.DB3$var.cut2)-1)*y.int
  z <- Temp.DB3$mosq.log
  }
  
  ##### GAM
  {
  fit <- gam(z ~ s(x, y), method="REML")
  grid.lines = 100
  x.pred <- seq(min(x), max(x), length.out = grid.lines)
  y.pred <- seq(min(y), max(y), length.out = grid.lines)
  xy <- expand.grid( x = x.pred, y = y.pred)
  z.pred <- matrix(predict(fit, newdata = xy), 
                   nrow = grid.lines, ncol = grid.lines)
  fitpoints <- predict(fit)
  }

  ##### 3D plot
  par(mar=c(2,2,2,2))
  ttheta = 230; pphi = 20
  scatter3D(x, y, z, 
            colvar = seq(min(z.pred), max(z.pred), length=length(z)),
            bty = "u", #b2
            col.grid = "grey", #darkgrey
            col.panel = "grey30",  #  black
            cex=0, pch=21, type = "p", bg="black",  
            theta = ttheta, #60 #좌(+), 우(-)회전 # 210
            phi = pphi, #위(+), 아래(-)회전
            xlab=x.var, ylab=y.var, zlab="log(Ind.+1,10)",
            #@ xlab="", ylab="", zlab="",
            zlim=c(range(z.pred)[1], range(z.pred)[2]),
            # zlim=c(0, range(z.pred)[2]),
            #@ zlim=c(range(z)[1], range(z)[2]),
            ticktype="detailed", cex.lab=0.8, cex.axis=1, 
            surf= list(x=x.pred, y=y.pred, z=z.pred, facets=NA),# fit=fitpoints,: point까지 거리
            clab = "mosq", # c("Occurrence","probability"),
            colkey = list(side = 4, length = 0.5)
  )
  
}

##### Description of input elements
# DB: Input DB 
# x.var: Variable name to be displayed on the x-axis in the 3D plot
# y.var: Variable name to be displayed on the y-axis in the 3D plot
# => Select one of the following for x.var and y.var: 
# "Temp_min_cum#", "Temp_av_cum#", "Temp_max_cum#", "Prec_#", or "Prec_day_cum#", where # denotes the cumulative number of days.


######################################################## Usage example
# load sample data and library
load("sample data_NH.RData") # DB.NH
library(plot3D)
library(mgcv)

# Run the function
fig.3D.f(DB.NH, x.var = "Temp_min_cum5", y.var = "Prec_cum5")

##### Description of output
# 3D plot with surface plane: (x-axis) x.var, (y-axis) y.var, (z-axis) Mosquito abundance (log10)


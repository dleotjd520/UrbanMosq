
######################################################## 2D plot with polynomial equation

# load library
library(reshape2)
library(mgcv)

# User-defined function
fig.2D.f <- function(DB.NH, DB.WH, T.var){
  ##### Preporcessing
  {
  DB.NH$cluster <- "NH"
  DB.WH$cluster <- "WH"
  total.DB <- rbind(DB.NH, DB.WH)
  total.DB$mosq.log <- log(total.DB$Mosq_abund+1, 10)
  
  Temp.DB <- total.DB[c("cluster", "mosq.log", T.var)]
  Temp.DB2 <- melt(Temp.DB, id.vars=c("cluster", "mosq.log"))
  
  if( length(grep("Prec_cum", T.var))>0 ){
    intv <- 10  # 10 mm
  } else {
    intv <- 1 # 1 celsius or 1 day
  }
  Temp.DB2$variable <- as.character(as.vector(Temp.DB2$variable))
  Temp.DB2$cut <- cut(Temp.DB2$value, seq(0,(max(Temp.DB2$value)+intv),by=intv), right=F)
  
  # Mean & Standard error
  Temp.DB3 <- aggregate(mosq.log ~ cut+variable+cluster, Temp.DB2, function(x){ c(m=mean(x), se=sd(x)/sqrt(length(x))) } )
  }
  
  ##### Division by habitat type 
  T.NH <- Temp.DB3[Temp.DB3$cluster == "NH", ]
  T.WH <- Temp.DB3[Temp.DB3$cluster == "WH", ]
  
  T.DB.NH <-data.frame("y" = T.NH$mosq.log[,"m"], "x" = (as.numeric(T.NH$cut)-1)*intv) # 
  T.DB.WH <-data.frame("y" = T.WH$mosq.log[,"m"], "x" = (as.numeric(T.WH$cut)-1)*intv) # 
  
  
  ##### Polynomial selection 
  {
  poly.f <- function(DB, i){
    if( i == 1 ){
      m <- lm(y ~ x, DB)
    } else{ 
      m <- lm(as.formula(paste0("y ~ x +", paste0("I(x^",2:i,")", collapse = "+"))), DB)
    }#
    return(m)
  }
  
  DB.m.NH2 <- DB.m.WH2 <- c()
  for(i in 1:4){
    # i=3
    NH.m <- poly.f(T.DB.NH, i)
    WH.m <- poly.f(T.DB.WH, i)

    # adj.R2
    NH.adj.r2 <- summary(NH.m)$adj.r.squared
    WH.adj.r2 <- summary(WH.m)$adj.r.squared

    DB.m.NH <- data.frame("i" = i, "adj.r2" =  NH.adj.r2)
    DB.m.WH <- data.frame("i" = i, "adj.r2" =  WH.adj.r2)
    
    DB.m.NH2 <- rbind(DB.m.NH2, DB.m.NH)
    DB.m.WH2 <- rbind(DB.m.WH2, DB.m.WH)
  }
  
  final.NH.i <- DB.m.NH2$i[which.max(DB.m.NH2$adj.r2)]
  final.WH.i <- DB.m.WH2$i[which.max(DB.m.WH2$adj.r2)]
  
  f.NH.m <- poly.f(T.DB.NH, final.NH.i)
  f.WH.m <- poly.f(T.DB.WH, final.WH.i)
  
  # model prediction
  NH.m.f <- predict(f.NH.m, se.fit = TRUE)
  WH.m.f <- predict(f.WH.m, se.fit = TRUE)
  
  NH.upr <- NH.m.f$fit + (1.96 * NH.m.f$se.fit) #95%
  NH.lwr <- NH.m.f$fit - (1.96 * NH.m.f$se.fit) #95%
  WH.upr <- WH.m.f$fit + (1.96 * WH.m.f$se.fit) #95%
  WH.lwr <- WH.m.f$fit - (1.96 * WH.m.f$se.fit) #95%
  }
  
  ##### Plot condition
  {
  y.range <- range(na.omit(c(Temp.DB3$mosq.log[,1], 
                             Temp.DB3$mosq.log[,1]+Temp.DB3$mosq.log[,2], 
                             Temp.DB3$mosq.log[,1]-Temp.DB3$mosq.log[,2]))
  )
  
  x.range <- na.omit(c(T.DB.NH$x, T.DB.WH$x))
  
  Ordinal.exp <- function(x) {
    suffix <- ifelse(x %% 100 %in% 11:13, "th",
                     ifelse(x %% 10 == 1, "st",
                            ifelse(x %% 10 == 2, "nd",
                                   ifelse(x %% 10 == 3, "rd", "th"))))
    paste0(x, suffix)
  }
  
  if( length(grep("Prec_cum", T.var))>0 ){
    xlabel <-"Precipitation (mm)"
  } else if( length(grep("Prec_day", T.var))>0 ){
    xlabel <-"Days"
  } else {
    xlabel <-"Temperature (celsius)"
  }
  
  par(mar=c(4,4,3,1))
  plot(T.DB.NH$x, T.DB.NH$y, las=1, type="n",
       main = paste0(T.var, " - NH:",Ordinal.exp(final.NH.i),", WH:", Ordinal.exp(final.WH.i)),
       xlab=xlabel, ylab="log10(mosq+1)",
       xlim=c(min(x.range), max(x.range)), 
       ylim=c(min(y.range), max(y.range))
  )
  
  # polygon 
  NH.i.for<-order(T.DB.NH$x)
  NH.i.back<-order(T.DB.NH$x, decreasing = T)
  NH.x.polygon <- c(T.DB.NH$x[NH.i.for],T.DB.NH$x[NH.i.back])
  NH.y.polygon <- c(NH.upr[NH.i.for], NH.lwr[NH.i.back])
  NH.y.polygon[NH.y.polygon<0] <- 0 

  WH.i.for<-order(T.DB.WH$x)
  WH.i.back<-order(T.DB.WH$x, decreasing = T)
  WH.x.polygon <- c(T.DB.WH$x[WH.i.for], T.DB.WH$x[WH.i.back])
  WH.y.polygon <- c(WH.upr[WH.i.for], WH.lwr[WH.i.back] )
  WH.y.polygon[WH.y.polygon<0] <- 0 

  grid()
  
  polygon(NH.x.polygon, NH.y.polygon, col= "gray90", border=NA)
  polygon(WH.x.polygon, WH.y.polygon, col= "gray90", border=NA)
  
  # points & error bar
  arrows(x0=T.DB.NH$x, 
         y0=T.NH$mosq.log[,"m"]-T.NH$mosq.log[,"se"], 
         y1=T.NH$mosq.log[,"m"]+T.NH$mosq.log[,"se"],  col="#009681", code=3, angle=90, length= 0.01 )
  arrows(x0=T.DB.WH$x, 
         y0=T.WH$mosq.log[,"m"]-T.WH$mosq.log[,"se"], 
         y1=T.WH$mosq.log[,"m"]+T.WH$mosq.log[,"se"], col="#CC476B", code=3, angle=90, length= 0.01 )

  points(T.DB.NH$x, T.NH$mosq.log[,"m"], pch=21, bg="white", col="#009681", cex=1.2)
  points(T.DB.WH$x, T.WH$mosq.log[,"m"], pch=21, bg="white", col="#CC476B", cex=1.2)
  
  # GAM regression curve
  lines(T.DB.NH$x, NH.m.f$fit, col="#009681", lwd=2)
  lines(T.DB.WH$x, WH.m.f$fit, col="#CC476B", lwd=2)
  }
}

##### Description of input elements
# DB.NH: DB for nonwaterside habitats 
# DB.WH: DB for waterside habitats
# T.var: Target variable name for analysis
# => Select one from: "Temp_min_cum#", "Temp_av_cum#", "Temp_max_cum#", "Prec_#", or "Prec_day_cum#", 
# where # denotes the cumulative number of days.


######################################################## Usage example

# load sample data and library
load("sample data_NH.RData") # DB.NH
load("sample data_WH.RData") # DB.WH
library(reshape2)
library(mgcv)

# Run the function
fig.2D.f(DB.NH, DB.WH, T.var = "Prec_cum30") # 

##### Description of output
# 2D plot: Mosquito abundance (log10) - Cumulative meteorology
# line: polynomial equation
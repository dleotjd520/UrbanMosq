
######################################################## Convergent Cross Mapping function
# Note: This function is a customized version of the CCM functionality from the rEDM package, tailored to the needs of this study.
# Reference: https://ha0ye.github.io/rEDM/articles/rEDM.html

# load library
# devtools::install_github("ha0ye/rEDM")
library(rEDM)

# User-defined function
CCM.f <- function(DB, T.var, col.forward = "#fb6330", col.reverse = "#138468"){
  
  ##### Data preprocessing
  {
  t.DB1 <- rbind(
    DB[as.character(DB$Date) %in% as.Date(as.Date("2012-05-01"):as.Date("2012-10-31")), ],
    DB[as.character(DB$Date) %in% as.Date(as.Date("2013-05-01"):as.Date("2013-10-31")), ],
    DB[as.character(DB$Date) %in% as.Date(as.Date("2015-05-01"):as.Date("2015-10-31")), ]
  )
  
  t.DB2 <- t.DB1[!is.na(t.DB1$Mosq_abund), ]
  t.DB3 <- t.DB2
  t.DB3$mosq.log <- log(t.DB3$Mosq_abund+1, 10)
  }

  ##### Preparation for CCM
  {
  vars <- c("mosq.log", T.var)
  composite_ts <- t.DB3[vars]
  composite_ts <- cbind("Date" = t.DB3$Date, composite_ts)
  
  # normalize each time series within a plot
  data_by_plot <- split(composite_ts, t.DB3$Site)
  
  segments_end <- cumsum(sapply(data_by_plot, NROW))
  segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
  segments <- cbind(segments_begin, segments_end)
  
  # Choose random segments for prediction
  set.seed(123)
  rndlib <- sample(1:NROW(segments), floor(NROW(segments) * 0.75))
  composite_lib <- segments[rndlib, ]
  composite_pred <- segments[-rndlib, ]
  
  simplex_out <- lapply(vars, function(var) {
    simplex(composite_ts[, c("Date", var)], E = 1:15, lib = composite_lib, pred = composite_pred)
  })
  
  best_E <- sapply(simplex_out, function(df) {
    df$E[which.max(df$rho)]
  })
  
  lib_sizes = c(100, round(seq(0, nrow(composite_ts), length.out= 15), 0)[-1])
  }
  
  ##### CCM
  {
  set.seed(456)
  var_mosq <- ccm(composite_ts, lib = segments, pred = segments, lib_column = T.var, 
                  target_column = "mosq.log", E = best_E[2], lib_sizes = lib_sizes, 
                  silent = TRUE)
  mosq_var <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "mosq.log",
                  target_column = T.var, E = best_E[1], lib_sizes = lib_sizes, 
                  silent = TRUE)
  var_mosq_means <- ccm_means(var_mosq) # forward
  mosq_var_means <- ccm_means(mosq_var) # reverse
  }
  
  ##### CCM plot
  {
  xlimit<- range(c(mosq_var_means$lib_size, var_mosq_means$lib_size))
  ylimit<- range(c(mosq_var_means$rho, var_mosq_means$rho))
  
  par(mar=c(4,4,3,3))
  plot(1, 1, type="n", xlim = c(xlimit[1], xlimit[2]), ylim = c(ylimit[1], ylimit[2]), las=1, 
       xlab="Library Size", ylab="Cross Map Skill (rho)")
  grid()
  lines(var_mosq_means$lib_size, var_mosq_means$rho, lwd=2, lty=1, col=col.forward)
  lines(mosq_var_means$lib_size, mosq_var_means$rho, lwd=2, lty=2, col=col.reverse)
  
  points(var_mosq_means$lib_size, var_mosq_means$rho, pch=21, cex=1.2, bg=col.forward)
  points(mosq_var_means$lib_size, mosq_var_means$rho, pch=19, cex=1.2, col=col.reverse)
  }
  
  ##### Output: CCM result
  total.DB <- rbind(mosq_var_means, var_mosq_means)
  return(total.DB)
}

##### Description of input elements
# DB: Input DB 
# T.var: Meteorological factor to be tested for its effect on mosquito abundance (forward direction)
# => Select one of them: "Temp_min_cum#", "Temp_av_cum#", "Temp_max_cum#", "Prec_#", or "Prec_day_cum#", where # denotes the cumulative number of days.
# col.forward: Color specification for the forward direction in the CCM plot
# col.reverse: Color specification for the reverse direction in the CCM plot


######################################################## Usage example
# load sample data and library
load("sample data_NH.RData") # DB.NH
load("sample data_WH.RData") # DB.WH
library(rEDM)

# Preparation of example data for CCM
DB.NH$Site <- "NH"
DB.WH$Site <- "WH"
target.DB = rbind(DB.NH, DB.WH)

# Run the function
CCM.DB <- CCM.f(target.DB, T.var = "Prec_cum5")
head(CCM.DB)

##### Description of output
# CCM plot: Plot showing both forward and reverse directions
# CCM.DB: Resulting DB generated through CCM
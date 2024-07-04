rm(list=ls(all=TRUE))

### LOAD PACKAGES ###
library(CoinCalc) #
library(ggplot2) # version 3.4.4
library(ggpubr) # version 0.6.0
library(patchwork) # version 1.2.0
library(grid) # 
library(gridExtra) # version 2.3
library(ggrepel)# version 0.9.4
library(raster) # version 3.6-26
library(latticeExtra) # 0.6-30
library(poweRlaw) #

######

### DIRECTORY SETTING ###
pc.name <- Sys.info()["user"]
dir.data <- paste("C:/Users/",pc.name,"/SomMigTime/data",sep = "")
dir.udf <- paste("C:/Users/",pc.name,"/SomMigTime/udf",sep = "")
dir.out <- paste("C:/Users/",pc.name,"/SomMigTime/out",sep = "")
dir.rdat <- paste("C:/Users/",pc.name,"/SomMigTime/Rdat",sep = "")

######

### DATA CLEANING ###
# Load user-defined functions
setwd(dir.udf)
dump("eca_udf", file="eca_udf");source("eca_udf.R")
dump("bursty_GB_udf", file="bursty_GB_udf");source("bursty_GB_udf.R")
dump("bursty_KJ_udf", file="bursty_KJ_udf");source("bursty_KJ_udf.R")
dump("BM_cal_udf", file="BM_cal_udf");source("BM_cal_udf.R")
###

# Load Data
setwd(dir.data)
distName <- read.csv("dist_name.csv",skip = 0,header = FALSE)
outidp_dat <- read.csv("Som16to20.csv")
conf_dat <- read.csv("alsha_dat.csv",skip = 0, header = FALSE)
rain_dat <- read.csv("no_rain.csv",skip = 0, header = FALSE)
tmax_dat <- read.csv("tmax2.csv",skip = 0, header = FALSE)
#pop <- read.csv("pop.csv",skip = 0, header = FALSE)
shp <- shapefile("Som_Admbnda_Adm2_UNDP_sort2.shp");
colnames(conf_dat) <- distName[[1]]
colnames(rain_dat) <- distName[[1]]
colnames(tmax_dat) <- distName[[1]]
###

# Create migration event time series
timeDistOther <- ts_crt(outidp_dat,0,"mig")
###

#Plot Pop
pop <-c()
for (i in 1:74){
  int <- sum(outidp_dat[outidp_dat$fromDist == i & outidp_dat$fromDist != outidp_dat$toDist,6])
  pop <- c(pop,int)
}

shp$pop <- log(pop)
trellis.par.set(axis.line=list(col=NA))
plot1 <- spplot(shp, "pop",
                col.regions = rev(heat.colors(101)),
                at = seq(min(log(pop))-0.01,max(log(pop))+0.01,length=100),
                lwd = 1
)
setwd(dir.out)
png(filename="pop.png",width = 1000, height = 1000, res = 300)
plot1 #Figure 1A
dev.off()
###
######

### BURSTINESS AND MEMORY ###
### B-M Calculation
BMOther <- BM_cal(timeDistOther)
###

# B-M plane visualization
setwd(dir.out)
p_BM <- BM_plot(BMOther)
p_BM
ggsave("BM_other.png", height = 7 , width = 7) # Figure 2A
###

# Event time series plot and inter-event time distribution
setwd(dir.out)
district <- "Baki"
p_ts1 <- ts_plot(district, distName,timeDistOther) #
p_ts1
ggsave("BhighMlow.png", height = 0.7, width = 3) # Figure 2B top

dd1 <- iet_plot(district, distName,timeDistOther,"power-law")
png(filename="tempplot_BhighMlow.png",width = 1200, height = 900, res = 300)
plot(dd1,xlab = expression(italic(tau)), ylab = expression(italic(P(T >= tau))))
lines(dd1, col = 2)
text(x = 8, y = 0.1, labels = bquote(italic("\u03B1") ~ "=" ~ .(round(dd1$pars,2))), col = "red", pos = 4)
dev.off() # Figure 2B bottom

district <- "Ceel Waaq"
p_ts2 <- ts_plot(district, distName,timeDistOther) #
p_ts2 # Figure 2C top
ggsave("BhighMhigh.png", height = 0.7, width = 3)

dd2 <- iet_plot(district, distName,timeDistOther,"power-law")
png(filename="tempplot_BhighMhigh.png",width = 1200, height = 900, res = 300)
plot(dd2,xlab = expression(italic(tau)), ylab = expression(italic(P(T >= tau))))
lines(dd2, col = 2)
text(x = 8, y = 0.25, labels = bquote(italic("\u03B1") ~ "=" ~ .(round(dd2$pars,2))), col = "red", pos = 4)
dev.off() # Figure 2C bottom

district <- "Ceel Dheer"
p_ts3 <- ts_plot(district, distName,timeDistOther) #
p_ts3 # Figure 2D top
ggsave("Bzero.png", height = 0.7, width = 3)

dd3 <- iet_plot(district, distName,timeDistOther,"exponential")
png(filename="tempplot_Bzero.png",width = 1200, height = 900, res = 300)
plot(dd3,xlab = expression(italic(tau)), ylab = expression(italic(P(T >= tau))))
lines(dd3, col = 2)
text(x = 2, y = 0.1, labels = bquote(italic("\u03B1") ~ "=" ~ .(round(dd3$pars,2))), col = "red", pos = 4)
dev.off() # Figure 2D bottom

district <- "Diinsoor"
p_ts4 <- ts_plot(district, distName,timeDistOther) #
p_ts4 # Figure 2E top
ggsave("Blow.png", height = 0.7, width = 3)

dd4 <- iet_plot(district, distName,timeDistOther,"power-law")
png(filename="tempplot_Blow.png",width = 1200, height = 900, res = 300)
plot(dd4, xlab = expression(italic(tau)), ylab = expression(italic(P(T >= tau))))
dev.off() # Figure 2E bottom
###

### ECA ### To save time, please run pre-run RData file below
# Find lags
lag_out1 <- lag_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,1)
lag_out4 <- lag_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,4)
lag_out13 <- lag_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,13)
###

# Run ECA
eca_pois1 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,1,lag_out1,"poisson") 
eca_pois4 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,4,lag_out4,"poisson") 
eca_pois13 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,13,lag_out13,"poisson")

eca_wtsur1 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,1,lag_out1,"wt.surrogate") # takes long time
eca_wtsur4 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,4,lag_out4,"wt.surrogate") # takes long time
eca_wtsur13 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,13,lag_out13,"wt.surrogate") # takes long time

# takes long time
eca_shuf1 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,1,lag_out1,"shuffle.surrogate") # takes long time
eca_shuf4 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,4,lag_out4,"shuffle.surrogate") # takes long time
eca_shuf13 <- real_eca(timeDistOther,conf_dat,rain_dat,tmax_dat,9,13,lag_out13,"shuffle.surrogate") # takes long time

#save.image(directory), for future use of the data to reduce running times
###

# Load ECA results 
setwd(dir.rdat)
load("dat_manu.RData")

# Clean ECA outputs
conf_idp_stat1 <- eca_clean(eca_pois1[[1]],eca_wtsur1[[1]],eca_shuf1[[1]])
conf_idp_stat4 <- eca_clean(eca_pois4[[1]],eca_wtsur4[[1]],eca_shuf4[[1]])
conf_idp_stat13 <- eca_clean(eca_pois13[[1]],eca_wtsur13[[1]],eca_shuf13[[1]])

rain_idp_stat1 <- eca_clean(eca_pois1[[2]],eca_wtsur1[[2]],eca_shuf1[[2]])
rain_idp_stat4 <- eca_clean(eca_pois4[[2]],eca_wtsur4[[2]],eca_shuf4[[2]])
rain_idp_stat13 <- eca_clean(eca_pois13[[2]],eca_wtsur13[[2]],eca_shuf13[[2]])

rain_conf_idp_stat1 <- eca_clean(eca_pois1[[3]],eca_wtsur1[[3]],eca_shuf1[[3]])
rain_conf_idp_stat4 <- eca_clean(eca_pois4[[3]],eca_wtsur4[[3]],eca_shuf4[[3]])
rain_conf_idp_stat13 <- eca_clean(eca_pois13[[3]],eca_wtsur13[[3]],eca_shuf13[[3]])

tmax_idp_stat1 <- eca_clean(eca_pois1[[4]],eca_wtsur1[[4]],eca_shuf1[[4]])
tmax_idp_stat4 <- eca_clean(eca_pois4[[4]],eca_wtsur4[[4]],eca_shuf4[[4]])
tmax_idp_stat13 <- eca_clean(eca_pois13[[4]],eca_wtsur13[[4]],eca_shuf13[[4]])

tmax_conf_idp_stat1 <- eca_clean(eca_pois1[[5]],eca_wtsur1[[5]],eca_shuf1[[5]])
tmax_conf_idp_stat4 <- eca_clean(eca_pois4[[5]],eca_wtsur4[[5]],eca_shuf4[[5]])
tmax_conf_idp_stat13 <- eca_clean(eca_pois13[[5]],eca_wtsur13[[5]],eca_shuf13[[5]])
###

# Save in shapefile
# Precursor coincidence rate
shp$conf_idp_pre_other1 <- eca_pois1[[1]]$`precursor coincidence rate`
shp$conf_idp_pre_other4 <- eca_pois4[[1]]$`precursor coincidence rate`
shp$conf_idp_pre_other13 <- eca_pois13[[1]]$`precursor coincidence rate`

shp$rain_idp_pre_other1 <- eca_pois1[[2]]$`precursor coincidence rate`
shp$rain_idp_pre_other4 <- eca_pois4[[2]]$`precursor coincidence rate`
shp$rain_idp_pre_other13 <- eca_pois13[[2]]$`precursor coincidence rate`

shp$rain_conf_idp_pre_other1 <- eca_pois1[[3]]$`precursor coincidence rate`
shp$rain_conf_idp_pre_other4 <- eca_pois4[[3]]$`precursor coincidence rate`
shp$rain_conf_idp_pre_other13 <- eca_pois13[[3]]$`precursor coincidence rate`

shp$tmax_idp_pre_other1 <- eca_pois1[[4]]$`precursor coincidence rate`
shp$tmax_idp_pre_other4 <- eca_pois4[[4]]$`precursor coincidence rate`
shp$tmax_idp_pre_other13 <- eca_pois13[[4]]$`precursor coincidence rate`

shp$tmax_conf_idp_pre_other1 <- eca_pois1[[5]]$`precursor coincidence rate`
shp$tmax_conf_idp_pre_other4 <- eca_pois4[[5]]$`precursor coincidence rate`
shp$tmax_conf_idp_pre_other13 <- eca_pois13[[5]]$`precursor coincidence rate`

# coincidence occurrence
shp$conf_idp_pre_n_other1 <- eca_pois1[[1]]$`precursor coincidence rate` * eca_pois1[[1]]$`N precursor`
shp$conf_idp_pre_n_other4 <- eca_pois4[[1]]$`precursor coincidence rate` * eca_pois4[[1]]$`N precursor`
shp$conf_idp_pre_n_other13 <- eca_pois13[[1]]$`precursor coincidence rate` * eca_pois13[[1]]$`N precursor`

shp$rain_idp_pre_n_other1 <- eca_pois1[[2]]$`precursor coincidence rate` * eca_pois1[[2]]$`N precursor`
shp$rain_idp_pre_n_other4 <- eca_pois4[[2]]$`precursor coincidence rate` * eca_pois4[[2]]$`N precursor`
shp$rain_idp_pre_n_other13 <- eca_pois13[[2]]$`precursor coincidence rate` * eca_pois13[[2]]$`N precursor`

shp$rain_conf_idp_pre_n_other1 <- eca_pois1[[3]]$`precursor coincidence rate` * eca_pois1[[3]]$`N precursor`
shp$rain_conf_idp_pre_n_other4 <- eca_pois4[[3]]$`precursor coincidence rate` * eca_pois4[[3]]$`N precursor`
shp$rain_conf_idp_pre_n_other13 <- eca_pois13[[3]]$`precursor coincidence rate` * eca_pois13[[3]]$`N precursor`

shp$tmax_idp_pre_n_other1 <- eca_pois1[[4]]$`precursor coincidence rate` * eca_pois1[[4]]$`N precursor`
shp$tmax_idp_pre_n_other4 <- eca_pois4[[4]]$`precursor coincidence rate` * eca_pois4[[4]]$`N precursor`
shp$tmax_idp_pre_n_other13 <- eca_pois13[[4]]$`precursor coincidence rate` * eca_pois13[[4]]$`N precursor`

shp$tmax_conf_idp_pre_n_other1 <- eca_pois1[[5]]$`precursor coincidence rate` * eca_pois1[[5]]$`N precursor`
shp$tmax_conf_idp_pre_n_other4 <- eca_pois4[[5]]$`precursor coincidence rate` * eca_pois4[[5]]$`N precursor`
shp$tmax_conf_idp_pre_n_other13 <- eca_pois13[[5]]$`precursor coincidence rate` * eca_pois13[[5]]$`N precursor`

# P-value
shp$conf_idp_pre_pv_pois1 <- eca_pois1[[1]]$`p-value precursor`
shp$conf_idp_pre_pv_pois4 <- eca_pois4[[1]]$`p-value precursor`
shp$conf_idp_pre_pv_pois13 <- eca_pois13[[1]]$`p-value precursor`
shp$conf_idp_pre_pv_wtsur1 <- eca_wtsur1[[1]]$`p-value precursor`
shp$conf_idp_pre_pv_wtsur4 <- eca_wtsur4[[1]]$`p-value precursor`
shp$conf_idp_pre_pv_wtsur13 <- eca_wtsur13[[1]]$`p-value precursor`
shp$conf_idp_pre_pv_shuf1 <- eca_shuf1[[1]]$`p-value precursor`
shp$conf_idp_pre_pv_shuf4 <- eca_shuf4[[1]]$`p-value precursor`
shp$conf_idp_pre_pv_shuf13 <- eca_shuf13[[1]]$`p-value precursor`

shp$rain_idp_pre_pv_pois1 <- eca_pois1[[2]]$`p-value precursor`
shp$rain_idp_pre_pv_pois4 <- eca_pois4[[2]]$`p-value precursor`
shp$rain_idp_pre_pv_pois13 <- eca_pois13[[2]]$`p-value precursor`
shp$rain_idp_pre_pv_wtsur1 <- eca_wtsur1[[2]]$`p-value precursor`
shp$rain_idp_pre_pv_wtsur4 <- eca_wtsur4[[2]]$`p-value precursor`
shp$rain_idp_pre_pv_wtsur13 <- eca_wtsur13[[2]]$`p-value precursor`
shp$rain_idp_pre_pv_shuf1 <- eca_shuf1[[2]]$`p-value precursor`
shp$rain_idp_pre_pv_shuf4 <- eca_shuf4[[2]]$`p-value precursor`
shp$rain_idp_pre_pv_shuf13 <- eca_shuf13[[2]]$`p-value precursor`

shp$rain_conf_idp_pre_pv_pois1 <- eca_pois1[[3]]$`p-value precursor`
shp$rain_conf_idp_pre_pv_pois4 <- eca_pois4[[3]]$`p-value precursor`
shp$rain_conf_idp_pre_pv_pois13 <- eca_pois13[[3]]$`p-value precursor`
shp$rain_conf_idp_pre_pv_wtsur1 <- eca_wtsur1[[3]]$`p-value precursor`
shp$rain_conf_idp_pre_pv_wtsur4 <- eca_wtsur4[[3]]$`p-value precursor`
shp$rain_conf_idp_pre_pv_wtsur13 <- eca_wtsur13[[3]]$`p-value precursor`
shp$rain_conf_idp_pre_pv_shuf1 <- eca_shuf1[[3]]$`p-value precursor`
shp$rain_conf_idp_pre_pv_shuf4 <- eca_shuf4[[3]]$`p-value precursor`
shp$rain_conf_idp_pre_pv_shuf13 <- eca_shuf13[[3]]$`p-value precursor`

shp$tmax_idp_pre_pv_pois1 <- eca_pois1[[4]]$`p-value precursor`
shp$tmax_idp_pre_pv_pois4 <- eca_pois4[[4]]$`p-value precursor`
shp$tmax_idp_pre_pv_pois13 <- eca_pois13[[4]]$`p-value precursor`
shp$tmax_idp_pre_pv_wtsur1 <- eca_wtsur1[[4]]$`p-value precursor`
shp$tmax_idp_pre_pv_wtsur4 <- eca_wtsur4[[4]]$`p-value precursor`
shp$tmax_idp_pre_pv_wtsur13 <- eca_wtsur13[[4]]$`p-value precursor`
shp$tmax_idp_pre_pv_shuf1 <- eca_shuf1[[4]]$`p-value precursor`
shp$tmax_idp_pre_pv_shuf4 <- eca_shuf4[[4]]$`p-value precursor`
shp$tmax_idp_pre_pv_shuf13 <- eca_shuf13[[4]]$`p-value precursor`

shp$tmax_conf_idp_pre_pv_pois1 <- eca_pois1[[5]]$`p-value precursor`
shp$tmax_conf_idp_pre_pv_pois4 <- eca_pois4[[5]]$`p-value precursor`
shp$tmax_conf_idp_pre_pv_pois13 <- eca_pois13[[5]]$`p-value precursor`
shp$tmax_conf_idp_pre_pv_wtsur1 <- eca_wtsur1[[5]]$`p-value precursor`
shp$tmax_conf_idp_pre_pv_wtsur4 <- eca_wtsur4[[5]]$`p-value precursor`
shp$tmax_conf_idp_pre_pv_wtsur13 <- eca_wtsur13[[5]]$`p-value precursor`
shp$tmax_conf_idp_pre_pv_shuf1 <- eca_shuf1[[5]]$`p-value precursor`
shp$tmax_conf_idp_pre_pv_shuf4 <- eca_shuf4[[5]]$`p-value precursor`
shp$tmax_conf_idp_pre_pv_shuf13 <- eca_shuf13[[5]]$`p-value precursor`

shp$conf_idp_pre_pv_other1 <- conf_idp_stat1*1
shp$conf_idp_pre_pv_other4 <- conf_idp_stat4*1
shp$conf_idp_pre_pv_other13 <- conf_idp_stat13*1
shp$conf_idp_pre_pv_other1[is.nan(shp$conf_idp_pre_pv_other1)] <- -1
shp$conf_idp_pre_pv_other4[is.nan(shp$conf_idp_pre_pv_other4)] <- -1
shp$conf_idp_pre_pv_other13[is.nan(shp$conf_idp_pre_pv_other13)] <- -1

shp$rain_idp_pre_pv_other1 <- rain_idp_stat1*1
shp$rain_idp_pre_pv_other4 <- rain_idp_stat4*1
shp$rain_idp_pre_pv_other13 <- rain_idp_stat13*1
shp$rain_idp_pre_pv_other1[is.nan(shp$rain_idp_pre_pv_other1)] <- -1
shp$rain_idp_pre_pv_other4[is.nan(shp$rain_idp_pre_pv_other4)] <- -1
shp$rain_idp_pre_pv_other13[is.nan(shp$rain_idp_pre_pv_other13)] <- -1

shp$rain_conf_idp_pre_pv_other1 <- rain_conf_idp_stat1*1
shp$rain_conf_idp_pre_pv_other4 <- rain_conf_idp_stat4*1
shp$rain_conf_idp_pre_pv_other13 <- rain_conf_idp_stat13*1
shp$rain_conf_idp_pre_pv_other1[is.nan(shp$rain_conf_idp_pre_pv_other1)] <- -1
shp$rain_conf_idp_pre_pv_other4[is.nan(shp$rain_conf_idp_pre_pv_other4)] <- -1
shp$rain_conf_idp_pre_pv_other13[is.nan(shp$rain_conf_idp_pre_pv_other13)] <- -1

shp$tmax_idp_pre_pv_other1 <- tmax_idp_stat1*1
shp$tmax_idp_pre_pv_other4 <- tmax_idp_stat4*1
shp$tmax_idp_pre_pv_other13 <- tmax_idp_stat13*1
shp$tmax_idp_pre_pv_other1[is.nan(shp$tmax_idp_pre_pv_other1)] <- -1
shp$tmax_idp_pre_pv_other4[is.nan(shp$tmax_idp_pre_pv_other4)] <- -1
shp$tmax_idp_pre_pv_other13[is.nan(shp$tmax_idp_pre_pv_other13)] <- -1

shp$tmax_conf_idp_pre_pv_other1 <- tmax_conf_idp_stat1*1
shp$tmax_conf_idp_pre_pv_other4 <- tmax_conf_idp_stat4*1
shp$tmax_conf_idp_pre_pv_other13 <- tmax_conf_idp_stat13*1
shp$tmax_conf_idp_pre_pv_other1[is.nan(shp$tmax_conf_idp_pre_pv_other1)] <- -1
shp$tmax_conf_idp_pre_pv_other4[is.nan(shp$tmax_conf_idp_pre_pv_other4)] <- -1
shp$tmax_conf_idp_pre_pv_other13[is.nan(shp$tmax_conf_idp_pre_pv_other13)] <- -1
###

# Plot ECA results
  # Produce four plots (c: conflict, d: dry weather, h: hot weather, m: migration)
  # 1. c-m coincidence rates, d-m coincidence rates, d-c-m coincidence rates (coinratemap_rain....png)
  # 2. c-m coincidence rates, h-m coincidence rates, h-c-m coincidence rates (coinratemap_tmax....png)
  # 3. c-m coincidence occurrences, d-m coincidence occurrences, d-c-m coincidence occurrences (coinnummap_rain....png)
  # 4. c-m coincidence occurrences, h-m coincidence occurrences, h-c-m coincidence occurrences (coinnummap_tmax....png)
setwd(dir.out)
eca_plot(shp,1,"other") 
eca_plot(shp,4,"other")
eca_plot(shp,13,"other") # Figures 3 and 4
###
######

### CORRELATION BETWEEN BURSTINESS AND COINCIDENCE ###
# Plot correlation between burstiness and coincidence rates
setwd(dir.out)
corr_plot(BMOther,eca_pois1,1,0.05,nan_ctr = FALSE,"other")
corr_plot(BMOther,eca_pois4,4,0.05,nan_ctr = FALSE,"other")
corr_plot(BMOther,eca_pois13,13,0.05,nan_ctr = FALSE,"other") # Figure 5

# to include nan coincidence rates into correlation calculation, do the following:
  # corr_plot(BMOther,eca_pois1,1,0.05,nan_ctr = TRUE,"other")
  # corr_plot(BMOther,eca_pois4,4,0.05,nan_ctr = TRUE,"other")
  # corr_plot(BMOther,eca_pois13,13,0.05,nan_ctr = TRUE,"other")
###

# Print aggregate values
agg_cal(shp,timeDistOther) # aggregate values in Figures 3 and 4
###
######

### SUPPLEMENTARY MATERIALS ###
# Data Display
setwd(dir.out)
data_plot("data_mig.png",timeDistOther,"green",distName) # Supplementary Figure 1
data_plot("data_conf.png",data.matrix(conf_dat),"red",distName) # Supplementary Figure 2
data_plot("data_rain.png",data.matrix(rain_dat),"blue",distName) # Supplementary Figure 3
data_plot("data_temp.png",data.matrix(tmax_dat),"orange",distName) # Supplementary Figure 4
###

# Different population size thresholds for B-M calculations
timeDistOther25 <- ts_crt(outidp_dat,0.25,"mig")
timeDistOther50 <- ts_crt(outidp_dat,0.50,"mig")
timeDistOther75 <- ts_crt(outidp_dat,0.75,"mig")

BMOther25 <- BM_cal(timeDistOther25)
BMOther50 <- BM_cal(timeDistOther50)
BMOther75 <- BM_cal(timeDistOther75)

setwd(dir.out)
p_BM25 <- BM_plot(BMOther25)
p_BM25
ggsave("BM_25.png", height = 7 , width = 7) # Supplementary Figure 5A

p_BM50 <- BM_plot(BMOther50)
p_BM50
ggsave("BM_50.png", height = 7 , width = 7) # Supplementary Figure 5B

p_BM75 <- BM_plot(BMOther75)
p_BM75
ggsave("BM_75.png", height = 7 , width = 7) # Supplementary Figure 5C
###

# Null model (random shuffling)
no_ev <- colSums(timeDistOther)
null_time <- timeDistOther
null_time[,] <- 0
for (i in 1:74){
  idx <- sample(c(1:365), no_ev[i])
  null_time[idx,i] <- 1 
}
BMNull <- BM_cal(null_time)

setwd(dir.out)
p_BM_null <- BM_plot(BMNull)
p_BM_null
ggsave("BM_null.png", height = 7 , width = 7) # Supplementary Figure 6

# Calculation of burstiness considering finite event size
setwd(dir.out)
p_AM <- BM_plot(BMOther,finite = TRUE)
p_AM
ggsave("AM_other.png", height = 7 , width = 7) # Supplementary Figure 15

# Maps of expected coincidence rates and p-values
  # Clean data
shp$expected_conf_idp_exp_pois1 <- eca_pois1[[1]]$`Expected precursor rate`
shp$expected_conf_idp_exp_pois4 <- eca_pois4[[1]]$`Expected precursor rate`
shp$expected_conf_idp_exp_pois13 <- eca_pois13[[1]]$`Expected precursor rate`
shp$expected_conf_idp_exp_wtsur1 <- eca_wtsur1[[1]]$`Expected precursor rate`
shp$expected_conf_idp_exp_wtsur4 <- eca_wtsur4[[1]]$`Expected precursor rate`
shp$expected_conf_idp_exp_wtsur13 <- eca_wtsur13[[1]]$`Expected precursor rate`
shp$expected_conf_idp_exp_shuf1 <- eca_shuf1[[1]]$`Expected precursor rate`
shp$expected_conf_idp_exp_shuf4 <- eca_shuf4[[1]]$`Expected precursor rate`
shp$expected_conf_idp_exp_shuf13 <- eca_shuf13[[1]]$`Expected precursor rate`

shp$expected_rain_idp_exp_pois1 <- eca_pois1[[2]]$`Expected precursor rate`
shp$expected_rain_idp_exp_pois4 <- eca_pois4[[2]]$`Expected precursor rate`
shp$expected_rain_idp_exp_pois13 <- eca_pois13[[2]]$`Expected precursor rate`
shp$expected_rain_idp_exp_wtsur1 <- eca_wtsur1[[2]]$`Expected precursor rate`
shp$expected_rain_idp_exp_wtsur4 <- eca_wtsur4[[2]]$`Expected precursor rate`
shp$expected_rain_idp_exp_wtsur13 <- eca_wtsur13[[2]]$`Expected precursor rate`
shp$expected_rain_idp_exp_shuf1 <- eca_shuf1[[2]]$`Expected precursor rate`
shp$expected_rain_idp_exp_shuf4 <- eca_shuf4[[2]]$`Expected precursor rate`
shp$expected_rain_idp_exp_shuf13 <- eca_shuf13[[2]]$`Expected precursor rate`

shp$expected_rain_conf_idp_exp_pois1 <- eca_pois1[[3]]$`Expected precursor rate`
shp$expected_rain_conf_idp_exp_pois4 <- eca_pois4[[3]]$`Expected precursor rate`
shp$expected_rain_conf_idp_exp_pois13 <- eca_pois13[[3]]$`Expected precursor rate`
shp$expected_rain_conf_idp_exp_wtsur1 <- eca_wtsur1[[3]]$`Expected precursor rate`
shp$expected_rain_conf_idp_exp_wtsur4 <- eca_wtsur4[[3]]$`Expected precursor rate`
shp$expected_rain_conf_idp_exp_wtsur13 <- eca_wtsur13[[3]]$`Expected precursor rate`
shp$expected_rain_conf_idp_exp_shuf1 <- eca_shuf1[[3]]$`Expected precursor rate`
shp$expected_rain_conf_idp_exp_shuf4 <- eca_shuf4[[3]]$`Expected precursor rate`
shp$expected_rain_conf_idp_exp_shuf13 <- eca_shuf13[[3]]$`Expected precursor rate`

shp$expected_tmax_idp_exp_pois1 <- eca_pois1[[4]]$`Expected precursor rate`
shp$expected_tmax_idp_exp_pois4 <- eca_pois4[[4]]$`Expected precursor rate`
shp$expected_tmax_idp_exp_pois13 <- eca_pois13[[4]]$`Expected precursor rate`
shp$expected_tmax_idp_exp_wtsur1 <- eca_wtsur1[[4]]$`Expected precursor rate`
shp$expected_tmax_idp_exp_wtsur4 <- eca_wtsur4[[4]]$`Expected precursor rate`
shp$expected_tmax_idp_exp_wtsur13 <- eca_wtsur13[[4]]$`Expected precursor rate`
shp$expected_tmax_idp_exp_shuf1 <- eca_shuf1[[4]]$`Expected precursor rate`
shp$expected_tmax_idp_exp_shuf4 <- eca_shuf4[[4]]$`Expected precursor rate`
shp$expected_tmax_idp_exp_shuf13 <- eca_shuf13[[4]]$`Expected precursor rate`

shp$expected_tmax_conf_idp_exp_pois1 <- eca_pois1[[5]]$`Expected precursor rate`
shp$expected_tmax_conf_idp_exp_pois4 <- eca_pois4[[5]]$`Expected precursor rate`
shp$expected_tmax_conf_idp_exp_pois13 <- eca_pois13[[5]]$`Expected precursor rate`
shp$expected_tmax_conf_idp_exp_wtsur1 <- eca_wtsur1[[5]]$`Expected precursor rate`
shp$expected_tmax_conf_idp_exp_wtsur4 <- eca_wtsur4[[5]]$`Expected precursor rate`
shp$expected_tmax_conf_idp_exp_wtsur13 <- eca_wtsur13[[5]]$`Expected precursor rate`
shp$expected_tmax_conf_idp_exp_shuf1 <- eca_shuf1[[5]]$`Expected precursor rate`
shp$expected_tmax_conf_idp_exp_shuf4 <- eca_shuf4[[5]]$`Expected precursor rate`
shp$expected_tmax_conf_idp_exp_shuf13 <- eca_shuf13[[5]]$`Expected precursor rate`

# Plot maps
  # expected coincidence rates of null hypothesis
setwd(dir.out)
si_exp_plot(shp,1) 
si_exp_plot(shp,4)
si_exp_plot(shp,13) # Supplementary Figures 7-9 -- results may be slightly different due to randomness

si_pval_plot(shp,1) 
si_pval_plot(shp,4)
si_pval_plot(shp,13) # Supplementary Figures 10-13 -- results may be slightly different due to randomness

# lags
setwd(dir.out)
# dt = 1
ggplot() + aes(lag_out1[[1]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out1[[1]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
  labs(title="",x="Time lag (week)", y = "Density")+
    scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("conf_mig_lag1.png", height = 7 , width = 7) # Supplementary Figure 14A
print(mean(lag_out1[[1]]$X10))
  
ggplot() + aes(lag_out1[[2]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out1[[2]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
labs(title="",x="Time lag (week)", y = "Density")+
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("rain_mig_lag1.png", height = 7 , width = 7) # Supplementary Figure 14B
print(mean(lag_out1[[2]]$X10))

ggplot() + aes(lag_out1[[5]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out1[[5]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
labs(title="",x="Time lag (week)", y = "Density")+
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("tmax_mig_lag1.png", height = 7 , width = 7) # Supplementary Figure 14C
print(mean(lag_out1[[5]]$X10))

# dt = 4
ggplot() + aes(lag_out4[[1]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out4[[1]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
labs(title="",x="Time lag (week)", y = "Density")+
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("conf_mig_lag4.png", height = 7 , width = 7) # Supplementary Figure 14D
print(mean(lag_out4[[1]]$X10))

ggplot() + aes(lag_out4[[2]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out4[[2]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
labs(title="",x="Time lag (week)", y = "Density")+
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("rain_mig_lag4.png", height = 7 , width = 7) # Supplementary Figure 14E
print(mean(lag_out4[[2]]$X10))

ggplot() + aes(lag_out4[[5]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out4[[5]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
labs(title="",x="Time lag (week)", y = "Density")+
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("tmax_mig_lag4.png", height = 7 , width = 7) # Supplementary Figure 14F
print(mean(lag_out4[[5]]$X10))

# dt = 13
ggplot() + aes(lag_out13[[1]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out13[[1]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
labs(title="",x="Time lag (week)", y = "Density")+
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("conf_mig_lag13.png", height = 7 , width = 7) # Supplementary Figure 14G
print(mean(lag_out13[[1]]$X10))

ggplot() + aes(lag_out13[[2]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out13[[2]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
labs(title="",x="Time lag (week)", y = "Density")+
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("rain_mig_lag13.png", height = 7 , width = 7) # Supplementary Figure 14H
print(mean(lag_out13[[2]]$X10))

ggplot() + aes(lag_out13[[5]]$X10)+ geom_histogram(aes(y=..density..),binwidth=1, colour="black", fill="white")+
  theme_minimal()+geom_vline(xintercept = mean(lag_out13[[5]]$X10),linetype="dashed",linewidth = 1,colour = "red") +
labs(title="",x="Time lag (week)", y = "Density")+
  scale_x_continuous(breaks = seq(0, 8, by = 1)) + ylim(0,1.0)+
  theme(axis.text = element_text(size = 15),
        axis.title  = element_text(size = 15))
ggsave("tmax_mig_lag13.png", height = 7 , width = 7) # Supplementary Figure 14I
print(mean(lag_out13[[5]]$X10))

# B-M calculation of Mobility
timeDistSelf <- ts_crt(outidp_dat,0,"mob")
BMSelf <- BM_cal(timeDistSelf)

setwd(dir.out)
p_BM2 <- BM_plot(BMSelf)
p_BM2
ggsave("BM_self.png", height = 7 , width = 7) # Supplementary Figure 16

# ECA for mobility
  # lag
lag_out_self1 <- lag_eca(timeDistSelf,conf_dat,rain_dat,tmax_dat,9,1)
lag_out_self4 <- lag_eca(timeDistSelf,conf_dat,rain_dat,tmax_dat,9,4)
lag_out_self13 <- lag_eca(timeDistSelf,conf_dat,rain_dat,tmax_dat,9,13)

  # Run ECA
eca_self1 <- real_eca(timeDistSelf,conf_dat,rain_dat,tmax_dat,9,1,lag_out_self1,"poisson")
eca_self4 <- real_eca(timeDistSelf,conf_dat,rain_dat,tmax_dat,9,4,lag_out_self4,"poisson")
eca_self13 <- real_eca(timeDistSelf,conf_dat,rain_dat,tmax_dat,9,13,lag_out_self13,"poisson")

shp$conf_idp_pre_self1 <- eca_self1[[1]]$`precursor coincidence rate`
shp$conf_idp_pre_self4 <- eca_self4[[1]]$`precursor coincidence rate`
shp$conf_idp_pre_self13 <- eca_self13[[1]]$`precursor coincidence rate`

shp$rain_idp_pre_self1 <- eca_self1[[2]]$`precursor coincidence rate`
shp$rain_idp_pre_self4 <- eca_self4[[2]]$`precursor coincidence rate`
shp$rain_idp_pre_self13 <- eca_self13[[2]]$`precursor coincidence rate`

shp$rain_conf_idp_pre_self1 <- eca_self1[[3]]$`precursor coincidence rate`
shp$rain_conf_idp_pre_self4 <- eca_self4[[3]]$`precursor coincidence rate`
shp$rain_conf_idp_pre_self13 <- eca_self13[[3]]$`precursor coincidence rate`

shp$tmax_idp_pre_self1 <- eca_self1[[4]]$`precursor coincidence rate`
shp$tmax_idp_pre_self4 <- eca_self4[[4]]$`precursor coincidence rate`
shp$tmax_idp_pre_self13 <- eca_self13[[4]]$`precursor coincidence rate`

shp$tmax_conf_idp_pre_self1 <- eca_self1[[5]]$`precursor coincidence rate`
shp$tmax_conf_idp_pre_self4 <- eca_self4[[5]]$`precursor coincidence rate`
shp$tmax_conf_idp_pre_self13 <- eca_self13[[5]]$`precursor coincidence rate`

setwd(dir.out)
eca_plot(shp,1,"self", occ = FALSE)
eca_plot(shp,4,"self", occ = FALSE)
eca_plot(shp,13,"self", occ = FALSE) # Supplementary Figure 17

corr_plot(BMSelf,eca_self1,1,0.05,nan_ctr = FALSE,"self")
corr_plot(BMSelf,eca_self4,4,0.05,nan_ctr = FALSE,"self")
corr_plot(BMSelf,eca_self13,13,0.05,nan_ctr = FALSE,"self") # Supplementary Figure 18


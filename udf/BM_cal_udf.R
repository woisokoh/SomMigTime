ts_crt <- function(outidp_dat,percentile,type){
  outidp_dat <- subset(outidp_dat,IDPPop > 0)
  
  timeDistOther <- data.frame(matrix(0, ncol = 74, nrow = 365))
  popDistOther <- data.frame(matrix(0, ncol = 74, nrow = 365))
  for (i in 1:nrow(outidp_dat)){
    if (type == "mig" & outidp_dat[i,2] != outidp_dat[i,3]){
      timeDistOther[outidp_dat[i,1],outidp_dat[i,2]] = 1
      popDistOther[outidp_dat[i,1],outidp_dat[i,2]] = popDistOther[outidp_dat[i,1],outidp_dat[i,2]]+ outidp_dat[i,6]
    }
    if (type == "mob" & outidp_dat[i,2] == outidp_dat[i,3]){
      timeDistOther[outidp_dat[i,1],outidp_dat[i,2]] = 1
      popDistOther[outidp_dat[i,1],outidp_dat[i,2]] = popDistOther[outidp_dat[i,1],outidp_dat[i,2]]+ outidp_dat[i,6]
    }
  }
  colnames(timeDistOther) <- distName[[1]]
  colnames(popDistOther) <- distName[[1]]
  cutPt <- quantile(popDistOther[popDistOther > 0],percentile) # change percentile over here
  timeDistOther <- (popDistOther > cutPt)*1
  return(timeDistOther)
}

BM_cal <- function(timeDistOther){
  BMOther <- data.frame(matrix(0, ncol = 5, nrow = 74))
  for (d in 1:74){
    if (sum(timeDistOther[,d]) > 3){
      BM <- burstiness_GB(data.frame(timeDistOther[,d]))
      A <- burstiness_KJ(data.frame(timeDistOther[,d]))
      BMOther[d,1]<-BM[2];BMOther[d,2]<-BM[1];BMOther[d,3]<-A
      BMOther[d,4]<- sum(timeDistOther[,d]);BMOther[d,5]<- 'pass'
    }else{
      BMOther[d,1]<-NA;BMOther[d,2]<-NA;BMOther[d,3]<-NA
      BMOther[d,4]<- sum(timeDistOther[,d]);BMOther[d,5]<-NA
    }
  }
  colnames(BMOther) <- c('M','B','A','eventno','boundCon')
  rownames(BMOther) <- distName[[1]]
  return(BMOther)
}

BM_plot <- function(BMOther, finite = FALSE){
  if (finite == FALSE){
    x <- BMOther[,c(1:2,5)]
  }else{
    x <- BMOther[,c(1,3,5)]
    colnames(x) <- c("M","B","boundCon")
  }
  
  x <- x[!is.na(x$B),]
  Mbig <- order(x$M,decreasing = TRUE)[1:7]
  Msmall <- order(x$M,decreasing = FALSE)[1:7]
  Bbig <- order(x$B,decreasing = TRUE)[1:7]
  Bsmall <- order(x$B,decreasing = FALSE)[1:7]
  
  idx1 <- Mbig[x$boundCon[Mbig] == 'pass']
  idx2 <- Msmall[x$boundCon[Msmall] == 'pass']
  idx3 <- Bbig[x$boundCon[Bbig] == 'pass']
  idx4 <- Bsmall[x$boundCon[Bsmall] == 'pass']
  idxl <- unique(c(idx1,idx2,idx3,idx4))
  x$lab <- ""
  x$lab[idxl] <- rownames(x)[idxl]
  
  Mavg <- mean(x$M[x$boundCon == 'pass'])
  Bavg <- mean(x$B[x$boundCon == 'pass'])
  
  plot1 <- ggplot(x, aes(x = M, y = B, color = boundCon, label = lab)) + 
    geom_point(aes(color = boundCon), size = 3) + 
    geom_point(shape = 1, color = "black", size = 3) + 
    geom_rug() + 
    scale_y_continuous(name = "B", 
                       limits = c(-1, 1), 
                       expand = c(0, 0)) + 
    scale_x_continuous(name = "M", 
                       limits = c(-1, 1), 
                       expand = c(0, 0)) + 
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 15)
    ) +
    geom_hline(yintercept = Bavg, linetype = "longdash",color = "red") +
    geom_vline(xintercept = Mavg, linetype = "longdash",color = "red") +
    #geom_text(data = x, color = "red",
    #          mapping = aes(x=-0.9, y=Bavg-0.05, 
    #                        label=paste("B_avg"," = ",
    #                                    format(round(Bavg, 2),
    #                                           nsmall=2))), 
    #          hjust=0) +
    #geom_text(data = x, color = "red", 
    #          mapping = aes(x=Mavg+0.03, y=-0.85, 
    #                        label=paste("M_avg"," = ",
    #                                    format(round(Mavg, 2),
    #                                           nsmall=2))), 
  #          hjust=0) + 
  geom_text_repel()
  
  theme0 <- function(...) theme( legend.position = "none",
                                 panel.background = element_blank(),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.margin = unit(0,"null"),
                                 axis.ticks = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_blank(),
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 axis.ticks.length = unit(0,"null"),
                                 axis.ticks.margin = unit(0,"null"),
                                 panel.border=element_rect(color=NA),...)
  
  dens1 <- ggplot(x, aes(x = M, fill = boundCon,color = boundCon)) +
    geom_density(alpha = 0.4) + 
    theme_bw() +
    theme0(plot.margin = unit(c(1,0,-0.48,2.2),"lines")) +
    theme(legend.position = "none")+ 
    scale_x_continuous(name = "M", 
                       limits = c(-1, 1), expand = c(0, 0))
  
  dens2 <- ggplot(x, aes(x = B, fill = boundCon,color = boundCon)) +
    geom_density(alpha = 0.4) + 
    theme_bw() +
    theme0(plot.margin = unit(c(0,1,1.2,-0.48),"lines")) +
    theme(legend.position = "none") + 
    coord_flip()+ 
    scale_x_continuous(name = "B", 
                       limits = c(-1, 1), expand = c(0, 0))
  
  p <- dens1 + plot_spacer() + plot1 + dens2 + 
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  # print(paste("B_avg"," = ",format(round(Bavg, 2),nsmall=2)))
  # print(paste("M_avg"," = ",format(round(Mavg, 2),nsmall=2)))
  
  return(p)
}

data_plot <- function(filename,dat,col,distName){
  png(filename=filename,width = 6000, height = 4000, res = 300)
  par(mar = c(2, 10, 2, 2)) # Set the margin on all sides to 2
  image.plot(1:nrow(dat), 1:ncol(dat), dat,col=c("white",col), axes = FALSE,
             xlab ="", ylab = "")
  for (i in 1:73){abline(h=i+0.5, col="black")}
  axis(2, at = seq(1, ncol(dat), by = 1), labels = t(distName), las = 1)
  axis(1, at = seq(1, nrow(dat), by = 10), labels = seq(1, nrow(dat), by = 10))
  dev.off()
}

ts_plot <- function(district, distName,timeDistOther){
  distNo <- which(distName == district)
  tsd <- which(timeDistOther[,distNo] == 1)
  p <- ggplot() + 
    geom_vline(xintercept = tsd,color = "blue",linewidth = 0.1) + 
    theme_bw() +
    scale_x_continuous(name = "weeks", 
                       limits = c(1, 365), 
                       expand = c(0, 0))
  return(p)
}

iet_plot <- function(district, distName,timeDistOther,fit){
  distNo <- which(distName == district)
  tsd <- which(timeDistOther[,distNo] == 1)
  if (fit == "power-law"){
    dd <- displ$new(diff(tsd))
  }
  if (fit == "exponential"){
    dd <- disexp$new(diff(tsd))
  }
  est <- estimate_pars(dd) # pars: alpha
  est <- estimate_xmin(dd) # xmin
  dd$setXmin(est)
  return(dd)
}
daysOfMonth <- function(year, month){
  days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  day <- days[month]
  if((year %% 4 == 0 && year %% 100 != 0 || year %% 400 ) && month ==2) day <- 29
  
  day 
}

dataArange_Before <- function(data, Segment = "04-03"){
  
  date_string <- sprintf("%d-%d-%d", data[, 1], data[, 2], data[, 3])
  date_daily <- as.Date(date_string)
  n <- nrow(data)
  
  if (Segment == "month"){
    beginMonth <- 1
    endMonth <- 12
    
    id_begin <- 1; id_end <- n
    
    
    date_daily <- date_daily[id_begin:id_end]
    data_trim <- cbind(data, labels = format(date_daily, "%Y-%m"))
  }
  
  clim.data <- data.frame(SpecialValue(data_trim), labels = data_trim$labels)
  
  clim.data
}

clim.quantile <- function(X, Period = "1960-2013"){
  data <- SpecialValue(X)
  date_string <- sprintf("%d-%d-%d", data[, 1], data[, 2], data[, 3])
  date_daily <- as.Date(date_string)
  data <- zoo(data, date_daily)
  BeginYear <- as.numeric(substr(Period, 1, 4)); EndYear <- as.numeric(substr(Period, 6, 9))
  
  standardPeriod <- seq(as.Date(sprintf("%d-01-01", BeginYear)), as.Date(sprintf("%d-12-31", EndYear)), by = "day")
  data_trim <- window(data, standardPeriod)
  if (length(data_trim) == 0){
    warning("Nonexistent data！")
    return()
  }
  
  filterSmooth <- function(x, width = 5){
    x <- as.numeric(x);n <- length(x)
    X_temp <- cbind(x[1:(n - width + 1)], x[2:(n - width + 2)], x[3:(n - width + 3)], x[4:(n - width + 4)], x[5:(n - width + 5)])
    apply(X_temp, 1, mean, na.rm = T)
  }
  
  Tmax.filter <- filterSmooth(data_trim[, 5])
  Tmin.filter <- filterSmooth(data_trim[, 6])
  
  
  Precp <- data_trim[, 7]; Precp_trim <- Precp[Precp >= 1.0]
  result <- c(quantile(Tmax.filter, c(0.1, 0.9)), quantile(Tmin.filter, c(0.1, 0.9)), quantile(Precp_trim, c(0.95, 0.99)))
}


ContinueTag <- function(X){
  Nind <- length(X)
  if (Nind == 0) return(0) 
  if (Nind == 1) return(1)  
  
  Tag <- numeric(Nind) 
  nEvent = 1;Tag[1] <- 1 
  
  for (i in 1:(Nind-1)){
    if (X[i+1]!=X[i]+1) nEvent<-nEvent+1
    Tag[i+1] <- nEvent
  }
  Tag##quickly return
}

clim.GSL <- function(Taver){
  if (length(Taver) < 365) {
    warning("year?")
    return(NA)
  }
  Nmid <- floor(length(Taver)/2); N <- length(Taver)
  
  Id_begin <- which(Taver > 5)
  Tag <- ContinueTag(Id_begin)
  segment.length <- sapply(1:Tag[length(Tag)], function(i) length(which(Tag == i)))
  TagId <- which(segment.length >= 6)[1]            
  point.begin <- Id_begin[which(Tag == TagId)[1]]   错
  
  Id_end <- which(Taver[(Nmid+1):N] < 5) + Nmid
  Tag <- ContinueTag(Id_end)
  segment.length <- sapply(1:Tag[length(Tag)], function(i) length(which(Tag == i)))
  TagId <- which(segment.length >= 6); TagId <- TagId[1]         
  point.end <- Id_end[which(Tag == TagId)]; point.end <- point.end[1]
  if ((length(point.begin)==0) | (length(point.end)==0 ))  return(NA)
  if (is.na(point.begin) | is.na(point.end))  return(NA)
  
  if (point.begin >= 183) point.begin <- NA
  if (point.end < 183) end.point <- NA
  if (is.na(point.begin) | is.na(point.end)) warning("Please check")
  point.end - point.begin
}


clim.RX <- function(precp){
  precp <- as.matrix(precp); Ndays <- length(precp)
  data_R5 <- cbind(precp[1:(Ndays-4)], precp[2:(Ndays-3)],precp[3:(Ndays-2)],precp[4:(Ndays-1)],precp[5:(Ndays)])
  Rx5 <- max(apply(data_R5, 1, sum, na.rm=T))
  Rx1 <- max(precp)
  data.frame(Rx1 = Rx1, Rx5 = Rx5)
}


clim.SDII <- function(precp) sum(precp[which(precp >= 1)])/length(which(precp >= 1))

clim.RRN <- function(precp, nm = c(10, 20, 25)){
  rrn <- sapply(nm, function(x) length(which(precp >= x)))
  rrn <- data.frame(t(rrn))
  colnames(rrn) <- paste("RR", nm, sep = "");rrn
}

clim.CDD <- function(precp, item = "drought"){
  if (!item %in% c("drought", "wet")) stop("item param must bu 'drought' or 'wet'!")
  if(item == "drought") cdd.id <- which(precp < 1) else cdd.id <- which(precp >= 1)
  
  Tag <- ContinueTag(cdd.id)
  cdd <- max(sapply(1:Tag[length(Tag)],function(i) length(which(Tag == i)))) 
  cdd##quickly return
}

clim.Rquantile <- function(precp, quantile.standard){
  R_days <- sapply(1:2, function(i) length(which(precp > quantile.standard[i])))#list1 for q95, list2 for q99, days
  R_precp <- sapply(1:2, function(i) sum(precp[which(precp > quantile.standard[i])]))#list1 for q95, list2 for q99, calculate
  data.frame(R95D = R_days[1], R99D = R_days[2], R95P = R_precp[1], R99p = R_precp[2])##quickly return 
}

clim.PRCPTOT <- function(precp) sum(precp[which(precp >= 1)])

precpIndice <- function(X, Quantile = clim_quantile){ 
  precp <- X$precp
  Rx <- clim.RX(precp)
  SDII <- clim.SDII(precp)
  RRN <- clim.RRN(precp, nm = c(10, 20, 25))
  CDD <- clim.CDD(precp, item = "drought")
  CWD <- clim.CDD(precp, item = "wet")
  Rquantile <- clim.Rquantile(precp, Quantile[c("RR.95th", "RR.99th")])
  PRCPTOT <- clim.PRCPTOT(precp)
  data.frame(Rx, SDII, RRN, CDD, CWD, Rquantile, PRCPTOT)##quickly return 13 extreme precp indice
}
allIndice <- function(X, Quantile = clim_quantile) data.frame(TIndice(X, Quantile), precpIndice(X, Quantile))
Clim.Indice <- function(X, quantile = clim_quantile, index = "all"){
  index <- c("FD", "SU", "ID", "TR", "GSL", "TXx","TNx", "TXn", "TNn", "Tthp", "WSDI", "CSDI") 
}

require(zoo)
require("PCICt")
require(xlsx)
require(climdex.pcic)
library(devtools)
library(jsonlite)
library(stringr)
library(tidyr)
library(dplyr)
library(RJSONIO)
library(gtools)

rm(list = ls())
setwd("")
source('Climate_Index.R', encoding = 'GB2312', echo=TRUE)
setwd("")
file_list <- list.files(pattern = ".txt")
mixedsort(file_list)

n <- length(file_list)
c <- read.table("",header = FALSE)
for (i in 1:1326){
  
  #data <- read.table(file_list[i],header = TRUE)
  
  #tmax.dates <- as.PCICt(do.call(paste, data[,1:3]), format = "%Y%m%d", cal = "gregorian")
  #tmin.dates <- as.PCICt(do.call(paste, data[,1:3]), format = "%Y%m%d", cal = "gregorian")
  #tavg.dates <- as.PCICt(do.call(paste, data[,1:3]), format = "%Y%m%d", cal = "gregorian")
  #prec.dates <- as.PCICt(do.call(paste, data[,1:3]), format = "%Y%m%d", cal = "gregorian")
  
  # year --------------------------------------------------------------------
  #ci <- climdexInput.raw(data[, 5], data[, 6], data[, 7], tmax.dates, tmin.dates, prec.dates, base.range = c(2015,2100),tavg = data[, 4], tavg.dates = tavg.dates)
  #CDD<-climdex.cdd(ci, spells.can.span.years = FALSE)
  #CWD<-climdex.cwd(ci, spells.can.span.years = FALSE)
  #R10mm<-climdex.r10mm(ci)
  #R20mm<-climdex.r20mm(ci)
  #R95ptot<-climdex.r95ptot(ci)
  #r99ptot<-climdex.r99ptot(ci)
  #Rx1day_y<-climdex.rx1day(ci,freq = "annual")
  #Rx5day_y<-climdex.rx5day(ci,freq = "annual")
  #SDII<-climdex.sdii(ci) 
  #PRCPTOT<-climdex.prcptot(ci)
  #prindex<-data.frame(CDD, CWD, R10mm, R20mm,R95ptot,r99ptot,Rx1day_y,Rx5day_y,SDII,PRCPTOT)
  
  #setwd("")
  #write.table(prindex,file=paste("year",i,".txt"),sep='\t',col.names = TRUE,quote=F,row.names = FALSE)
  #setwd("")
  #fnames <- dir("", pattern = file_list[i], full.names = T)
  clim.OriginalData <- read.table(file_list[c[i,1]],header = TRUE)
  # month ------------------------------------------------------------------
  
  clim_data <- dataArange_Before(clim.OriginalData, Segment = "month")
  clim_dataLs <- split(clim_data, clim_data$labels) 
  clim_quantile <- clim.quantile(clim.OriginalData, Period = "2015-2100")    
  AIndice <-lapply(clim_dataLs, precpIndice)        
  AIndice <- do.call(rbind.data.frame, AIndice)                   
  write.table(AIndice,file=paste("","month",c[i,1],".txt"),sep='\t',col.names = TRUE,quote=F,row.names = FALSE)
  # summer ------------------------------------------------------------------
  #clim_data <- dataArange_Before(clim.OriginalData, Segment = "06-08")
  #clim_dataLs <- split(clim_data, clim_data$labels) 
  #clim_quantile <- clim.quantile(clim.OriginalData, Period = "2015-2100")   
  #AIndice <-lapply(clim_dataLs, precpIndice)        
  #AIndice <- do.call(rbind.data.frame, AIndice)                   
  
  #write.table(AIndice,file=paste("","summer",i,".txt"),sep='\t',col.names = TRUE,quote=F,row.names = FALSE)
  # autumn ------------------------------------------------------------------
  #clim_data <- dataArange_Before(clim.OriginalData, Segment = "09-11")
  #clim_dataLs <- split(clim_data, clim_data$labels) 
  #clim_quantile <- clim.quantile(clim.OriginalData, Period = "2015-2100")  
  #AIndice <-lapply(clim_dataLs, precpIndice)       
  #AIndice <- do.call(rbind.data.frame, AIndice)                   
  #write.table(AIndice,file=paste("","autumn",i,".txt"),sep='\t',col.names = TRUE,quote=F,row.names = FALSE)
}
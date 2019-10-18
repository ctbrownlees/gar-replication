rm(list=ls())
library(lubridate)
library(readxl)

d     <- read_excel('../../data/raw/wui.xlsx',sheet='T2')
dates <- parse_date_time(d$year,orders = 'y!*q!*')
cnames <- names(read.csv('../../data/raw/keep_countries.csv'))

dmatched <- data.frame(d[,colnames(d) %in% cnames])
cnames[which(!(cnames %in% colnames(d) ))] # To see whos out

dout <- data.frame(matrix(NA,nrow=nrow(d),ncol=length(cnames)+1))
colnames(dout) <- c('Dates',cnames)
global <- rowMeans(dmatched,na.rm = T)
dout[ , colnames(dout) %in% colnames(dmatched) ] <- dmatched
dout[ , !(colnames(dout) %in% colnames(dmatched)) ] <- global
dout$Dates <- dates

write.csv(dout,'../../data/clean/wui.csv')

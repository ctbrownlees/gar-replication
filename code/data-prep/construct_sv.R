rm(list=ls())
library(lubridate)

files <- dir('../../data/final/realized_vol/')
files <- files[grepl('_rv_',files)]

dates <- parse_date_time(as.character(d$Dates),orders = '%d.%m.%Y')
cnames <- names(read.csv('../../data/names.csv'))

dmatched <- data.frame(d[,colnames(d) %in% cnames])
cdmatched <- colnames(dmatched)
cnames[which(!(cnames %in% colnames(d) ))] # To see whos out

dout <- data.frame(matrix(NA,nrow=nrow(d)-1,ncol=length(cnames)+1))
colnames(dout) <- c('Dates',cnames)
dmatched <- sapply(1:ncol(dmatched),FUN = function(x) diff(log(dmatched[,x]),na.rm = T)*100 )
global <- rowMeans(dmatched,na.rm=T)
dout[ , colnames(dout) %in% cdmatched ] <- dmatched
dout$ISL <- global
dout$Dates <- dates[2:length(dates)] +days(1)
dout[,-1] <- scale(dout[,-1])
for (i in 2:ncol(dout)){
 dout[is.na(dout[,i]),i] <- global[is.na(dout[,i])]
}
write.csv(dout,'../../data/screening/creditgdp.csv')

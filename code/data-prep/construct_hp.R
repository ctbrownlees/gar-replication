rm(list=ls())
library(lubridate)

d     <- read.csv('../../data/raw/hp.csv')
dates <- parse_date_time(as.character(d$X),orders = '%d.%m.%Y')
cnames <- names(read.csv('../../data/raw/keep_countries.csv'))

dmatched <- data.frame(d[,colnames(d) %in% cnames])
cdmatched <- colnames(dmatched)
cnames[which(!(cnames %in% colnames(d) ))] # To see whos out

# fix japan
nas    <- which(is.na(d$JPN))
interp <- (d$JPN[nas-1] +  d$JPN[nas+1] )/2
d$JPN[nas] <- interp
#

dout <- data.frame(matrix(NA,nrow=nrow(d)-1,ncol=length(cnames)+1))
colnames(dout) <- c('Dates',cnames)
dmatched <- sapply(1:ncol(dmatched),FUN = function(x) diff(log(dmatched[,x]),na.rm = T)*100 )
global <- scale( rowMeans(scale(dmatched),na.rm=T) )
dout[ , colnames(dout) %in% cdmatched ] <- dmatched
dout$ISL <- global
dout$Dates <- dates[2:length(dates)] +days(1)
# dout[,-1] <- scale(dout[,-1])
for (i in 2:ncol(dout)){
  if (sum(is.na(dout[,i])) == length(dout[,i])){
    dout[is.na(dout[,i]),i] <- global[is.na(dout[,i])]
  }else{
    dout[is.na(dout[,i]),i] <- global[is.na(dout[,i])] * sd(dout[,i], na.rm = TRUE)
  }
}
dout <- dout[rowSums(is.na(dout))==0,]
write.csv(dout,'../../data/clean/hp.csv')

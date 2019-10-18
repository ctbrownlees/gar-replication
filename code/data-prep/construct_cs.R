rm(list=ls())
library(lubridate)
library(dplyr)

d <-  read.csv('../../data/raw/BAA10Y.csv',stringsAsFactors = F)

# Interpolate missing
baa      <- as.numeric(d$BAA10Y)
i=0
while( sum(is.na(baa))>0 ){
  i        = i+1
  nas      <- which(is.na(baa))
  baa[nas] <- (baa[nas-i]  + baa[nas+i])/2 
}

dframe <- cbind.data.frame('date'= parse_date_time(as.character(d$DATE[1:length(d$DATE)]),orders = 'Y!-m!*-d!'),
                             'ret' = (baa) )
quarterly <- dframe %>% group_by(dates=ceiling_date(date, "quarter")) %>%
    summarize(ret=tail(ret,1))
  
write.csv(quarterly,sprintf('../../data/clean/cs.csv'))

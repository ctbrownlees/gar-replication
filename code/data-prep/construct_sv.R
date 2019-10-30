rm(list=ls())
library(lubridate)
library(dplyr)

d    <- read.csv(sprintf('../../data/raw/GSPC.csv',file),stringsAsFactors = F)
dframe <- cbind.data.frame('date'= parse_date_time(as.character(d$Date[2:length(d$Date)]),orders = 'Y!-m!*-d!'),
                             'ret' = diff(log(as.numeric(d$Adj.Close)))*100 )
quarterly <- dframe %>% group_by(dates=floor_date(date, "quarter")) %>%
  summarize(ret=mean((ret^2),na.rm=TRUE))
quarterly <- quarterly[sapply(1:nrow(quarterly),function(w) !any(is.na(quarterly[w,]))),]

write.csv(quarterly,'../../data/clean/sv.csv')

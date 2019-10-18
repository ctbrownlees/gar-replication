rm(list=ls())
library(dplyr)
library(lubridate)
library(readxl)
d <- read_excel('../../data/raw/gpr.xlsx')

dates <- as.Date(as.numeric( d$Date[1:(length(d$Date)-3)]) ,origin = "1899-12-30")
gpr   <- d$GPR

dframe <- cbind.data.frame('date' = dates,'gpr' = unname(d[1:(nrow(d)-3),2]) )
  
quarterly <- dframe %>% group_by(dates=floor_date(date, "quarter")) %>%
    summarize(gpr=mean(gpr)/100)
quarterly <- quarterly[quarterly$dates > "1972-11-01" & quarterly$dates < "2017-01-01" ,]
quarterly <- quarterly[sapply(1:nrow(quarterly),function(w) !any(is.na(quarterly[w,]))),]

  
write.csv(quarterly,'../../data/clean/gpr.csv')

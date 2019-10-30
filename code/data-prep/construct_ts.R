rm(list=ls())
# Require lubridate to parse dates easier.
library(lubridate) 

# Set path to save data, if no saving required, set to NULL
save.path = '../../data/clean/'

# Load data from OECD.
lt <- read.csv('../../data/raw/lt_interest.csv')
st <- read.csv('../../data/raw/st_interest.csv')

lt <- lt[lt$FREQUENCY=='Q',]
st <- st[st$FREQUENCY=='Q',]

# Get country names
ctr_lt  <- unique(as.character(lt$LOCATION))
ctr_st  <- unique(as.character(st$LOCATION))
ctr_smp <- names(read.csv('../../data/raw/keep_countries.csv'))
# Take US to set lengths
d_st      <- st[st$LOCATION=="USA",]
d_lt      <- lt[lt$LOCATION=="USA",]

# Create data frame with N+1 to keep dates.
dframe_st <- data.frame( matrix(NA, nrow = nrow(d_st) ,ncol=length(ctr_smp)+1)) 
colnames(dframe_st) <- c('Dates',ctr_smp)
dframe_st[,'Dates'] <- as.character(d_st$TIME) 

dframe_lt <- data.frame( matrix(NA, nrow = nrow(d_st) ,ncol=length(ctr_smp)+1)) 
colnames(dframe_lt) <- c('Dates',ctr_smp)
dframe_lt[,'Dates'] <- as.character(d_st$TIME) 

for( i in 1:( length(ctr_smp) )){
  d_st   <- st[ st$LOCATION == ctr_smp[ i ] , ]
  date   <- as.character(d_st$TIME,orders="y-q") 
  common.dates <- intersect(date,dframe_st[,1])
  dframe_st[ dframe_st[ , 'Dates'] %in% common.dates, ctr_smp[ i ] ] <- d_st$Value[ date %in% common.dates ] 
  
  d_lt   <- lt[ lt$LOCATION == ctr_smp[ i ] , ]
  date   <- as.character(d_lt$TIME,orders="y-q") 
  common.dates <- intersect(date,dframe_st[,1])
  dframe_lt[ dframe_lt[ , 'Dates'] %in% common.dates, ctr_smp[ i ] ] <- d_lt$Value[ date %in% common.dates ] 

}
# Spread
d_spread <- dframe_lt[,-1] - dframe_st[,-1]
d_spread <- cbind.data.frame('Dates'=dframe_lt[,1],d_spread)
# Filling out NAs
global   <- scale(rowMeans(scale(d_spread[,-1]),na.rm=TRUE))
for (i in 1:length(ctr_smp)){
  d_spread[is.na(d_spread[,1+i]),1+i] <- scale(global[is.na(d_spread[,1+i])]) * sd(d_spread[,1+i],na.rm=T)
}

# Dates
dates_out <-  parse_date_time( dframe_st[,1] , order='y-q')
d_spread[,1]  <- dates_out

write.csv(d_spread,'../../data/clean/ts.csv')

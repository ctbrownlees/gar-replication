rm(list=ls())
# Require lubridate to parse dates easier.
library(lubridate) 
# Set path to save data, if no saving required, set to NULL
save.path = '../../data/clean/'

# Load data from OECD.
data <- read.csv('../../data/raw/gdp.csv')
data <- data[data$FREQUENCY=='Q',] # Keep quarterly Data
data <- data[data$MEASURE=='PC_CHGPP',]       # Keep quarter on quarter 
# Get country names
countries <- unique(as.character(data$LOCATION))

# Take US to set lengths
d      <- data[data$LOCATION=="USA",]

# Create data frame with N+1 to keep dates.
dframe <- data.frame( matrix(NA, nrow = nrow(d) ,ncol=length(countries)+1))
colnames(dframe) <- c('Dates',countries)
dframe[,'Dates'] <- as.character(d$TIME) 

for( i in 1:( length(countries) )){
  d      <- data[ data$LOCATION == countries[ i ] , ]
  date   <- as.character(d$TIME,orders="y-q") 
  common.dates <- intersect(date,dframe[,1])
  dframe[ dframe[ , 'Dates'] %in% common.dates, i + 1 ] <- d$Value 
}

cnames <- read.csv('../../data/raw/keep_countries.csv')
cnames <- colnames(cnames)
dframe <- dframe[,colnames(dframe) %in% c('Dates',cnames)]
if(!is.null(save.path)) write.csv(dframe, paste0(save.path , 'gdp.csv'))


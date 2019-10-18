rm(list=ls())
library(readxl)
library(lubridate)
library(countrycode)
## First, we do "Rest of the countries"
name.range <- "A1:CU1"
names <- read_excel("../../data/raw/nfci.xlsx",sheet = 3,range = name.range)
names.idx <- seq(1,ncol(names),7)
names <- colnames(names)[names.idx]
data.range <- "A2:CY178"
data  <- read_excel("../../data/raw/nfci.xlsx",sheet = 3,range = data.range,)
idx <- names.idx + 4
dates <- data[,names.idx]
fci   <-  data[, idx]
dframe <- data.frame(dates[,ncol(dates)])
dframe <- cbind.data.frame(dframe,matrix(NA,ncol=ncol(fci),nrow=nrow(dframe)))

for (i in 1:ncol(fci)){
  dframe[dframe[,1] %in% dates[[i]] , i+1] <- fci[[i]][dates[[i]] %in% dframe[,1]]
}
colnames(dframe) <- c("Dates",names)

## Next, we do "Countries included in fig_3.2"
name.range <- "A1:AJ1"
names <- read_excel("../../data/raw/nfci.xlsx",sheet = 2,range = name.range)
names.idx <- seq(1,ncol(names),7)
names <- gsub("^[0-9]. ",x=colnames(names)[names.idx],"")
data.range <- "A2:AN178"
data  <- read_excel("../../data/raw/nfci.xlsx",sheet = 2,range = data.range,)
idx <- names.idx + 4
dates <- data.frame( data[,names.idx] )
dates <-  lapply(1:length(dates),FUN=function(w) parse_date_time(dates[,w],orders="y!*q!*") ) 
fci   <-  data[, idx]
ncol.pre <- ncol(dframe)
dframe <- cbind.data.frame(dframe,matrix(NA,ncol=ncol(fci),nrow=nrow(dframe)))
ncol.pos <- ncol(dframe)

for (i in 1:ncol(fci)){
  dframe[parse_date_time(dframe[,1],'yq') %in% dates[[i]] , ncol.pre + i] <-
    fci[[i]][dates[[i]] %in% parse_date_time(dframe[,1],'yq')]
}
colnames(dframe) <- c(colnames(dframe)[1:ncol.pre] , names)

## Sorting by names in our sample.
cnames <- (read.csv('../../data/raw/keep_countries.csv',header = T))
cnames <- colnames(cnames)

# Translate IMF names and find intersection
imf.names <- countrycode(colnames(dframe[,-1]),origin = 'country.name',destination = 'iso3c')
nfci      <- data.frame(matrix(NA,nrow=nrow(dframe),ncol=length(cnames)))
colnames(nfci) <- cnames
for (j in 2:ncol(dframe)){
  country       <- countrycode(colnames(dframe)[j],origin = 'country.name', destination = 'iso3c')
  if (country %in% cnames) nfci[country] <- dframe[,j]
}

global.nfci  <- scale(rowMeans(dframe[,-1],na.rm=TRUE)) # To use information from all NFCIs, not just the ones we keep.
nfci[colSums(is.na(nfci))==nrow(nfci)] <- rowMeans(dframe[,-1],na.rm = TRUE)
for (i in 1:ncol(nfci)){
  nfci[is.na(nfci[,i]),i] <- global.nfci[is.na(nfci[,i])] * sd(nfci[,i],na.rm = TRUE)
}

write.csv(nfci,'../../data/clean/nfci.csv')
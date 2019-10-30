rm(list=ls())
library(dplyr)
library(lubridate)
library(readxl)
files <- dir('../../data/raw/epu/')
files <- files[grepl('.xlsx',files)]

# Set dataframe
dfout <- data.frame(read_excel('../../data/raw/epu/CAN.xlsx'))
dates <- paste0('01/',dfout[,2],'/',dfout[,1])
dates <- dates[1:(length(dates))]
dfout <- cbind.data.frame('date' = parse_date_time(as.character(dates),orders = 'd!/m!*/Y!')
                           ,'epu' = unname(dfout[1:(nrow(dfout)),3]) )

dfout <- data.frame(dfout %>% group_by(dates=floor_date(date, "quarter")) %>%
                          summarize(epu=mean(epu)))
cnames <- names(read.csv('../../data/raw/keep_countries.csv'))
df <- data.frame(matrix(NA,nrow=nrow(dfout),ncol=length(cnames)+1))
colnames(df) <- c('Dates',cnames)
df[,1] <- dfout$dates

for(j in 1:length(files)){
  file=files[j]
  d <- data.frame( read_excel(sprintf('../../data/raw/epu/%s',file)) )
  name <- gsub('.xlsx','',file)
  if (ncol(d) == 3){
    dates <- paste0('01/',d[,2],'/',d[,1])
    dates <- dates[1:(length(dates))]
    
    dframe <- cbind.data.frame('date' = parse_date_time(as.character(dates),orders = 'd!/m!*/Y!')
                               ,'epu' = unname(d[1:(nrow(d)),3]) )
    
    quarterly <- data.frame(dframe %>% group_by(dates=floor_date(date, "quarter")) %>%
                              summarize(epu=mean(epu)))
  } else if (ncol(d) == 2){
    dates <- d[,1]
    
    dframe <- cbind.data.frame('date' = parse_date_time(as.character(dates),orders = 'Y!/m!*/d!')
                               ,'epu' = unname(d[1:(nrow(d)),2]) )
    
    quarterly <- data.frame(dframe %>% group_by(dates=floor_date(date, "quarter")) %>%
                              summarize(epu=mean(epu)))
    }
  if (sum((df[,1] %in% quarterly$dates )) ==0 ) cat(sprintf('error at %s',name))
  df[ which(df[,1] %in% quarterly$dates ), 1 + which(cnames==name)] <- quarterly$epu[which(quarterly$dates %in% df[,1] )]
}

global <- scale( rowMeans( scale( df[,-1] ) ,na.rm = TRUE) )
for (i in 2:ncol(df)){
  if (sum(is.na(df[,i])) == nrow(df) ) {
    df[is.na(df[,i]),i] <-global[is.na(df[,i])]
  }else{
    df[is.na(df[,i]),i] <- global[is.na(df[,i])] * sd( df[,i],na.rm = TRUE )
  }
} 

df <- df[rowSums(is.na(df))==0,]
write.csv(df,sprintf('../../data/clean/epu.csv'))

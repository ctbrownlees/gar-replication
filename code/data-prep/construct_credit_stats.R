rm(list=ls())
library(lubridate)
library(countrycode)
library(readxl)

d     <- read_xlsx('../../data/raw/credit_stats.xlsx',sheet = 3)
cr    <- d[ , c(1, grep("ratios" , colnames(d))) ]
cg    <- d[ , c(1, grep("gaps" , colnames(d)))]
countries <- countrycode(as.character(unname(cr[2,])),'country.name','iso3c')
drop      <- which(is.na(countries))
countries <- countries[-drop]

cr   <- cr[4:nrow(cr),-drop[2]]
cg   <- cg[4:nrow(cg),-drop[2]]
colnames(cr) <- colnames(cg) <- c('Dates',countries)
cr[,1] <- as.Date(as.double(cr$Dates),origin="1899-12-30")
cg[,1] <- as.Date(as.double(cg$Dates),origin="1899-12-30")

cnames <- names(read.csv('../../data/raw/keep_countries.csv'))

cg.matched <- data.frame(cg[,colnames(cg) %in% c('Dates',cnames)])
cg.matched[,2:ncol(cg.matched)]   <- sapply(2:ncol(cg.matched),function(w) as.numeric(cg.matched[,w]))
cg.diff  <- sapply(2:ncol(cg.matched),function(w) diff((cg.matched[,w]),na.rm = T) )
cg <- cbind.data.frame('Dates'=cg.matched[2:nrow(cg.matched),1],cg.diff)
colnames(cg) <- colnames(cg.matched)
cg.global        <- scale( rowMeans(scale( cg.diff ),na.rm = TRUE) )

cr.matched <- data.frame(cr[,colnames(cr) %in% c('Dates',cnames)])
cr.matched[,2:ncol(cr.matched)]   <- sapply(2:ncol(cr.matched),function(w) as.numeric(cr.matched[,w]))
cr.diff  <- sapply(2:ncol(cr.matched),function(w) diff(log(cr.matched[,w]),na.rm = T)*100 )
cr <- cbind.data.frame('Dates'= cr.matched[2:nrow(cr.matched),1],cr.diff)
colnames(cr)    <- colnames(cr.matched)
cr.global       <- scale( rowMeans( scale( cr.diff ),na.rm = TRUE) )

cr.out <- data.frame(matrix(NA,nrow=nrow(cr),ncol=ncol(cr)+1))
colnames(cr.out) <- c('Dates',cnames)
for (name in colnames(cr.out)){
  if (name %in% colnames(cr))
  cr.out[ , name] <- cr[ ,name ]
}

cr.out$ISL   <- cr.global
cr.out$Dates <- cr.out$Dates + days(1)

cg.out <- data.frame(matrix(NA,nrow=nrow(cg),ncol=ncol(cg)+1))
colnames(cg.out) <- c('Dates',cnames)
for (name in colnames(cg.out)){
  if (name %in% colnames(cg))
    cg.out[ , name] <- cg[ ,name ]
}
cg.out$ISL <- cg.global
cg.out$Dates <- cg.out$Dates + days(1)

for (i in 2:ncol(cg.out)){
  cg.out[is.na(cg.out[,i]),i] <- cg.global[is.na(cg.out[,i])] * sd(cg.out[,i] , na.rm=T)
  cr.out[is.na(cr.out[,i]),i] <- cr.global[is.na(cr.out[,i])] * sd(cr.out[,i] , na.rm=T)
}
# Keeping only non NAs
cg.out <- cg.out[rowSums(is.na(cg.out))==0,]
cr.out <- cr.out[rowSums(is.na(cr.out))==0,]

write.csv(cg.out,'../../data/clean/cg.csv')
write.csv(cr.out,'../../data/clean/cr.csv')

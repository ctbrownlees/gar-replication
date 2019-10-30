rm(list=ls())
library(xtable)
########## QR SCREENING ##############
varlist <- c('nfci','ts','gf','cs','hp','sv','epu','cg','cr','wui','gpr')
vartex <- list()
for (vname in 1:length(varlist)){
  d <-  read.csv(sprintf('raw/insample/%s.csv',varlist[vname]))
  if( vname != 1){
    d[2,3:14] <- d[2,3:14]*100 # Tick loss as percentages
  }
  underset <- list()
  for (i in 2+seq(3,14,3)){
    underset[[i%/%3]] <- c(sprintf("$\\underset{[%.2f  \\, %.2f]}{%.2f}$",d[,i-2],d[,i],d[,i-1])[1:2],
                           sprintf("$\\underset{[%.2f  \\, %.2f]}{%.2f}$",d[,i-2],d[,i],d[,i-1])[3])
  }
  dtab <- d[,c(1,2,2+seq(2,12,3))]
  for(i in 1:(ncol(dtab)-2)){
    dtab[,2+i] <- underset[[i]]
  }
  if (vname==1) {
    tmp              <- c(sprintf("\\multirow{4}{*}{%s} &
                                  \\multirow{4}{*}{\\shortstack{%s \\\\ %s}}
                                ",
                                  toupper(varlist[vname]), strsplit(as.character(dtab[2,2]),"to")[[1]][1],strsplit(as.character(dtab[2,2]),"to")[[1]][2] ),
                          
                          sprintf("& %s & %s &  %s & %s & %s \\\\ ",
                                  c('\\multirow{2}{*}{Coef.}','&$TL$','& Sig.'),dtab[,3],dtab[,4],dtab[,5],dtab[,6]))
    
  } else {
    tmp              <- c(sprintf("\\multirow{4}{*}{%s} &
                                  \\multirow{4}{*}{\\shortstack{%s \\\\ %s}}
                                ",
                                  toupper(varlist[vname]), strsplit(as.character(dtab[2,2]),"to")[[1]][1],strsplit(as.character(dtab[2,2]),"to")[[1]][2] ),
                          
                          sprintf("& %s & %s &  %s & %s & %s \\\\ ",
                                  c('\\multirow{2}{*}{Coef.}','& $\\Delta TL$','& Sig.'),dtab[,3],dtab[,4],dtab[,5],dtab[,6]))
    
  }
  tmp <- c(tmp[1],tmp[2],tmp[4],tmp[3],"[1em]")
  
  vartex[[vname]]  <- tmp
}

### Header and footer 
t1 <- xtable(dtab)
lines <- print(t1, include.rownames = FALSE , sanitize.text.function = identity,booktab=TRUE)
lines <- strsplit(lines,'\n')[[1]]
lines[5] <- "\\begin{tabular}{lcccccc}" 
lines[7] <- " & Estimation  & & \\multicolumn{4}{c}{$h$}  \\\\ \\cmidrule(lr){4-7} & Window & & 1 & 2 & 3 & 4 \\\\
 \\cmidrule(lr){2-2}  \\cmidrule(lr){4-4}
\\cmidrule(lr){5-5} \\cmidrule(lr){6-6} \\cmidrule(lr){7-7} 
 "
header   <- lines[5:7]
footer   <- lines[12:13]

final.table <- c(header,unlist(vartex),footer)
final.table <- gsub("$\\underset{[0.000  \\, 0.000]}{0.000}$","$-$",final.table,fixed = TRUE)
final.table <- gsub("$\\underset{[0.00  \\, 0.00]}{-0.00}$","$0.00$",final.table,fixed = TRUE)
final.table <- gsub("\\underset{[0.000  \\, 0.000]}","",final.table,fixed = TRUE)
final.table <- gsub("\\underset{[0.00  \\, 0.00]}","",final.table,fixed = TRUE)
final.table <- gsub("HOUSEPRICES","HP",final.table,fixed=TRUE)
final.table <- gsub("CREDITGDP","CR",final.table,fixed=TRUE)
final.table <- gsub("CREDITGAP","CG",final.table,fixed=TRUE)
final.table <- gsub("FACTOR","GF",final.table,fixed=TRUE)
final.table <- gsub("SPREAD","TS",final.table,fixed=TRUE)
final.table <- gsub("BAATS","CS",final.table,fixed=TRUE)
final.table <- gsub("RV","SV",final.table,fixed=TRUE)


save.path <- 'tex'
name     <- 'qr_screen'
fileConn <- sprintf('%s/%s.tex',save.path,name)
write(final.table,fileConn)
############### QR MODELS  ##################
d <-  read.csv('raw/insample/qr_models.csv')
underset <- list()
for (i in seq(3,nrow(d),3)){
  if ( (i%/%3)%%2 ){
    underset[[i%/%3]] <- sprintf("$\\underset{[%.2f  \\, %.2f]}{%.2f}$",d[i-2,-1],d[i,-1],d[i-1,-1]) 
  } else {
    underset[[i%/%3]] <- sprintf("${%.2f}$",d[i-1,-1]*100)
  }
  
}
dtab <- rbind.data.frame( d[seq(2,nrow(d)-1,3),])
for(i in 1:nrow(dtab)){
  dtab[i,-1] <- underset[[i]]
}
dtab  <- rbind.data.frame(dtab,c('TL',sprintf('%0.4f',d[nrow(d)-1,-1])),c('DQ',sprintf('%0.2f',d[nrow(d),-1])))
t1 <- xtable(dtab)
lines <- print(t1, include.rownames = FALSE , sanitize.text.function = identity,booktabs = TRUE)
lines <- strsplit(lines,'\n')[[1]]
lines[5]  <- "\\begin{tabular}{lccccccccccccccccc}" 
lines[7]  <- "\\multicolumn{1}{c}{$h$} & & \\multicolumn{4}{c}{$1$} & \\multicolumn{4}{c}{$2$} & \\multicolumn{4}{c}{$3$}
              & \\multicolumn{4}{c}{$4$} \\\\
              \\cmidrule(lr){1-1}  \\cmidrule(lr){3-6} \\cmidrule(l){7-10} \\cmidrule(l){11-14} \\cmidrule(l){15-18}
              & & (1) & (2) & (3) & (4)
              & (1) & (2) & (3) & (4)
              & (1) & (2) & (3) & (4)
              & (1) & (2) & (3) & (4) \\\\ " 
lines[8]  <- '\\cmidrule(lr){3-6} \\cmidrule(l){7-10} \\cmidrule(l){11-14} \\cmidrule(l){15-18}'

lines[9]  <- gsub("AR_2","\\\\multirow{3}{*}{$Y_{it}$} &\\\\multirow{2}{*}{Coef.}",lines[9])
lines[9]  <- paste(lines[9],'[0.8em]')
lines[10] <- gsub("TAR_2","& Sig.",lines[10])
lines[10]  <- paste(lines[10],'[1em]')

lines[11] <- gsub("NFCI_2","\\\\multirow{3}{*}{NFCI} &\\\\multirow{2}{*}{Coef.}",lines[11])
lines[11]  <- paste(lines[11],'[0.8em]')
lines[12] <- gsub("TNFCI_2","& Sig.",lines[12])
lines[12]  <- paste(lines[12],'[1em]')

lines[13] <- gsub('TSP_2','\\\\multirow{3}{*}{TS} &\\\\multirow{2}{*}{Coef.}',lines[13])
lines[13]  <- paste(lines[13],'[0.8em]')

lines[14] <- gsub("TTSP_2","& Sig.",lines[14])
lines[14]  <- paste(lines[14],'[1em]')

lines[15] <- gsub("FAC_2","\\\\multirow{3}{*}{GF} &\\\\multirow{2}{*}{Coef.}",lines[15])
lines[15]  <- paste(lines[15],'[0.8em]')

lines[16] <- gsub("TFAC_2","& Sig.",lines[16])
lines[16]  <- paste(lines[16],'[1em]')

lines[17] <- gsub('CRV_2','\\\\multirow{3}{*}{SV} &\\\\multirow{2}{*}{Coef.}',lines[17])
lines[17]  <- paste(lines[17],'[0.8em]')

lines[18] <- gsub("TRV_2","& Sig.",lines[18])
lines[18]  <- paste(lines[18],'[1em]')


lines[19] <- gsub('CHP_2','\\\\multirow{3}{*}{HP} &\\\\multirow{2}{*}{Coef.}',lines[19])
lines[19]  <- paste(lines[19],'[0.8em]')

lines[20] <- gsub("TCHP_2","& Sig.",lines[20])
lines[20]  <- paste(lines[20],'[1em]')

lines[20] <- paste(lines[20],
                   '\\cmidrule(lr){3-6} \\cmidrule(l){7-10} \\cmidrule(l){11-14} \\cmidrule(l){15-18}')

lines[21] <- gsub('TL','TL &',lines[21])
lines[22] <- gsub('DQ','DQ Hits & ',lines[22])

final.table <- lines[5:24]
final.table <- gsub("$\\underset{[0.00  \\, 0.00]}{0.00}$","\\multirow{3}{*}{$-$}",final.table,fixed = TRUE)
final.table <- gsub("NaN","",final.table,fixed = TRUE)
save.path <- 'tex'
name     <- 'qr_models'
fileConn <- sprintf('%s/%s.tex',save.path,name)
write(final.table,fileConn)

############### GARCH IS   ##################
body <- list()
d <- read.csv('raw/insample/garch_iterated_is.csv')

# GARCH
under_setG <- sprintf("$\\underset{[%.3f  \\, %.3f]}{%.3f}$",d$GARCH025,d$GARCH075,d$GARCH050)
under_setG <- c( under_setG[1:9] , sprintf("$\\underset{[%.4f  \\, %.4f]}{%.4f}$",d$GARCH025[10],d$GARCH075[10],d$GARCH050[10]))
# TARCH
under_setT <- sprintf("$\\underset{[%.3f  \\, %.3f]}{%.3f}$",d$TARCH025,d$TARCH075,d$TARCH050)
under_setT <-  c( under_setT[1:9] , sprintf("$\\underset{[%.4f  \\, %.4f]}{%.4f}$",d$TARCH025[10],d$TARCH075[10],d$TARCH050[10]))
# FGARCH
under_setFG <- sprintf("$\\underset{[%.3f  \\, %.3f]}{%.3f}$",d$FGARCH025,d$FGARCH075,d$FGARCH050)
under_setFG <- c( under_setFG[1:9] , sprintf("$\\underset{[%.4f  \\, %.4f]}{%.4f}$",d$FGARCH025[10],d$FGARCH075[10],d$FGARCH050[10]))
dmed <- d[,c(1,3,6,9)]
dmed[,2] <- c(under_setG,  sprintf('%.2f' , dmed[11:13 ,2]) )
dmed[,3] <- c(under_setT , sprintf('%.2f' , dmed[11:13 ,3]) )
dmed[,4] <- c(under_setFG, sprintf('%.2f' , dmed[11:13 ,4]) )
dmed     <- dmed[-2,]
t1 <- xtable(dmed)
# Header
lines <- print(t1, include.rownames = FALSE , sanitize.text.function = identity,booktabs = TRUE)
lines <- strsplit(lines,'\n')[[1]]
lines[5] <- "\\begin{tabular}{lccc}"
lines[7] <-  sprintf(' & GARCH & GJR-GARCH & F-GARCH  \\\\ ')
lines[8] <-  "\\cmidrule(lr){2-2} \\cmidrule(lr){3-3} \\cmidrule(lr){4-4}"

lines[9]  <- sprintf(' $ \\phi$ & %s & %s & %s  \\\\[1em]  ',dmed$TARCH050[1],dmed$TARCH050[1],dmed$TARCH050[1] )
lines[11] <- sprintf('$\\sigma^2_u$ & %s & %s & %s \\\\[1em] ', dmed$TARCH050[3],dmed$TARCH050[3],dmed$TARCH050[3] )
lines[10] <- gsub('Loadings',' $\\\\lambda$',lines[10])
lines[10] <- paste(lines[10],'[1em]')
lines[12] <- gsub('Persistence',' Pers.',lines[12])
lines[12] <- paste(lines[12],'[1em]')

lines[13] <- gsub('Beta',' $\\\\beta$',lines[13])
lines[13] <- paste(lines[13],'[1em]')

lines[14] <- gsub('Gamma',' $\\\\gamma$',lines[14])
lines[14] <- paste(lines[14],"[1em]")

lines[15] <- gsub('Skewness','Skew.',lines[15])
lines[15] <- paste(lines[15],'[1em]')
lines[16] <- gsub('Kurtosis','Kurt.',lines[16])
lines[16] <- paste(lines[16],'\\cmidrule(lr){2-2} \\cmidrule(lr){3-3} \\cmidrule(lr){4-4}')

lines[17] <- gsub('TickLoss',' TL',lines[17])
lines[17] <- paste(lines[17],'[1em]')

lines[18] <- gsub("LR"," LR Test",lines[18])
lines[18] <- paste(lines[18],'[1em]')

lines[19] <- gsub("ARCH","ARCH-LM",lines[19])
lines[19] <- paste(lines[19],'[1em]')

lines[20] <- gsub("DQ"," DQ Hits",lines[20])

lines_rep_1 <- lines[17]
lines[17:18] <- lines[18:19]
lines[19]   <- lines_rep_1 

# ARCH - LM AND DQ HITS
header    <- lines[5:7]
body      <- lines[8:20]
tail      <- lines[21:22]

final.table <- c(header,unlist(body),tail)
final.table <- gsub("$\\underset{[0.000  \\, 0.000]}{0.000}$","$-$",final.table,fixed = TRUE)
final.table <- gsub("\\underset{[0.0000  \\, 0.0000]}","",final.table,fixed = TRUE)
final.table <- gsub("\\underset{[0.000  \\, 0.000]}","",final.table,fixed = TRUE)
save.path <- 'tex'
name     <- 'garch_table_iterated'
fileConn <- sprintf('%s/%s.tex',save.path,name)
write(final.table,fileConn)

############################ OOS TABLES
body <- list()
for (cover in c(0.95)){
  cov.lev <- gsub("\\.","",sprintf('%.2f',cover))
  # Uniform tables
  read.path = sprintf('raw',cov.lev)
  tables <- dir( read.path )
  tables.u <- tables[grepl('Uniform', tables)]
  save.path = 'tex'
  body <- list()
  for (t in 1:length(tables.u)){
    
    d <- read.csv(sprintf("%s/%s",read.path,tables.u[t]))
    dtmp <- d
  
      # Adjusting statistics
    d$DQ[d$DQ>=1-0.0001] <- 0.9990
    d$DQReal[d$DQReal>=1-0.0001] <- 0.9990
    d$DQFin[d$DQFin>=1-0.0001] <- 0.9990
    #
  
    name <- strsplit(tables.u[t],".csv")[[1]]
    l1    <- xtable(d, caption=sprintf("%s", gsub("_|\\d+_ahead"," ",name)  ))
    digits(l1) <- c(0,0,2,3,3,3,3,3)
    lines <- print( l1,include.rownames=FALSE,caption.placement = 'top',scalebox = 0.55,NA.string = '-',booktabs = TRUE)
    lines <- strsplit(lines,'\n')[[1]]
    lines[7]  <- "\\begin{tabular}{llp{10em}cccccc}"
    lines[9]  <- "{$h$} & {Method} & {Model} & {Cov.} &
  {Length} & {Unc.} & {Hits} & {NFCI} & {Real} \\\\
    \\cmidrule{1-1} \\cmidrule(lr){2-3}  \\cmidrule(lr){4-5} \\cmidrule(lr){6-9}"
    lines[10] <- ""
    lines[11] <- paste(sprintf('\\multirow{15}{*}{{%s}} &
                      \\multirow{1}{*}{{Benchmark}} &',strsplit(name,'_')[[1]][4])
                       ,gsub('Historical \\+ BJPR','Benchmark',lines[11]),
                       '\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-9}')
    ## QR
    lines[12] <- gsub('QR-NFCI','NFCI',lines[12])
    lines[13] <- gsub('QR-NFCI\\+TS','NFCI + TS',lines[13])
    lines[14] <- gsub('QR-NFCI\\+GF\\+TS','NFCI + TS + GF ',lines[14])
    
    lines[12] <- paste('& \\multirow{5}{*}{{QR + Bonf.}} &',lines[12])
    lines[13:16] <- paste('& &',lines[13:16])
    lines[16] <- gsub('Lasso','LASSO',lines[16]) 
    lines[16] <- paste(lines[16],
                       '\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-9}')
    ## Copy to move
    tmpM      <- lines[23:25]
    tmpM      <- gsub("\\+ Marg",'',tmpM)
    tmpBJPR   <- lines[17:19]
    tmpBJPR   <- gsub("\\+ BJPR",'',tmpBJPR)
    tmpBonf   <- lines[20:22] 
    tmpBonf   <- gsub("\\+ Bonf.",'',tmpBonf)
    ## Marg
    lines[17] <- paste('& \\multirow{3}{*}{{GARCH + Marg.}} & ',tmpM[1])
    lines[18:19] <- paste('& &',tmpM[2:3])
    lines[19] <- paste(lines[19],
                       '\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-9}')
    ## Bonf
    lines[20] <- paste('& \\multirow{3}{*}{{GARCH + Bonf.}} & ',tmpBonf[1])
    lines[21:22] <- paste('& &',tmpBonf[2:3])
    lines[22] <- paste(lines[22],
                       '\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-9}')
    ## BJPR
    lines[23] <- paste('& \\multirow{3}{*}{{GARCH + BJPR}} & ',tmpBJPR[1])
    lines[24:25] <- paste('& &',tmpBJPR[2:3])  
    lines[26] <- "\\toprule"
    body[[t]]<- lines[11:26]
  }
  name  <- sprintf('uniform_tables_%s',cov.lev)
  header <- lines[7:10]
  tail   <- lines[27]
  final.table <- c(header , unlist(body), tail)
  fileConn <- sprintf('%s/%s.tex',save.path,name)
  write(final.table,fileConn)
  
  # Marginal tables
  tables.m <- tables[grepl('Marginal', tables)]
  body <- list()
  for (t in 1:length(tables.m)){
    
    d <- read.csv(sprintf("%s/%s",read.path,tables.m[t]))
    name <- strsplit(tables.m[t],".csv")[[1]]
    dtmp <- d
    d$Tloss <- (1-d$Tloss/d$Tloss[1])*100
    dtmp$Tloss[-1] <- d$Tloss[-1]
    
    tloss.min  <- which(d$Tloss == max(d$Tloss))
    d$Tloss[tloss.min]    <- 1007.11
    d$Tloss[1]    <- dtmp$Tloss[1]
    
    l1    <- xtable(d, caption=sprintf("%s", gsub("_|\\d+_ahead"," ",name)  ))
    dmat    <- matrix( rep( c(0,0,2,3,2,2,2,2,2), 9) ,nrow = 9 ,byrow = TRUE)
    dmat[1,9] <- 4
    digits(l1) <- dmat
    lines <- print( l1,include.rownames=FALSE,caption.placement = 'top',scalebox = 0.6,booktabs = TRUE)
    lines <- strsplit(lines,'\n')[[1]]
    lines[c(1:6,25,26)] <- "" 
    lines[7]  <- "\\begin{tabular}{llp{10em}ccccccc}" 
    lines[9]  <- "{$h$} & {Method} & {Model} & 
  {Cov.} & {Length} & {Unc.}  & {Hits} & {NFCI} 
  & {Real}& {TL} \\\\
  \\cmidrule{1-1} \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}
  \\cmidrule(lr){6-9} \\cmidrule(lr){10-10}
    "
    lines[10] <- ""
    lines[11] <- paste(sprintf('\\multirow{9}{*}{{%s}} &
                      \\multirow{1}{*}{{Benchmark}} &',strsplit(name,'_')[[1]][4])
                       ,lines[11],
                       '\\cmidrule(lr){2-3}\\cmidrule(lr){4-5} \\cmidrule(lr){6-9} \\cmidrule(lr){10-10}')
    lines[12:13] <- paste('& &',lines[12:13])
    lines[15:16] <- paste('& &',lines[15:16])
    lines[14]    <- paste('& \\multirow{1}{*}{{QR}} & ' , lines[14])
    lines[16] <- gsub('Lasso','LASSO',lines[16]) 
    lines[16] <- paste((lines[16]),
                       '\\cmidrule(lr){2-3}\\cmidrule(lr){4-5} \\cmidrule(lr){6-9} \\cmidrule(lr){10-10}') 
    lines[17] <- paste('& \\multirow{3}{*}{{GARCH}} & ',lines[17])
    lines[18:19] <- paste('& &',lines[18:19])
    lines[19]    <- paste(lines[19],'\\toprule')
    
    lines[12] <- gsub('QR-NFCI','NFCI',lines[12])
    lines[13] <- gsub('QR-NFCI\\+TS','NFCI + TS',lines[13])
    lines[14] <- gsub('QR-NFCI\\+TS\\+GF','NFCI + TS + GF ',lines[14])
    
    header <- lines[7:10]
    body[[t]] <- lines[11:19]
    
    body[[t]]<- gsub("1007.11",sprintf("\\\\textbf{%.2f}",as.numeric(dtmp$Tloss[tloss.min])),body[[t]])
    tail     <-  lines[21]
    final.table <- c(header , unlist(body),tail)
  }
  name  <- sprintf('marginal_tables_%s',cov.lev)
  fileConn <- sprintf('%s/%s.tex',save.path,name)
  write(final.table,fileConn)
  
  ########## Make IMF tables #############
  tables.imf <- tables[grepl('IMF', tables)]
  body <-dtmp <-  list()
  for (t in 1:length(tables.imf)){
    d <- read.csv(sprintf("%s/%s",read.path,tables.imf[t]))
    name <- strsplit(tables[t],".csv")[[1]]
    dtmp[[t]] <- d
    d    <- d[-c(6:8,13:15,20:22)]
    cov <- cov.lev
    
    l1    <- xtable(d, caption=sprintf("%s", gsub("_|\\d+_ahead"," ",name)  ))
    digits(l1) <- c(0,0,rep(c(2,3,3,4),3))
    lines <- print( l1,include.rownames=FALSE,caption.placement = 'top',scalebox = 0.6,booktabs=TRUE)
    lines <- strsplit(lines,'\n')[[1]]
    lines[c(1:6,25,26)] <- "" 
    lines[7]  <- "\\begin{tabular}{llcccccccccccc}"
    lines[9]  <- "\\multicolumn{1}{c}{\\multirow{2}{*}{$h$}} & \\multicolumn{1}{c}{\\multirow{2}{*}{Country}} &
  \\multicolumn{4}{c}{{Historical}} &   \\multicolumn{4}{c}{{QR NFCI}} & 
  \\multicolumn{4}{c}{{GARCH}} \\\\ "
    lines[9] <- paste(lines[9],
                      "& &     Cov. & Length & Unc. & TL &
                               Cov. & Length & Unc. & TL &
                               Cov. & Length & Unc. & TL \\\\ ")
    lines[10] <- " \\cmidrule{1-1}  \\cmidrule(lr){2-2} \\cmidrule(lr){3-6} \\cmidrule(lr){7-10} \\cmidrule(lr){11-14}" 

    lines[11] <- paste(sprintf('\\multirow{12}{*}{{%s}}
                       &',strsplit(name,'_')[[1]][4])
                       ,lines[11])
    lines[12:22] <- paste('& ',lines[12:22])
    lines[23]    <- "\\toprule"
    
    header <- lines[1:10]
    body[[t]] <- lines[11:23]
    tail   <- lines[24:26]
    
    final.table <- c(header , unlist(body), tail)
  }
  
  boldface <- function(tb){
    nr <- length(tb)
    for(i in 1:nr){
      obj <-strsplit(tb[i],"&")[[1]] 
      if (length(obj)==14){
        parsed_obj  <- as.numeric(unlist(strsplit(gsub("[^[:digit:].]", "", obj), " +")))
        if(length(parsed_obj>12)) parsed_obj <- tail(parsed_obj,12)
        #to_replace  <- sprintf("%.4f",c( parsed_obj[c(1,5,9)][which.min(abs(parsed_obj[c(1,5,9)]-cov))] ,
        #                 min(parsed_obj[2:6]) , max(parsed_obj[7:9]) , min(parsed_obj[10:12])))
        to_replace  <- sprintf("%.4f",c(min(parsed_obj[c(4,8,12)])))
        for (w in 1:length(to_replace)){
          obj <- gsub(to_replace[w],sprintf("\\\\textbf{%s}",to_replace[w]),obj)
        }
        tbindex <- seq(10,62,13)
        for (t in 2:length(tbindex)){
          if (i > tbindex[t-1]  & i<=tbindex[t]){
            # QR stars
            # How many stars
            nstar <- sum(dtmp[[t-1]][i-tbindex[t-1],13:15])
            if( dtmp[[t-1]][i-tbindex[t-1],15] ) obj[10] <- sprintf('$%s^{%s}$',obj[10],paste(rep('*',nstar),collapse=''))
            # GARCH stars
            nstar <- sum(dtmp[[t-1]][i-tbindex[t-1],20:22])
            if( dtmp[[t-1]][i-tbindex[t-1],ncol(dtmp[[t-1]])] ) obj[14] <- sprintf('$%s^{%s}$ \\\\',strsplit(obj[14]," ")[[1]][2],
                                                                                   paste(rep('*',nstar),collapse=''))  
          } 
        }
        
        tb[i] <- paste(obj[1],paste("&",obj[2:length(obj)],collapse = " "))
      } else{
        next
      }
    }
    tb
  }
  
  final.table <- boldface(final.table)
  name  <- sprintf('IMF_table_%s',cov.lev)
  fileConn <- sprintf('%s/%s.tex',save.path,name)
  write(final.table,fileConn)
}
  

#' @title SCRAPI v1.0
#'
#' @description The OUTDATED juvenile companion to SCOBI
#'
#' @param smoltData the name of the file containing the smolt biological data
#' @param Dat the column in smoltData containing the sampling date
#' @param Rr column heading containing code for rearing type
#' @param Primary This is the primary classification variable for grouping
#' @param Secondary Use NA if no secondary composition factor is requested
#' @param passageData the name of the file containing the smolt passage data
#' @param strat This is the column name for stratification (original weeks)
#' @param dat Column heading containing sampling date
#' @param tally Column heading containing counts at the dam
#' @param samrate Column heading containing trap sampling rate
#' @param guidance Column heading containing value of guidance efficiency
#' @param collaps column containing the collapsing scheme
#' @param Run synopsis of the run being conducted. This will be used as the prefix for all of your output files
#' @param RTYPE which rear type (wild [H], hatchery [H], or hatchery unclipped [HNC]) would you like to perform analysis on
#' @param REARSTRAT would you like pWild calculated on a time-stratified basis (TRUE) or across the entire emigration (FALSE)
#' @param alph the desired alpha level to calculate confidence intervals
#' @param B the desired number of bootstraps
#'
#' @author Kirk Steinhorst and Mike Ackerman
#'
#' @import stats plyr
#' @return NULL

SCRAPIv1.0 <- function(smoltData = NULL, Dat = "CollectionDate", Rr = "Rear", Primary = "GenStock", Secondary = NA, passageData = NULL,
                       strat = "Week", dat = "SampleEndDate", tally = "SampleCount", samrate = "SampleRate", guidance = "GuidanceEfficiency",
                       collaps = "Collapse", Run = "output", RTYPE = "W", REARSTRAT = TRUE, alph = 0.1, B = 5000)
{

# Import data
All  <- read.csv(file = smoltData, header = TRUE)
pass <- read.csv(file = passageData, header = TRUE)

# Write header
cat("\nStart time: ",date(),"\n")
cat("\nThis is a run of ", Run, "\n")
cat("\nFocus is on fish of type ",RTYPE,"\n")
cat("\nPrimary composition variable is ",Primary," and the secondary variable is ",Secondary, "\n")
cat("\nNumber of bootstrap iterations: B = ",B,"\n")
cat("\nAlpha = ",alph,"\n")

# Function to compute the overall average composition for SECONDARY at each level of PRIMARY
getAvgProp <- function(Fh){
  Fh <- droplevels(Fh[which(Fh$SGrp!="NA"),])
  Freqs <- mApply(1/Fh$SR,list(factor(Fh$PGrp,levels=Pgrps),factor(Fh$SGrp,levels=Sgrps)),sum)
  Freqs[is.na(Freqs)] <- 0
  AP <- prop.table(Freqs,margin=2)
  return( AP )
}

# The following function computes estimates of total smolt passage and passage by category
# Fish should have collapsed week, PRIMARY category, SECONDARY category (if needed), and
# realized sampling rate (true sampling rate x sampling proportion for PRIMARY category) for each fish.
# passage <- data.frame(Stratum=pass[,PASScollaps],Tally=pass[,PASScounts],Ptrue=pass$true)
thetahat <- function(Fish,passage){
  dailypass <- passage$Tally/passage$Ptrue
  bystrata <- mApply(dailypass,passage$Stratum,sum)
  bystrata <- ProWild*bystrata # This adjusts the (statistical) weekly passage to just wild smolts
  total <- sum(bystrata)
  Primarystrata <- mApply(1/Fish$SR,list(Fish$Week,Fish$PGrp),sum)
  Primarystrata[is.na(Primarystrata)] <- 0
  Primaryproportions <- prop.table(Primarystrata,margin=2)
  Primaryests <- Primaryproportions%*%bystrata
  if( !is.na(Secondary) ){
    SecondAbund <- array(numeric(nPgrps*nSgrps*nweaks),dim=c(nPgrps,nSgrps,nweaks))
# Delete fish with NA for SECONDARY categorical variable
    Fish <- droplevels(Fish[which(Fish$SGrp!="NA"),])
    Freqs <- mApply(1/Fish$SR,list(factor(Fish$PGrp,levels=Pgrps),factor(Fish$SGrp,levels=Sgrps),factor(Fish$Week,levels=weaks)),sum)
    Freqs[is.na(Freqs)] <- 0
    Props <- prop.table(Freqs,margin=c(1,3))
# Check to see if the proportions for the Secondary category are NaNs because the corresponding Freqs are all 0.
# If so, put in the average proportions.
    if( any(is.nan(Props)) ) {
      for( h in 1:nweaks ) {
        for( i in 1:nPgrps ) {
        if( any(is.nan(Props[i,,h])) ) Props[i,,h] <- AvgProp[,i]
        }
      }
    }
    for( h in 1:nweaks ) {
      ThisPrime <- bystrata[h]*Primaryproportions[,h]
      SecondAbund[,,h] <- as.vector(t( c(ThisPrime * Props[,,h] )))
    }
    PrimeBySecond <- apply(SecondAbund,c(1,2),sum)
    PrimeBySecond[is.na(PrimeBySecond)] <- 0
    assign("sTable", PrimeBySecond, envir = .GlobalEnv)
    Second <- apply(PrimeBySecond,2,sum)
    return( list(total,bystrata,Primaryproportions,Primaryests,Second,PrimeBySecond) )
  } else   return( list(total,bystrata,Primaryproportions,Primaryests) )

} # end of thetahat function

bootsmolt <- function(FishDat,LGDdaily){

# Set up storage for bootstrap results
  theta.b <- matrix(numeric(p*B),ncol=p)
  if( !is.na(Secondary) )
    theta.b[1,] <- c(ests[[1]],t(ests[[4]]),ests[[5]],as.vector(t(ests[[6]])))  # First row contains estimates from original data
    else  theta.b[1,] <- c(ests[[1]],t(ests[[4]]))  # First row contains estimates from original data
  dailypass <- round(LGDdaily$Tally/LGDdaily$Ptrue)

# bootstrap loop
  for (b in 2:B) {
# Find bootstrap values for number of smolts passing LGR each day
    cntstar = numeric(ndays)
    for (i in 1:ndays) {
      if(dailypass[i] != 0 ) cntstar[i] <- rbinom(1,dailypass[i],LGDdaily$Ptrue[i])
    } # end of for loop on daily smolts
    dailyStar <- data.frame(Stratum=LGDdaily$Stratum,Tally=cntstar,Ptrue=LGDdaily$Ptrue)

# Now find a bootstrap sample of AllData (i.e. individual fish) BY statistical WEEK using a weighted bootstrap
     H <- 0
     for ( h in weaks) {
        justwk <- FishDat[FishDat$Week==h,]    # grab All smolts of a given statistical week
        nwk <- nrow(justwk)
        i <- sample.int(nwk,replace=TRUE,prob=unlist(justwk$SR))
        wkstar <- justwk[i,]
        if ( H == 0 ) {indivStar <- wkstar
                       H <- 1
        } else { indivStar <- rbind(indivStar,wkstar)  }
     }
     eststar <- thetahat(indivStar,dailyStar)
     if( !is.na(Secondary) )
       theta.b[b,] <- c(eststar[[1]],t(eststar[[4]]),eststar[[5]],as.vector(t(eststar[[6]])))
       else  theta.b[b,] <- c(eststar[[1]],t(eststar[[4]]))

  } # end of bootstrap loop

# Find confidence intervals for each statistic
  CI <- matrix(numeric((ncol(theta.b)*3)),ncol=3)
  for  (j in 1:p) {
      CIj <- quantile(theta.b[,j],c(alph/2,1-alph/2))
      CI[j,] <- c(theta.b[1,j],CIj)
  }
  CI <- round(CI)

# Find simultaneous confidence intervals for smolt composition via Mandel/Betensky 2008
# MBCI <- SCSrank(theta.b[,-1],conf.level=1-alph)
  return( CI)
} # end of bootsmolt

#------------  MAIN ---------------------------------------------

# Set column numbers for fish data
# FISHstrat <- which(Strat == names(All))
FISHdate <- which(Dat == names(All))
FISHrear <- which(Rr == names(All))
FISHpndx <- which(Primary == names(All))
if( !is.na(Secondary) ) { FISHsndx <- which(Secondary == names(All))
    FISHsgrps <- unique(All[,FISHsndx])
    nSgrps <- length(FISHsgrps) }

# Set column numbers for passage data
PASSstrat   <- which(strat == names(pass))
strata      <- unique(pass[,PASSstrat])
PASSdate    <- which(dat == names(pass))
PASSrate    <- which(samrate == names(pass))
PASScounts  <- which(tally == names(pass))
PASSguideff <- which(guidance == names(pass))
PASScollaps <- which(collaps == names(pass))

# Number of days for which there are passage counts
ndays <- nrow(pass)

# Capture collapsing pattern
Cpattern <- unique(cbind(pass[,PASSstrat],pass[,PASScollaps]))
cat("\nStrata are collapsed according to: \n")
temp <- t(Cpattern)
rownames(temp) <- c("Week","Strata")
colnames(temp) <- rep("",ncol(temp))
print(temp)
weeks = Cpattern[,1]
nwks <- length(weeks) # Number of weeks where passage counts are made

# Calculations for passage data
pass$true      <- pass[,PASSrate]*pass[,PASSguideff]  # "True" sample rate = Fred's rate x guidance efficiency
pass$estimated <- pass[,PASScounts]/pass$true # This is the estimate of fish passing on each day.
totalpassage   <- round(sum(pass$estimated))
cat("\nEstimate of total smolts: ",totalpassage,"\n")

# passage should have daily values for collapsed week, smolt counts and true sampling rate (trap rate x guidance efficiency).
passdata <- data.frame(Stratum=pass[,PASScollaps],Tally=pass[,PASScounts],Ptrue=pass$true)

# Calculate smolt passage estimates by original weeks
passweeks  <- mApply(pass$estimated,pass[,PASSstrat],sum)
rpassweeks <- round(passweeks)
cat("\nEstimate of total smolts by week: \n")
print(rpassweeks)

# Calculate smolt passage estimates by collapsed weeks
passcollaps  <- mApply(pass$estimated,pass[,PASScollaps],sum)
rpasscollaps <- round(passcollaps)
cat("\nEstimate of total smolts by statistical week: \n")
print(rpasscollaps)

# Sort out the fish composition data

nAll <- nrow(All)
# Create All$Collaps by matching dates in All and pass
All$Collaps <- numeric(nAll)
AllDates    <- unique(All[,FISHdate])
for( d in AllDates ) {
  CollStrat <- pass[which(pass[,PASSdate] == d),PASScollaps]
  All$Collaps[which(All[,FISHdate] == d)] <- CollStrat
}

# Calculate the proportion of smolts that are wild and hatchery
Tabl        <- table(All[,FISHrear],All$Collaps)
PropWild    <- prop.table(Tabl,margin=2)
PercentWild <- 100*round(PropWild,3)
ColTotals   <- apply(Tabl,1,sum)
cat("\nTrapped hatchery and wild by statistical week\n")
print( cbind(Tabl,ColTotals) )
Proportions <- ColTotals/sum(ColTotals)
Ovrallpct <- 100*round(Proportions,3)
cat("\nPercent hatchery and wild by statistical week\n")
print( cbind(PercentWild,Ovrallpct) )
if(REARSTRAT == TRUE) ProWild  <- PropWild[2,] else ProWild <- Proportions[2]
WildCollaps <- ProWild*passcollaps
TotalWild   <- sum(WildCollaps)
cat("\nWild smolts by statistical week\n")
print( round(c(WildCollaps,TotalWild)) )

# Select those fish of RTYPE
AllRTYPE   <- droplevels(All[which(All[,FISHrear]==RTYPE),])
# Delete fish with NA for PRIMARY categorical variable
AllPrimary <- droplevels(AllRTYPE[which(AllRTYPE[,FISHpndx]!="NA"),])  # CHECK TO SEE IF THIS WORKS FOR CHARACTER and NUMERIC VARIABLES
nFISH <- nrow(AllPrimary)

# Find dates where smolt counts were made AND PRIMARY fish samples were taken.
set    = intersect(AllPrimary[,FISHdate],pass[,PASSdate])
ndates <- length(set)

Pgrps  <- sort(unique(AllPrimary[,FISHpndx]))
nPgrps <- length(Pgrps)
p      <- 1 + nPgrps
weaks  <- unique(Cpattern[,2])
nweaks <- length(weaks)

# Find number of fish sampled for the PRIMARY categorical variable by date
tabl   <- table(AllPrimary[,FISHdate],AllPrimary[,FISHpndx])
nPrime <- apply(tabl,1,sum)

# The realized sampling rate for a date is pass$true x proportion of fish sampled for the PRIMARY categorical variable on that date.
# Put the realized sampling rate on each smolt record, i.e.,
# ID smolts of a given date and attach the realized sampling rate to each.
#AllPrimary$SR = numeric(nFISH)
#for ( nn in 1:ndates ) {
#  passtemp <- pass[pass[,PASSdate]==as.Date(set[nn],origin="1970-01-01"),]
#  Primecount <- nPrime[which( as.Date(names(nPrime)) ==  as.Date(set[nn],origin="1970-01-01"))]
#  AllPrimary$SR[AllPrimary[,FISHdate] == as.Date(set[nn],origin="1970-01-01")] <- passtemp$true*Primecount/passtemp[,PASScounts]
#}
AllPrimary$SR = numeric(nFISH)
for ( nn in 1:ndates ) {
  passtemp <- pass[as.Date(pass[,PASSdate],format="%m/%d/%Y")==as.Date(set[nn],format="%m/%d/%Y"),]
  Primecount <- nPrime[which( as.Date(names(nPrime),format="%m/%d/%Y") ==  as.Date(set[nn],format="%m/%d/%Y"))]
  AllPrimary$SR[as.Date(AllPrimary[,FISHdate],format="%m/%d/%Y") == as.Date(set[nn],format="%m/%d/%Y")] <- passtemp$true*Primecount/passtemp[,PASScounts]
}

# Keep week, primary (and if needed secondary) group(s), and realized sampling rate
if( !is.na(Secondary) ){
      Allfish <- data.frame(Week=AllPrimary$Collaps,PGrp=AllPrimary[,FISHpndx],SGrp=AllPrimary[,FISHsndx],SR=AllPrimary$SR)
      Sgrps <- unique(as.character(Allfish$SGrp))
      Sgrps <- sort(Sgrps[Sgrps!="NA"])
      nSgrps <- length(Sgrps)
      p <- p + nSgrps + nPgrps*nSgrps
      AvgProp <- getAvgProp(Allfish)
      } else
      Allfish <- data.frame(Week=AllPrimary$Collaps,PGrp=AllPrimary[,FISHpndx],SR=AllPrimary$SR)

# Run thetahat() function
ests <- thetahat(Allfish,passdata)
pPropTable <- t(ests[[3]])
if(!is.na(Secondary)) sAbunTable <- ests[[6]]

# Run bootsmolt() function
answer <- bootsmolt(Allfish,passdata)
answer <- cbind(answer,round((((answer[,3]-answer[,2])/answer[,1])/2)*100, 1))
colnames(answer) <- c("Estimate","LCI","UCI","P1")

if( !is.na(Secondary) ){
# Now concatenate the primary and secondary names
      grpnams <- ""
      for( prim in Pgrps ) {
        for( nam in Sgrps ) grpnams = c(grpnams,paste0(nam,prim))
      }
      grpnams <- grpnams[-1]
      rownames(answer) <- c("WildSmolts",as.character(Pgrps),Sgrps,grpnams)
} else  rownames(answer) <- c("WildSmolts",as.character(Pgrps))
cat("\n")
print(answer)

# WRITE OUTPUTS
# Write Total Smolts, p(Wild), and Wild Smolts by Strata
tSmolts <- t(rbind(rpasscollaps,round(ProWild,4),round(WildCollaps,0)))
colnames(tSmolts) <- c("TotalSmolts","p(Wild)","WildSmolts")
write.csv(tSmolts, file = paste(Run,"Rear.csv",sep=""))#, append = FALSE, row.names = TRUE, col.names = TRUE, quote = FALSE)

# Write CIs Output
ciFile <- paste(Run,"CIs.csv",sep="")
header <- if(is.na(Secondary)) { paste(RTYPE,"-",Primary) } else { paste(RTYPE,"-",Primary,"-",Secondary) }
write.table(header, file = ciFile, append = FALSE, row.names = FALSE, col.names = FALSE, quote = FALSE)
suppressWarnings(write.table(answer, file = ciFile, col.names = NA, sep =",", append = TRUE))
write.table(paste("RearSampleSize =",nAll), file = ciFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(paste("PrimeSampleSize =",nFISH), file = ciFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
if(!is.na(Secondary)) { write.table(paste("SecondSampleSize =",sum(!is.na(Allfish$SGrp))), file = ciFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE) }

# Write Primary Results
primeFile   <- paste(Run,"Prime.csv",sep="")
pFreqTable  <- table(AllPrimary$Collaps,AllPrimary[,Primary])
pAbunTable  <- round(sweep(pPropTable, MARGIN = 1, WildCollaps, '*'),0)
PrimeTotals <- round(apply(pAbunTable,2,sum),0)
write.table(pFreqTable, file = primeFile, col.names = NA, sep = ",", append = FALSE)
write.table(pPropTable, file = primeFile, row.names = TRUE, col.names = FALSE, append = TRUE, sep = ",")
write.table(pAbunTable, file = primeFile, row.names = TRUE, col.names = FALSE, append = TRUE, sep = ",")
write.table(t(PrimeTotals), file = primeFile, row.names = "PrimeTotals", col.names = FALSE, append = TRUE, sep = ",")

# Write Primary x Secondary Results
if(!is.na(Secondary)) {
  sAbunTable <- cbind(sAbunTable,apply(sAbunTable,1,sum))
  sAbunTable <- rbind(sAbunTable,apply(sAbunTable,2,sum))
  rownames(sAbunTable) <- c(levels(Pgrps),"sTotals")
  colnames(sAbunTable) <- c(Sgrps,"pTotals")
  sAbunTable <- round(sAbunTable,0)
  write.table(sAbunTable, file = paste(Run,"PxS.csv",sep=""), col.names = NA, append = FALSE, sep = ",")
}

# END
cat("\nEnd time: ",date(),"\n")
}


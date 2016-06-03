#' @title SCOBI: Salmonid COmposition Boostrap Intervals (v2.0)
#'
#' @description Perform compositional analyses of adults at Lower Granite Dam.
#'
#' @param adultData the .csv file containing the biological data for the adults to be analyzed. The file should contain data for all fish trapped in
#' a given spawn year and for a given species. The function \code{lgr2SCOBI()} can be used to format data exported from the LGTrappingDB to make it
#' ready for the \code{SCOBI()} function
#' @param windowData the .csv file containing the window count data for the adults to be analyzed. \code{windowData} should contain 3 columns: Strata,
#' Count, and Collaps
#' @param Run synopsis of the run being conducted. \code{Run} will be used as the prefix for all of your output files. \code{Run} should generally
#' contain the spawn year, rear, species, primary, and (if desired) secondary categories. For example, "sy2015wSthdStockSex"
#' @param RTYPE which rear type (wild [W], hatchery [H], or hatchery unclipped [HNC]) would you like to perform the current analysis on
#' @param Primary the primary category in the \code{adultData} to be estimated. i.e., after \code{RTYPE} what category would you like to decompose next?
#' @param Secondary the secondary category in the \code{adultData} to be estimated. i.e. after \code{Primary} what category would you like to
#' decompose next? Use \code{Secondary = NA} if no secondary decomposition is desired
#' @param SizeCut if \code{Secondary} = "LGDFLmm", what FL would you like to use to separate Sm versus Lg fish. For steelhead, set \code{SizeCut = 780}
#' (the default) to separate A-run and B-run fish. For Chinook, set \code{SizeCut = 570} to separate Jacks and Adults
#' @param alph the alpha used for confidence intervals (e.g., \code{alph} = 0.10 results in 90 percent CIs)
#' @param B the number of bootstraps to be performed
#' @param writeThetas would you like to write the thetas to a .csv? These are the \code{B} bootstrap estimates of the Rear, Primary, and
#' Secondary (if set) categories.
#' @param writeOutput would you like to write output files?
#' @param pbtExpand should only be used is \code{RTYPE} is set to "H" or "HNC". If this is an analysis of PBT assigned fish, would you like to expand
#' the frequencies within each time strata by the PBT tag rates defines in \code{pbtRates}? H stock proportions within each strata will then be
#' calculated using the expanded frequencies.
#' @param pbtRates if \code{pbtExpand = TRUE, SCOBI()} expects a .csv file containing the PBT tag rates for any H stock that occurs in
#' the \code{adultData}. the \code{pbtRates} should contain 2 columns; the first column should contain the tag rates that correspond exactly
#' to the H stocks in the \code{adultData}, the second column should contain the PBT tag rate for each H stock
#'
#' @seealso \code{\link[MCPAN]{SCSrank}}
#' @author Kirk Steinhorst and Mike Ackerman
#'
#' @examples SCOBI(adultData = sthdScobiInput, windowData = sthdWindowCounts, Run = "sthdDemo",
#' RTYPE = "W", Primary = "GenStock", Secondary = "GenSex", alph = 0.1, B = 100)
#'
#' SCOBI(adultData = pbtTestFish, windowData = sthdWindowCounts, Run = "sthdPbtDemo",
#' RTYPE = "H", Primary = "PbtStock", Secondary = NA, alph = 0.1, B = 100)
#'
#' @references Steinhorst, K., T. Copeland, M. W. Ackerman, W. C. Schrader, E. C. Anderson (In review) Estimates and Confidence Intervals for Run Composition
#' of Returning Salmonids. Fishery Bulletin.
#'
#' @importFrom Hmisc mApply
#' @export
#' @return NULL

SCOBI <- function(adultData = NULL, windowData = NULL, Run = "output", RTYPE = "W", Primary = "GenStock",
                  Secondary = NA, SizeCut = 780, alph = 0.1, B = 5000, writeThetas = FALSE, writeOutput = TRUE,
                  pbtExpand = FALSE, pbtRates = NULL)
{
  # THIS IS SCOBIv2.0
  if ( is.character(adultData) == TRUE )  { Fishdata <- read.csv(file = adultData, header = TRUE, na.strings = c("NA","")) } else { Fishdata <- adultData }
  if ( is.character(windowData) == TRUE ) { Windata  <- read.csv(file = windowData, header = TRUE) } else { Windata <- windowData }
  if ( pbtExpand == TRUE ) {
    if ( is.character(pbtRates) == TRUE ) { pbtRate <- read.csv(file = pbtRates, header = TRUE) } else { pbtRate <- pbtRates }
    }
  collaps <- Windata$Collaps

# Write header
cat("\nStart time: ",date(),"\n")
cat("\nThis is a run of ", Run, "\n")
cat("\nFocus is on fish of type ",RTYPE,"\n")
cat("\nPrimary composition variable is ",Primary," and the secondary variable is ",Secondary, "\n")
cat("\nParametric bootstrap: B = ",B,"\n")
cat("\nAlpha = ",alph,"\n")

## Bootstrap function
AllBoot <- function(Tr,Pr,Gr) {

# Set up storage for bootstrap results
if( length(Gr) == 1 ) { theta.b <- matrix(numeric((3+nPgrp)*(B+1)),ncol=3+nPgrp)
} else { theta.b <- matrix(numeric((3+nPgrp+n2+ngrps)*(B+1)),ncol=3+nPgrp+n2+ngrps) }
# bootstrap loop
  for (b in 1:(B+1)) {
    if ( b == 1 ) {
      tstar <- Tr
      pstar <- Pr
    }
    else {
      tstar <- matrix(numeric(3*nw),ncol=3)
      for ( h in 1:nw ) {
        Eq0 <- which(Tr[h,] == 0)
        if( length(Eq0) == 3 )  tstar[h,] == c(0,0,0)  # if all probabilities are 0 then this stratum has no fish
        if( length(Eq0) == 2 )  { tstar[h,] <- Tr[h,] # if there are 2 0's, then probabilities are (1,0,0) or (0,1,0) or (0,0,1)
        } else {  # if there are at least 2 nonzero entries then we can generate bootstrap fish
        tstar[h,] <- t(rmultinom(1,nT[h],Tr[h,])) # These are number of wild fish for week h
        tstar[h,] <- tstar[h,]/nT[h] # Convert multinomial count into a proportion
        }  # end of else
      }  # end of for
        pstar <- matrix(numeric(nw*nPgrp),ncol=nPgrp)
      for ( h in 1:nw ) {
        Eq0 <- which(Pr[h,] == 0)
        if( length(Eq0) == nPgrp ) pstar[h,] <- rep(0,nPgrp) # if all probabilities are 0 then this stratum has no fish of type RTYPE
        else if( length(Eq0) == nPgrp-1 ) { pstar[h,] <- Pr[h,] # There is only one group with fish. No randomness here.
        } else pstar[h,] <- rmultinom(1,nPrime[h],Pr[h,])/nPrime[h]   # These are the bootstrap proportions by Wgrp for this stratum
      }  # end of for
    } # end of else from b==1 if statement
# Calculate RTYPE by week for the bootstrap data
    RearByWe <- c(WinData[,2])*tstar  # This is numbers by week for all rearing types
    RBW <- RearByWe  # This is redundant for the moment -- no second level collapsing
    NhatRear <- apply(RearByWe,2,sum)
    PrimeCats <- RBW[,Tcolumn] %*% pstar  # Bootstrap numbers of wild fish by group
# Now add PRIMARY by SECONDARY if needed
    if( length(Gr) != 1 ) {
      GrpAbund <- array(numeric(nPgrp*n2*nw),dim=c(nPgrp,n2,nw))
      if( b==1 ) {
        gstar <- Gr
      } else {
        gstar <- array(numeric(nPgrp*n2*nw),dim=c(nPgrp,n2,nw))
        for( h in 1:nw ) {
          ThisPrime <- RBW[h,Tcolumn]*pstar[h,]
          NE0 <- which(ThisPrime != 0)
          if( length(NE0) != 0 ) {
            for( jj in NE0 ) {
              Eq0 <- which(Gr[jj,,h] == 0)
              if( length(Eq0) == n2 ) gstar[jj,,h] <- rep(0,n2) # if all probabilities are 0 then this stratum has no fish from this grouping
              else if( length(Eq0) == n2-1 ) { gstar[jj,,h] <- Gr[jj,,h] # There is only one group with fish. No randomness here.
              } else { gstar[jj, ,h] <- rmultinom(1,ThisPrime[jj],Gr[jj,,h])/ThisPrime[jj] } } } } }
# gstar is set, now find group abundances
      for( h in 1:nw ) {
        ThisPrime <- RBW[h,Tcolumn]*pstar[h,]
        GrpAbund[,,h] <- as.vector(t( c(ThisPrime * gstar[,,h] )))

      } # end of h in 1:nw
      GrpTotals <- apply(GrpAbund,c(1,2),sum)
# rks added sum over Secondary variable so that two-way sums add to the correct number in each direction 2/16/2015
      theta.b[b,] <- c(NhatRear,PrimeCats,apply(GrpTotals,2,sum),as.vector(t(GrpTotals)))
      thetaNames <- c(levels(Types),levels(PrimaryNames),SecondaryNames,grpnams)
      colnames(theta.b) <- thetaNames
     } else {
      theta.b[b,]  <- c(NhatRear,PrimeCats)
      thetaNames <- c(levels(Types),PrimaryNames)
      colnames(theta.b) <- thetaNames }  # end of if length Gr
  if(b == (B+1) & writeThetas == TRUE) {write.csv(theta.b,file = paste(Run,"Thetas.csv",sep=""))}
  } # End of bootstrap loop

  Estimates <- theta.b[1,]
  theta.tots <- theta.b[,1:3]
  theta.Prime <- theta.b[,4:(3+nPgrp)]

# Find one-at-a-time confidence intervals for each statistic
  CI <- matrix(numeric((ncol(theta.b)*2)),ncol=2)
  for  (j in 1:ncol(theta.b)) {
    CI[j,] <- quantile(theta.b[,j],c(alph/2,1-alph/2))
  }

# Find simultaneous rectangular confidence intervals via Mandel/Betensky 2008
  MBCItot <- SCSrank(theta.tots,conf.level=1-alph)
  MBCItot <- MBCItot$conf.int
  MBCIprime <- SCSrank(theta.Prime,conf.level=1-alph)
  MBCIprime <- MBCIprime$conf.int
  if( length(Gr) > 1 ) {
   theta.Secondary <- theta.b[,(4+nPgrp):(3+nPgrp+n2)]
   MBCIsecondary <- SCSrank(theta.Secondary,conf.level=1-alph)
   MBCIsecondary <- MBCIsecondary$conf.int
   theta.grps <- theta.b[,(4+nPgrp+n2):ncol(theta.b)]
   MBCIgrps <- SCSrank(theta.grps,conf.level=1-alph)
   MBCIgrps <- MBCIgrps$conf.int
  } else {
    MBCIsecondary <- c(0.0)
    MBCIgrps <- c(0,0) }

  return( list(Estimates,CI,MBCItot,MBCIprime,MBCIsecondary,MBCIgrps) )
}
#################### End of bootstrap CI function ##############################
#                                                                              #
####################      MAIN                    ##############################
#                                                                              #

# Note.  Collapsed statistical weeks should be numbered 1,2,3,...
cat("\nStrata are collapsed according to: \n")
temp <- rbind(Windata$Strata,collaps)
rownames(temp) <- c("Week","Strata")
print(temp)
WinCounts <- mApply(Windata[,2],collaps,sum)
FinalStrata <- unique(collaps)  # These are the statistical weeks (i.e. strata) used for all subsequent analysis
WinData <- data.frame(FinalStrata,WinCounts) # This is the window data reduced to strata defined by collaps
nw <- nrow(WinData)   # This is the (perhaps) reduced number of statistical weeks (strata)

# Set up Trap data frame containing strata and rearing type
Fishdata <- droplevels(Fishdata[which(Fishdata$Rear!="NA"),])  # Drop those fish with no rearing (H,HNC,W) entry

# Recode WeekNumber in fish data to collapsed strata
jj  <- 0
chk <- 0
nullstrat <- NA
for( strat in Windata$Strata) {  # Note that the original Windata is used here
  jj <- jj+1
  juststrat <- Fishdata[Fishdata$WeekNumber ==  strat,]
  if( nrow(juststrat) == 0 ) {
            if( is.na(nullstrat[1]) ) nullstrat <- strat
            else nullstrat <- c(nullstrat,strat)
    } else { juststrat$Strata <- Windata[which(Windata$Strata == strat),3]
      if( chk == 0 ) {FishData <- juststrat
                      chk <- 1
      } else FishData <- rbind(FishData,juststrat) }
} # end of for strat
# Check to see if any of the original weeks are missing trapped fish
if( is.na(nullstrat[1]) ) { cat("\nThere are trapped fish for every week\n")
} else { cat("\nWeeks ",nullstrat," have no trapped fish.  Make sure that the window counts")
       cat(" for these weeks are collapsed with weeks having trapped fish.\n") }

rm( Windata )  # original Windata no longer needed
rm( Fishdata ) # Get rid of original Fishdata; new FishData has recoded statistical weeks
Trap <- data.frame(FishData$Strata,FishData$Rear) # Keep just the two columns of interest
names(Trap) <- c("Strata","Rear")
Tf <- table(Trap$Strata,Trap$Rear) # Frequency of each rearing type each week
rearN <- sum(Tf)
nT <- apply(Tf,1,sum) # This is the number of fish trapped each week
Trp <- Tf/nT   # Proportion of each rearing type each week

# Estimate the numbers for each rearing type (H,HNC,W) for each stratum
Rear.week <- c(WinData[,2])*Trp    # Note.  The c() function is a trick to get R to get the dimensions correct.
Rear.week.rounded <- round(Rear.week,0)
Totals <- round(apply(Rear.week,2,sum))
GrandTotal <- sum(Totals)
OverallEstimates <- c(Totals, GrandTotal)
names(OverallEstimates) <- c(names(Totals),"Total")
cat("\nEstimates of totals by rearing type and total count: \n")
print(OverallEstimates)
# and write rearing tables to csv file
RearingTables <- rbind(Tf,Trp,Rear.week.rounded,Totals)
if(writeOutput == TRUE) { write.csv(RearingTables,file = paste(Run,"Rearing.csv",sep="")) }

# Now determine which rearing type is of interest
Types <- unique(Trap$Rear)
Tcolumn <- which(levels(Types) == RTYPE)

# Grab fish of rearing type RTYPE
TheseFish <- FishData[FishData$Rear == RTYPE,]  # Note. If other columns are used, this must be changed.
rm(FishData) # We have everything we need out of the fish data

# Run the Primary analysis
#####
  if( Primary == "LGDFLmm") {  # Define large and small if Primary is fork length
    TheseFish$Size <- rep("Lg",nrow(TheseFish))
    TheseFish$Size[TheseFish$LGDFLmm < SizeCut] <- "Sm"
    TheseFish$LGDFLmm <- TheseFish$Size
  }
#####
  Pndx <- which(Primary == names(TheseFish))
  TheseFish <- droplevels(TheseFish[which(TheseFish[Pndx]!="NA"),])
  PN <- unique(TheseFish[Pndx])
  PrimaryNames <- PN[order(PN),]
  nPgrp <- length(PrimaryNames)
  FPrime <- table(factor(TheseFish$Strata, levels=FinalStrata),factor(TheseFish[,Pndx],levels=PrimaryNames))
  nPrime<- apply(FPrime,1,sum)
  primaryN <- sum(nPrime)
  PropPrime <- FPrime/nPrime
  PropPrime[is.na(PropPrime)] <- 0
  if ( pbtExpand == TRUE ) {
  FPrimeExp <- FPrime
    for (p in colnames(FPrime)[!colnames(FPrime) %in% "Unassigned"]) {
      if (p %in% pbtRate[,1]) {
        pbtndx <- which(pbtRate == p)
        FPrimeExp[,p] <- FPrime[,p] / pbtRate[pbtndx,2]
      } else { print(paste(p,"has no tag rates defined in tagging rate table")) }
      nPrimeExp <- apply(FPrimeExp,1,sum)
      totExp    <- nPrimeExp - nPrime
      FPrimeExp[,"Unassigned"] <- FPrimeExp[,"Unassigned"] - totExp
      }
    for (z in 1:nrow(FPrimeExp)) {
      if (FPrimeExp[z,"Unassigned"] < 0) {
        FPrimeExp[z,"Unassigned"] <- 0
        }
    }
  nPrimeExp <- apply(FPrimeExp,1,sum)
  PropPrimeExp <- FPrimeExp/nPrimeExp
  }
  if ( pbtExpand == TRUE ) {
    PrimaryAbundances <- c(Rear.week[,Tcolumn]) * PropPrimeExp
  } else {
    PrimaryAbundances <- c(Rear.week[,Tcolumn]) * PropPrime
  }
  PrimeEstimates <- round(apply(PrimaryAbundances,2,sum))
#  cat("\nEstimates for top level grouping are: \n")
#  print(PrimeEstimates)
  if ( pbtExpand == TRUE ) {
    forPrime <- rbind(FPrime,FPrimeExp,PropPrimeExp,round(PrimaryAbundances,0),PrimeEstimates)
  } else {
    forPrime <- rbind(FPrime,PropPrime,round(PrimaryAbundances,0),PrimeEstimates)
  }
  cat("\nFrequencies, proportions, and estimates by statistical week have been written to a .csv file for the top level grouping\n")
  if(writeOutput == TRUE) { write.csv(forPrime,file = paste(Run,"Prime.csv",sep="")) }

# Now add a secondary grouping variable if specified above
    if( !is.na(Secondary) ) {
      if( Secondary == "LGDFLmm") {  # Define large and small if Secondary is fork length
        TheseFish$Size <- rep("Lg",nrow(TheseFish))
        TheseFish$Size[TheseFish$LGDFLmm < SizeCut] <- "Sm"
        TheseFish$LGDFLmm <- TheseFish$Size
      }
      ndx <- which(Secondary == names(TheseFish))
      TheseFish <- droplevels(TheseFish[which(TheseFish[ndx]!="NA"),])
      secondaryN <- nrow(TheseFish)
      TheseFish$Combined <- paste(TheseFish[,Primary],TheseFish[,Secondary], sep = "")
      SecFreq <- table(TheseFish$Strata,TheseFish$Combined)
      SN <- unique(TheseFish[ndx])
      SecondaryNames <- SN[order(SN),1]
      n2 <- length(SecondaryNames)
      Props <- array(numeric(nPgrp*n2*nw),dim=c(nPgrp,n2,nw))
# Now concatenate the primary and secondary names
      grpnams <- NULL
      for( prim in PrimaryNames ) {
        for( nam in SecondaryNames ) grpnams = c(grpnams,paste0(nam,prim  ))
      }
      #grpnams <- grpnams[-1]
      ngrps <- length(grpnams)
# Figure out which primary groups for a week had positive estimates
       NE0.List <- lapply(1:length(FinalStrata), function(x) {
       NE0 <- which(PrimaryAbundances[x,] != 0) # NE0 is a vector of subscripts where the stocks have a positive estimate
         list(NE0 = NE0)  # return this into a list
         }
       )  # end of the lapply function
# Set up storage to accumulate primary group by secondary group tables over strata (statistical weeks)
       sumtable <- matrix(numeric(nPgrp*n2),ncol=n2)
# For each stratum, find the table of proportions of Primary by Secondary groups
       PGf <- array(numeric(nPgrp*n2*nw),dim=c(nPgrp,n2,nw))
       for( jj in FinalStrata ) {
         ThisStrata <- TheseFish[TheseFish$Strata == jj,]
         NDXtable <- table(factor(ThisStrata[,Pndx],levels=PrimaryNames),factor(ThisStrata[,ndx],levels=SecondaryNames))
         PGf[,,jj] <- NDXtable
         sumtable <- sumtable + NDXtable
         NDXrowtotals <- apply(NDXtable,1,sum)
         NDXrowtotals[NDXrowtotals == 0] <- 1
         Props[,,jj] <-  NDXtable/NDXrowtotals
       }
# Now go back and fix Props where there are fish for a stock-week, but none of those fish had sex (or age or size)
# so all proportions are zero. In that case, put in the average composition for the season.
         avgcomp <- apply(sumtable,2,sum)
         avgcomp <- avgcomp/sum(avgcomp)
         for( jj in FinalStrata ) {
           NE0ones <- unlist(NE0.List[[jj]])
           for( kk in NE0ones ) {
             if(max(Props[kk,,jj]) == 0) Props[kk,,jj] <- avgcomp } # This is the interesting line...there are projected
# fish that week for that stock, but that row of the stock by grp table is all zero.  If left zero, those fish would not be
# distributed and the totals over grp would not add to the stock estimate.
             if( jj == 1) {
               GrpAbundances <- as.vector(t( c(PrimaryAbundances[1,]) * Props[,,1] ))
             } else {
               GrpAbundances <- rbind(GrpAbundances,as.vector(t(c(PrimaryAbundances[jj,]) * Props[,,jj])))
             }
         }
         colnames(GrpAbundances) <- grpnams
         rownames(GrpAbundances) <- FinalStrata
         GrpAbundancesRnd  <- round(GrpAbundances)
         SecondaryTotals   <- round(apply(GrpAbundances,2,sum))
         SecondaryEsts     <- rbind(GrpAbundancesRnd,SecondaryTotals)
         for( jj in FinalStrata ) {
           if( jj == 1) {
             forSubgroups <- cbind(rep(jj,nPgrp),PGf[,,jj],Props[,,jj])
           } else{
             forSubgroups <- rbind(forSubgroups,cbind(rep(jj,nPgrp),PGf[,,jj],Props[,,jj]))
           }
         }
         colnames(forSubgroups) <- c("Strata",as.character(SecondaryNames),as.character(SecondaryNames))
         rownames(forSubgroups) <- rep(PrimaryNames,nw)
         if(writeOutput == TRUE) { write.csv(forSubgroups,file = paste(Run,"Subgroups.csv",sep="")) }
         if(writeOutput == TRUE) { write.csv(SecondaryEsts,file = paste(Run,"GrpEsts.csv",sep="")) }
         cat("\nSubgroup frequencies, proportions, and estimates by statistical week have been written to a .csv file\n")

         GrpTotals <- apply(GrpAbundances,2,sum)
         names(GrpTotals) <- grpnams
         SecondaryMatrix <- matrix(GrpTotals,ncol=nPgrp)
         colnames(SecondaryMatrix) <- as.character(PrimaryNames)
         rownames(SecondaryMatrix) <- as.character(SecondaryNames)
         SecondEsts <- apply(SecondaryMatrix,1,sum)
#         cat( "\nEstimates for ",Secondary," are: \n")
#         print(round(SecondEsts,1))
#         cat( "\nEstimates for ",Primary," and  ",Secondary," are: \n")
#         print(round(GrpTotals,1))
    }

# Now call the bootstrap routine which returns confidence intervals for everything
if( B > 0) {
  if(  is.na(Secondary) ) {
    if( pbtExpand == TRUE ) { AllCIs <- AllBoot(Trp,PropPrimeExp,0)
    } else { AllCIs <- AllBoot(Trp,PropPrime,0) }
  }
  if( !is.na(Secondary) ) {
    if( pbtExpand == TRUE ) { AllCIs <- AllBoot(Trp,PropPrimeExp,Props)
    } else { AllCIs <- AllBoot(Trp,PropPrime,Props) }
  }
    cat("\nEstimates and confidence intervals for numbers by rearing type\n")
    Rear <- round(AllCIs[[2]][1:3,])
    P1Rear <- round(100*(Rear[,2] - Rear[,1])/Totals/2,1)
    RearMB <- round(AllCIs[[3]])
    PsimRear <- round(100*(RearMB[,2] - RearMB[,1])/Totals/2,1)
    Rearing <- cbind(Totals,Rear,P1Rear,RearMB,PsimRear)
    colnames(Rearing) <- c("Estimates","L", "U","P1","Lsim","Usim","Psim")
    rownames(Rearing) <- as.character(levels(Types))
    print(Rearing)
    cat("\nP1 is the one-at-a-time percent half confidence width; Psim is the same for Mandel-Betensky intervals\n")
    cat("\nResults for ",Primary," for rearing type ",RTYPE,"\n")
    PrimeCIs <- round(AllCIs[[2]][4:(3+nPgrp),])
    P1prime <- round(100*(PrimeCIs[,2] - PrimeCIs[,1])/PrimeEstimates/2,1)
    PrimeMB <- round(AllCIs[[4]])
    PsimPrime <- round(100*(PrimeMB[,2] - PrimeMB[,1])/PrimeEstimates/2,1)
    PrimeResults <- cbind(PrimeEstimates,PrimeCIs,P1prime,PrimeMB,PsimPrime)
    colnames(PrimeResults) <-colnames(Rearing)
    rownames(PrimeResults) <- as.character(PrimaryNames)
    print(PrimeResults)
    if( !is.na(Secondary) ) {
      SecondaryTotal <- SecondEsts
      SecondaryTotal[SecondaryTotal == 0] <- 1
      SecondCIs <- round(AllCIs[[2]][(4+nPgrp):(3+nPgrp+n2),])
      P1second <- round(100*(SecondCIs[,2] - SecondCIs[,1])/SecondaryTotal/2,1)
      SecondMB <- round(AllCIs[[5]])
      PsimSecond <- round(100*(SecondMB[,2] - SecondMB[,1])/SecondaryTotal/2,1)
      SecondResults <- cbind(round(SecondaryTotal),SecondCIs,P1second,SecondMB,PsimSecond)
      colnames(SecondResults) <-colnames(Rearing)
      rownames(SecondResults) <- SecondaryNames
      cat("\nResults for ",Secondary," for rearing type ",RTYPE,"\n")
      print(SecondResults)
      GrpTot <- GrpTotals
      GrpTot[GrpTot == 0] <- 1
      #GrpCIs <- round(AllCIs[[2]][(4+nPgrp):(3+nPgrp+ngrps),])
      GrpCIs <- round(AllCIs[[2]][(4+nPgrp+n2):nrow(AllCIs[[2]]),])
      P1grp <- round(100*(GrpCIs[,2] - GrpCIs[,1])/GrpTot/2,1)
      GrpMB <- round(AllCIs[[6]])
      PsimGrp <- round(100*(GrpMB[,2] - GrpMB[,1])/GrpTot/2,1)
      GrpResults <- cbind(round(GrpTotals),GrpCIs,P1grp,GrpMB,PsimGrp)
      colnames(GrpResults) <-colnames(Rearing)
      rownames(GrpResults) <- grpnams
      cat("\nResults for ",Secondary," by ",Primary," for rearing type ",RTYPE,"\n")
      print(GrpResults)
    }
   }
   if(is.na(Secondary)) {CIs <- rbind(Rearing,PrimeResults)} else {CIs <- rbind(Rearing,PrimeResults,SecondResults,GrpResults)}
   if(writeOutput == TRUE){
   write.csv(CIs,file = paste(Run,"CIs.csv",sep=""))
   write.table(paste("RearSampleSize =",rearN),file = paste(Run,"CIs.csv",sep=""),append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
   write.table(paste("PrimarySampleSize =",primaryN),file = paste(Run,"CIs.csv",sep=""),append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
   if(!is.na(Secondary)) { write.table(paste("SecondarySampleSize =",secondaryN),file = paste(Run,"CIs.csv",sep=""),append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE) }
   }
   cat("\nEnd time: ",date(),"\n")

   #########################################################################################
   #                            SCOBI.r                                 Febuary 2015       #
   #  Copyright                                                         Kirk Steinhorst    #
   #                                                                                       #
   #  "SCOBI.r" is used to find estimates and confidence intervals for two composition     #
   #  factors--sex or age or genetic stock or size (for steelhead.  One is listed as the   #
   #  PRIMARY factor and the other is listed as the SECONDARY factor.  Estimates are given #
   #  for the PRIMARY factor and the estimates for the SECONDARY factor add to the PRIMARY #
   #  factor estimates.                                                                    #
   #                                                                                       #
   #  Estimates are given for H, HNC, and W escapement.  Composition estimates are         #
   #  given for fish of a selected type -- H or HNC or W.  Additionally, the PRIMARY       #
   #  factor can be broken down one more level if requested.                               #
   #                                                                                       #
   # Permissible column names in the excel window counts tab are Strata, a count variable, #
   # and Collaps.  On the fish data tab, you can have any number of columns with any       #
   # number of names, but the names must include WeekNumber and the names of the PRIMARY   #
   # and SECONDARY composition factors from GenSex, BY, GenStock and, for steelhead,       #
   # LGDFLmm.                                                                              #
   #                                                                                       #
   #  Window counts are collapsed into "statistical weeks" = strata defined by the user.   #
   #  The collapsing pattern must collapse window counts in weeks where no trapped fish    #
   #  are available into adjacent weeks where trapped fish are available.  Clearly we      #
   #  would like to have a good composition sample for each statistical week, but          #
   #  practically speaking we don't need to know composition well for weeks when few fish  #
   #  of the type of interest (H, HNC, or W) are passing.                                  #
   #                                                                                       #
   #  Outputs include estimates and confidence intervals for numbers of H, HNC, and W fish #
   #  passing Lower Granite dam.  You also get estimates and confidence intervals for      #
   #  the PRIMARY (and SECONDARY, if requested) composition factor(s).                     #
   #                                                                                       #
   #  Three to five .csv files are produced -- *Rearing.csv, *Prime.csv, *CIs.csv, and     #
   #  (if subgroups are requested) *GrpEsts.csv and *Subgroups.csv. The rearing file       #
   #  contains frequencies, proportions, and estimates by statistical week for H, HNC, and #
   #  W. The prime file contains frequencies, proportions, and estimates by statistical    #
   #  week for the PRIMARY composition factor. The CIs file contains abundance estimates   #
   #  and confidence intervals (both one-at-a-time and simultaneous) for the PRIMARY (and  #
   #  SECONDARY if requested) composition factors. The GrpEsts file contains abundance     #
   #  estimates for the SECONDARY composition groups by statistical week. The Subgroups    #
   #  file contains the frequency and proportions of the SECONDARY groups in a             #
   #  Primary x Secondary format.                                                          #
   #                                                                                       #
   #  Be sure to change the names of these files or change directories before making       #
   #  another run or these .csv files will be overwritten.                                 #
   #                                                                                       #
   #  Copyright February 16, 2015  Kirk Steinhorst                                         #
   #                                                                                       #
   #########################################################################################
}

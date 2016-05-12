#' @title SCRAPI v2.1 The smolt companion to SCOBI
#'
#' @description Perform compositional analyses of smolts at Lower Granite Dam.
#'
#' @param smoltData the .csv file containing the biological data for the smolts to be analyzed. The file should contain data for all smolts
#' trapped in a given migratory year and for a given 'species' (sthd, ch0, or ch1). The function \code{lgr2SCRAPI()} can be used to format
#' raw smolt data exported from the LGTrappingDB to make it ready for the \code{SCRAPI()} function
#' @param Dat speciy the column in \code{smoltData} that contains the sample date for each smolt
#' @param Rr specify the column in \code{smoltData} that contains the rear type (W or HNC) for each smolt
#' @param Primary specify the primary category in \code{smoltData} to be estimated. i.e., after \code{RTYPE} what category would you like to decompose?
#' @param Secondary the secondary category in the \code{smoltData} to be estimated. i.e. after \code{Primary} what category would you like to
#' decompose next? Use \code{Secondary = NA} if no secondary decomposition is desired
#' @param passageData  the .csv file containg the smolt passage data. The file should contain the sampling date for each day of the season and
#' for each date the 1) calendar week, 2) sampling rate, 3) count of smolts in the trap, and 4) the juvenile bypass system guidance efficiency. \code{passageData}
#' should also contain the collapsing scheme to be used.
#' @param strat specify the column in \code{passageData} that contains the calendar weeks
#' @param dat specify the column in \code{passageData} that contains the sampling dates
#' @param tally specify the column in \code{passageData} containg the daily count of smolts in the trap
#' @param samrate specify the column in \code{passageData} containing the daily trap sample rate
#' @param guidance specify the column in \code{passageData} containing the daily guidance efficiency estimates
#' @param collaps specify the column in \code{passageData} containing the collapsing scheme to be used
#' @param Run synopsis of the run being conducted. This will be used as the prefix for all of your output files
#' @param RTYPE which rear type (wild [H] or hatchery unclipped [HNC]) would you like to perform analysis on
#' @param REARSTRAT would you like pWild calculated on a time-stratified basis (TRUE) or across the entire emigration (FALSE)
#' @param alph the alpha used for confidence intervals (e.g., \code{alph} = 0.10 results in 90 percent CIs)
#' @param B the desired number of bootstraps
#' @param dateFormat What format are the dates in your fish and passage data? The default format used in the LGTrappingDB
#' is 01/01/1900
#'
#' @author Kirk Steinhorst and Mike Ackerman
#'
#' @import stats plyr
#' @export
#' @return NULL

#            SCRAPI 2.1
# This R program calculates bootstrap confidence intervals for total wild (or hatchery)
# and wild (or hatchery) by group (of your choice) for smolts arriving at Lower Granite dam.
# Groups can be made up of a single catetgorical variable (sex, age, genstock, or size) or
# two categorical variables--one nested within the other.
#
# Bootstrap confidence intervals account for variability in counts at Lower Granite Dam
# and variability in sampling smolts at the dam.
#
# The excel input file has 2 tabs--fish data and passage data
# The column headings are arbitrary--you can use any names you want (without spaces though).
# The column names will be supplied in the "Things that might change" section.
# The fish tab should include columns for stratum, date, rearing condition (W,H), a PRIMARY categorization variable
# such as sex, age, genstock... and (perhaps) a SECONDARY categorization variable.
# The passage tab should include stratum, date, smolt count, trap sampling rate, guidance efficiency, and collapsing pattern
#
# Version 2.0 adds variability for proportion wild and hatchery to the bootstrap process
# Previous versions calculated proportion wild versus hatchery just once and left it fixed.
# Version 2.1 weights rearing status (HNC and W) by true sampling rate
#
# Copyright Kirk Steinhorst May 9, 2016

SCRAPIv2.1 <- function(smoltData = NULL, Dat = "CollectionDate", Rr = "Rear", Primary = "GenStock", Secondary = NA, passageData = NULL,
                       strat = "Week", dat = "SampleEndDate", tally = "SampleCount", samrate = "SampleRate", guidance = "GuidanceEfficiency",
                       collaps = "Collapse", Run = "output", RTYPE = "W", REARSTRAT = TRUE, alph = 0.1, B = 5000, dateFormat = "%m/%d/%Y")
{
  # Import data
  All  <- read.csv(file = smoltData, header = TRUE)
  pass <- read.csv(file = passageData, header = TRUE)

  # Write header
  cat("\nStart time: ",date(),"\n")
  cat("\nThis is a run of ", Run, "\n")
  cat("\nFocus is on fish of type ",RTYPE,"\n")
  cat("\nPrimary composition variable is ",Primary," and the secondary variable is ",Secondary, "\n")
  cat("\nWild adjustment by week is ", REARSTRAT,"\n")
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
  # passage <- data.frame(Stratum=pass[,PASScollaps],Tally=pass[,PASScounts],Ptrue=pass$true)
  # RearDat should have Stratum, Rear, and True sampling rate
  # Fish should have stratum, PRIMARY category, SECONDARY category (if needed), and
  # realized sampling rate (true sampling rate x sampling proportion for PRIMARY category) for each fish.

  thetahat <- function(passage,RearDat,Fish){
    dailypass <- passage$Tally/passage$Ptrue
    bystrata <- mApply(dailypass,passage$Stratum,sum) # Estimates of total smolts by strata
    # Calculate proportion HNC and W by strata
    HNCWstrat <- mApply(1/RearDat$True,list(RearDat$Stratum,RearDat$Rear),sum)
    HNCWstrat[is.na(HNCWstrat)] <- 0
    if( ncol(as.data.frame(HNCWstrat)) == 1 ) PWild <- 1 else {
      HNCWprop <- prop.table(HNCWstrat,margin=2)
      if(REARSTRAT == FALSE) {
        ColTotals   <- apply(HNCWstrat,1,sum)
        Proportions <- ColTotals/sum(ColTotals) # This has pooled values for proportion HNC and W
        PWild <- Proportions[2]
      } else { PWild  <- HNCWprop[2,] }
    }
    WildStrata <- PWild*bystrata  # This is the estimate of wild smolts by strata (or pooled)
    # bystrata <- pWild*bystrata # This adjusts the (statistical) weekly passage to just wild smolts
    # I changed to WildStrata from bystrata so as not to use bystrata in two ways
    total <- sum(WildStrata)
    Primarystrata <- mApply(1/Fish$SR,list(Fish$Strat,Fish$PGrp),sum)
    Primarystrata[is.na(Primarystrata)] <- 0
    Primaryproportions <- prop.table(Primarystrata,margin=2)
    Primaryests <- Primaryproportions%*%WildStrata
    if( !is.na(Secondary) ){
      SecondAbund <- array(numeric(nPgrps*nSgrps*nstrats),dim=c(nPgrps,nSgrps,nstrats))
      # Delete fish with NA for SECONDARY categorical variable
      Fish <- droplevels(Fish[which(Fish$SGrp!="NA"),])
      nsFish <<- nrow(Fish)
      Freqs <- mApply(1/Fish$SR,list(factor(Fish$PGrp,levels=Pgrps),factor(Fish$SGrp,levels=Sgrps),factor(Fish$Strat,levels=strats)),sum)
      Freqs[is.na(Freqs)] <- 0
      Props <- prop.table(Freqs,margin=c(1,3))
      # Check to see if the proportions for the Secondary category are NaNs because the corresponding Freqs are all 0.
      # If so, put in the average proportions.
      if( any(is.nan(Props)) ) {
        AvgProp <- getAvgProp(Fish)
        for( h in 1:nstrats ) {
          for( i in 1:nPgrps ) {
            if( any(is.nan(Props[i,,h])) ) Props[i,,h] <- AvgProp[,i]
          }
        }
      }
      for( h in 1:nstrats ) {
        ThisPrime <- WildStrata[h]*Primaryproportions[,h]
        SecondAbund[,,h] <- as.vector(t( c(ThisPrime * Props[,,h] )))
      }
      PrimeBySecond <- apply(SecondAbund,c(1,2),sum)
      PrimeBySecond[is.na(PrimeBySecond)] <- 0
      assign("sTable", PrimeBySecond, envir = .GlobalEnv)
      Second <- apply(PrimeBySecond,2,sum)
      return( list(total,WildStrata,Primaryproportions,Primaryests,Second,PrimeBySecond) )
    } else   return( list(total,WildStrata,Primaryproportions,Primaryests) )

  } # end of thetahat function

  bootsmolt <- function(FishWH,FishDat,LGDdaily){

    # Set up storage for bootstrap results and plug in estimates for line 1
    theta.b <- matrix(numeric(p*B),ncol=p)
    # bootstrap loop
    for (b in 1:B) {
      if( b == 1 ){
        dailyStar <- passdata
        RearStar <- RearData
        indivStar <- FishDat
      } else {
        dailypass <- round(LGDdaily$Tally/LGDdaily$Ptrue)
        # Find bootstrap values for number of smolts passing LGR each day
        cntstar = numeric(ndays)
        for (i in 1:ndays) {
          if(dailypass[i] != 0 ) cntstar[i] <- rbinom(1,dailypass[i],LGDdaily$Ptrue[i])
        } # end of for loop on daily smolts
        dailyStar <- data.frame(Stratum=LGDdaily$Stratum,Tally=cntstar,Ptrue=LGDdaily$Ptrue)

        # Find a weighted bootstrap sample of FishWH BY stratum
        H <- 0
        for ( h in strats) {
          justwk <- FishWH[FishWH$Stratum==h,]    # grab All smolts of a given statistical week
          nwk <- nrow(justwk)
          i <- sample.int(nwk,replace=TRUE,prob=unlist(justwk$True))
          wkstar <- justwk[i,]
          if ( H == 0 ) {WHstar <- wkstar
          H <- 1
          } else { WHstar <- rbind(WHstar,wkstar)  }
        } # end of for h

        # Now find a bootstrap sample of FishData (i.e. PRIMARY fish) BY strata using a weighted bootstrap
        H <- 0
        for ( h in strats) {
          justwk <- FishDat[FishDat$Strat==h,]    # grab All smolts of a given statistical week
          nwk <- nrow(justwk)
          i <- sample.int(nwk,replace=TRUE,prob=unlist(justwk$SR))
          wkstar <- justwk[i,]
          if ( H == 0 ) {indivStar <- wkstar
          H <- 1
          } else { indivStar <- rbind(indivStar,wkstar)  }
        }
      }

      #  Now we have (real or) bootstrap data for passage, rearing data, and PRIMARY fish data--find estimates
      eststar <- thetahat(dailyStar,RearStar,indivStar)
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

  # Add true sampling rate for dates to each fish
  # Find dates where smolt counts were made AND fish samples were taken.
  set    = intersect(All[,FISHdate],pass[,PASSdate])
  ndates <- length(set)
  All$true = numeric(nAll)
  for ( nn in 1:ndates ) {
    passtemp <- pass[as.Date(pass[,PASSdate], format = dateFormat) == as.Date(set[nn], origin="1970-01-01", format = dateFormat),]
    All$true[as.Date(All[,FISHdate], format = dateFormat) == as.Date(set[nn], origin="1970-01-01" , format = dateFormat)] <- passtemp$true
  }

  # Set up data frame for proportion wild
  # Note that we assume each smolt has a rear type of H or W.  If any are NA, then we need to add a "droplevels".
  RearData <- data.frame(Rear=All[,FISHrear],Stratum=All$Collaps,True=All$true)

  # Calculate the proportion of smolts that are wild and hatchery.  Note.  This section could be deleted.  I kept it to conform to the
  #                                                                        previous printouts.  Down to 'print( round(c(WildCollaps,TotalWild)) )'.
  if(nlevels(RearData$Rear) == 1) ProWild <- 1 else {
    HNCWstrata <- mApply(1/RearData$True,list(RearData$Stratum,RearData$Rear),sum)
    HNCWstrata[is.na(HNCWstrata)] <- 0
    HNCWproportions <- prop.table(HNCWstrata,margin=2)
    HNCWests <- HNCWproportions%*%passcollaps
    #Tabl  <- table(RearData$Rear,RearData$Stratum)
    #PropWild    <- prop.table(Tabl,margin=2)
    PercentHNCW <- 100*round(HNCWproportions,3)
    ColTotals   <- apply(HNCWstrata,1,sum)
    cat("\nTrapped hatchery and wild by stratum\n")
    print( cbind(HNCWstrata,ColTotals) )
    Proportions <- ColTotals/sum(ColTotals) # This has pooled values for proportion H and W
    Ovrallpct <- 100*round(Proportions,3)
    cat("\nPercent hatchery and wild by stratum\n")
    print( cbind(PercentHNCW,Ovrallpct) )
    if(REARSTRAT == TRUE) ProWild  <- HNCWproportions[2,] else ProWild <- Proportions[2]
  }
  WildCollaps <- ProWild*passcollaps
  TotalWild   <- sum(WildCollaps)
  cat("\nWild smolts by stratum\n")
  print( round(c(WildCollaps,TotalWild)) )

  # Select those fish of RTYPE.  Note.  We assume that all trapped smolts have an HNC or W designation.
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
  strats  <- unique(Cpattern[,2])
  nstrats <- length(strats)

  # Find number of fish sampled for the PRIMARY categorical variable by date
  tabl   <- table(AllPrimary[,FISHdate],AllPrimary[,FISHpndx])
  nPrime <- apply(tabl,1,sum)

  # The realized sampling rate for a date is pass$true x proportion of fish sampled for the PRIMARY categorical variable on that date.
  # Put the realized sampling rate on each smolt record, i.e.,
  # ID smolts of a given date and attach the realized sampling rate to each.
  AllPrimary$SR = numeric(nFISH)
  for ( nn in 1:ndates ) {
    passtemp <- pass[as.Date(pass[,PASSdate], format = dateFormat) == as.Date(set[nn], origin="1970-01-01", format = dateFormat),]
    Primecount <- nPrime[which( as.Date(names(nPrime), format = dateFormat) ==  as.Date(set[nn], origin="1970-01-01", format = dateFormat))]
    AllPrimary$SR[as.Date(AllPrimary[,FISHdate], format = dateFormat) == as.Date(set[nn], origin="1970-01-01", format = dateFormat)] <- passtemp$true*Primecount/passtemp[,PASScounts]
  }

  # Keep week, primary (and if needed secondary) group(s), and realized sampling rate
  if( !is.na(Secondary) ){
    AllPrime <- data.frame(Strat=AllPrimary$Collaps,PGrp=AllPrimary[,FISHpndx],SGrp=AllPrimary[,FISHsndx],SR=AllPrimary$SR)
    Sgrps <- unique(as.character(AllPrime$SGrp))
    Sgrps <- sort(Sgrps[Sgrps!="NA"])
    nSgrps <- length(Sgrps)
    p <- p + nSgrps + nPgrps*nSgrps
    AvgProp <- getAvgProp(AllPrime)
  } else {
    AllPrime <- data.frame(Strat=AllPrimary$Collaps,PGrp=AllPrimary[,FISHpndx],SR=AllPrimary$SR)
  }
  # Run thetahat() function  Note. This is not necessary now that bootsmolt has a 'if b == 1' step
  ests <- thetahat(passdata,RearData,AllPrime)
  pPropTable <- t(ests[[3]])
  if(!is.na(Secondary)) sAbunTable <- ests[[6]]

  # Run bootsmolt() function
  answer <- bootsmolt(RearData,AllPrime,passdata)
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
  if(!is.na(Secondary)) { write.table(paste("SecondSampleSize =",nsFish), file = ciFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE) }

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


#' @title lgr2SCOBI: Format adult data from LGTrappingDB for SCOBI
#'
#' @description lgr2SCOBI() is a function to convert raw data exported or 'dumped' out of the LGTrappingDB to a format that is ready
#' to be analyzed by the \link{SCOBI}() function. It is to be used on one species/spawn year combination (e.g., SY2015 steelhead) at a time. It will
#' write a csv file ready to be used as the \code{adultData} argument of the SCOBI() function.
#'
#' @param input the name of the csv file containing adult data exported from the LGTrappingDB. Can also accept an \code{object} of a similar
#' format. See \code{\link[SCOBI]{exRawSthdAdultData}}
#' @param species what is the species of the data you are attempting to format for SCOBI. The default option is \code{"chnk"} for
#' Chinook salmon. Simply set as \code{"sthd"} for steelhead.
#' @param exportFile What would you like to name your exported file containing your data formatted for SCOBI analysis? By
#' default \code{lgr2SCOBI()} exports the formatted data as a csv file which is the preferred format for SCOBI.
#'
#' @author Mike Ackerman
#'
#' @examples lgr2SCOBI(input = exRawSthdAdultData, species = "sthd", exportFile = "SthdScobiInput")
#'
#' @import stringr car
#' @export
#' @return NULL

lgr2SCOBI <- function(input = NULL, species = "chnk", exportFile = NULL)
{
  # IMPORT UNFORMATTED DATA DUMPED FROM THE LGRTrappingDB
  if(is.character(input) == TRUE) { rawData <- read.table(file = input, header = TRUE, sep = ",", na.strings = c("","NA"), comment.char = "")
  } else { rawData <- input }
  data   <- subset(rawData, select = c("WeekNumber","CollectionDate","SpawnYear","MasterID","BioSamplesID","SRR","LGDMarkAD","LGDFLmm","GenSex",
                                       "GenStock","GenStockProb","BioScaleFinalAge","GenPBT_ByHat","GenPBT_RGroup","GenParentHatchery","GenBY","BiosamplesValid","LGDValid",
                                       "LGDNumPIT"))

  # CREATE REAR COLUMN TO POPULATE WITH W, H, or HNC
  data[,"Rear"] <- NA
  w <- substr(data$SRR,3,3)=="W" & data$LGDMarkAD=="AI"
  data[w,"Rear"] <- "W"
  h <- substr(data$SRR,3,3)=="H" & data$LGDMarkAD=="AD"
  data[h,"Rear"] <- "H"
  hnc <- substr(data$SRR,3,3)=="H" & data$LGDMarkAD=="AI"
  data[hnc,"Rear"] <- "HNC"

  # IN THE EXISTING GenSex COLUMN, REPLACE NG, U, and (blank) with NA
  data$GenSex[data$GenSex == ""]   <- NA
  data$GenSex[data$GenSex == "NG"] <- NA
  data$GenSex[data$GenSex == "U"]  <- NA
  data$GenSex <- factor(data$GenSex)

  # IN THE EXISTING GenStock COLUMN, REPLACE NG and (blank) with NA
  data$GenStock[data$GenStock == ""]   <- NA
  data$GenStock[data$GenStock == "NG"] <- NA
  data$GenStock <- factor(data$GenStock)

  # CREATE MPG COLUMN. For steelhead, replace UPSALM, MFSALM, SFSALM and LOSALM w/ SALMON & replace UPCLWR, SFCLWR, and LOCLWR w/ CLRWTR.
  # For Chinook, replace CHMBLN w/ MFSALM
  data[,"MPG"] <- data[,"GenStock"]
  if(species == "sthd") {
    data$MPG <- recode(data$MPG,"c('UPSALM','MFSALM','SFSALM','LOSALM')='SALMON'")
    data$MPG <- recode(data$MPG,"c('UPCLWR','SFCLWR','LOCLWR')='CLRWTR'")
  }
  if(species == "chnk") {
    data$MPG[data$MPG == "CHMBLN"] <- "MFSALM"
  }

  # CREATE fwAge COLUMN. THIS COLUMN WILL INCLUDE THE FRESHWATER AGE ASSIGNED TO EACH FISH. FRESHWATER AGE IS LEFT OF THE COLON IN THE
  # BioScaleFinalAge COLUMN.
  data[,"fwAge"] <- substr(data$BioScaleFinalAge,1,1)
  data$fwAge <- recode(data$fwAge,"c('N','?','')=NA")

  # CREATE swAge COLUMN. THIS COLUMN WILL INCLUDE THE SALTWATER AGE ASSIGNED TO EACH FISH. SALTWATER AGE IS EVERYTHING RIGHT OF THE COLON
  # IN THE BioScaleFinalAge COLUMN.
  data[,"swAge"] <- gsub(".*:","",data$BioScaleFinalAge)
  if(species == "chnk") { data$swAge[data$swAge == "MJ"] <- 0 }
  data$swAge <- recode(data$swAge,"c('A','?','')=NA")
  data$swAge[is.na(data$fwAge)] <- NA    # if fwAge == NA, change swAge to NA
  data$swAge <- gsub("S","R",data$swAge) # replace "S" with "R" in swAge column
  data$swAge <- gsub("s","R",data$swAge) # replace "S" with "R" in swAge column

  # CREATE NEW Age COLUMN THAT CONCATENATES THE fwAge & swAge
  data[,"Age"] <- NA
  for (i in 1:nrow(data)){
    if( is.na(data$fwAge[i]) ) data[i,"Age"] <- NA else data[i,"Age"] <- paste("F",data$fwAge[i],"S",data$swAge[i],sep="")
  }

  # CREATE NEW totalAge COLUMN THAT WILL BE USED TO DETERMINE THE BROOD YEAR (BY) FOR EACH FISH
  data[,"totalAge"] <- NA
  if(species == "chnk") {
    for (i in 1:nrow(data)){
      if( is.na(data$Age[i]) ) data[i,"totalAge"] <- NA else data[i,"totalAge"] <- sum(as.numeric(data$fwAge[i]),as.numeric(data$swAge[i]),1)
    }
  }
  if(species == "sthd") {
    data[,"totalAge"] <- NA
    data$totalAge <- gsub("R",1,data$swAge)
    for (i in 1:nrow(data)){
      if( is.na(data$Age[i]) ) data[i,"totalAge"] <- NA else data[i,"totalAge"] <- sum(as.numeric(unlist(strsplit(data$totalAge[i],""))),as.numeric(data$fwAge[i]),1)
    }
  }

  # We’re now done using the steelhead repeat spawner ages in the swAge column. Any repeat spawners (includes an ‘R’)
  # now need to have their swAge changed to ‘Repeat’
  if(species == "sthd") {
    for (i in 1:nrow(data)){
      if(grepl("^[^_]+R",data$swAge)[i] == TRUE) data$swAge[i] <- "Repeat"
    }
  }

  # CREATE NEW BY COLUMN. THIS WILL CONTAIN THE BY FOR EACH FISH.
  data[,"BY"] <- as.numeric(substr(data$SpawnYear,3,7)) - as.numeric(data$totalAge)
  #data$BY[!is.na(data$GenBY)] <- data$BY
  for(i in 1:nrow(data)){
    if( !is.na(data$GenBY[i]) ) data$BY[i] <- data$GenBY[i] else data$BY[i] <- data$BY[i]
    if( !is.na(data$GenBY[i]) ) data[i,c("fwAge","swAge","Age","totalAge")] <- NA
  }

  # For any fish where BioSamplesValid = 0, (blank), or NA we need to replace the cell contents with NA for the following columns:
  # fwAge, swAge, Age, totalAge, BY. This removes fish that were detected at PIT tag arrays that are not needed for aggregate age composition.
  for(i in 1:nrow(data)){
    if( is.na(data$BiosamplesValid[i]) ) data[i,c("fwAge","swAge","Age","totalAge","BY")] <- NA
  }

  # BEGINNING AND END OF STEELHEAD SEASON. Fish collected on or after 7/1 should be changed to 27A, the beginning of the season.
  # Fish collected before 7/1 should be changed to 27B, the end of the season.
  if(species == "sthd") {
    sy <- as.numeric(substr(data$SpawnYear[1],3,7))
    a <- data$WeekNumber == 27 & as.Date( as.character(data$CollectionDate), "%m/%d/%Y") >= paste(sy-1,"-07-01", sep = "")
    data[a,"WeekNumber"] <- "27A"
    b <- data$WeekNumber == 27 & as.Date( as.character(data$CollectionDate), "%m/%d/%Y") <= paste(sy,"-06-30", sep = "")
    data[b,"WeekNumber"] <- "27B"
  }

  write.table(data, file = paste(exportFile,".csv",sep=""), append = FALSE, quote = FALSE, sep = ",", row.names = FALSE)

}

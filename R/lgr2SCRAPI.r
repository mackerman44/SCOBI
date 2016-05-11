#' @title lgr2SCRAPI
#'
#' @description lgr2SCRAPI() is a function to convert raw data exported or 'dumped' out of the LGTrappingDB to a format that is ready
#' to be analyzed by the \link{SCRAPI}() function. It is to be used on one species/migratory year combination (e.g., MY2015 steelhead) at
#' a time. It will write a csv file ready to be used as the \code{smoltData} argument of the SCRAPI() function.
#'
#' @param input the name of the csv file containing smolt data exported from the LGTrappingDB. Can also accept an \code{object} of a similar
#' format. See \code{\link[SCOBI]{exRawSthdSmoltData}}
#' @param species what is the species of the data you are attempting to format for SCRAPI. The default option is \code{"chnk"} for
#' Chinook salmon. Simply set as \code{"sthd"} for steelhead.
#' @param exportFile What would you like to name your exported file containing your data formatted for \code{SCOBI} analysis? By
#' default \code{lgr2SCRAPI()} exports the formatted data as a csv file which is the preferred format for SCRAPI.
#'
#' @author Mike Ackerman
#'
#' @examples lgr2SCRAPI(input = exRawSthdSmoltData, species = "sthd", exportFile = "SthdScrapiInput")
#'
#' @import stringr car
#' @export
#' @return NULL

lgr2SCRAPI <- function(input = NULL, species = "chnk", exportFile = NULL)
{
  # IMPORT UNFORMATTED SMOLT DATA DUMPED FROM THE LGTrappingDB
  if(is.character(input) == TRUE) { rawData <- read.table(file = input, header = TRUE, sep = ",", na.strings = c("","NA"), comment.char = "")
  } else { rawData <- input }
  data   <- subset(rawData, select = c("MasterID","WeekNumber","CollectionDate","SRR","GenSex","GenStock","BioScaleFinalAge","BioSamplesID",
                                       "SpawnYear","LGDFLmm","GenRear","LGDLifeStage","GenStockProb","GenParentHatchery","GenBY","GenRun",
                                       "LGDMarkAD"))

  # CREATE REAR COLUMN TO POPULATE WITH W or HNC
  data[,"Rear"] <- NA
  w <- substr(data$SRR,3,3)=="W"
  data[w,"Rear"] <- "W"
  hnc <- substr(data$SRR,3,3)=="H"
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
  if(species == "sthd") {
    data[,"fwAge"] <- substr(data$BioScaleFinalAge,1,1)
    data$fwAge <- recode(data$fwAge,"c('N','?','')=NA")
  }

  write.table(data, file = paste(exportFile,".csv",sep=""), append = FALSE, quote = FALSE, sep = ",", row.names = FALSE)
}

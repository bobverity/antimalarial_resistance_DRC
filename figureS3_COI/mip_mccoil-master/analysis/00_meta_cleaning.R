# DHS/DRC

# grab our samples data
miss <- readRDS("analysis/data/raw/samples_MC.rds")
non_miss <- readRDS("analysis/data/raw/samples_nonMC.rds")
ID_split <- rbind(miss, non_miss)

# # grab the meta and first make it unique (there was some duplication with a bad merge)
old_meta <- data.table::fread("analysis/data/raw/meta_old.txt")
old_meta <- unique(old_meta)
#
# # subset for the missing children
missing <- old_meta[grepl("Missing",old_meta$`Stock Plate`),]

# sample names to barcode
ID_split$Barcode <- gsub("([0-9]{1,4})([A-Z][0-9][A-Z][0-9][A-Z])","\\2",ID_split$ID)
ID_split$DHS <- grepl("[A-Z][0-9][A-Z][0-9][A-Z]",ID_split$Barcode)

# barcode to lat/long/cluster/year
library(rdhs)
drcs <- dhs_datasets(countryIds = "CD",fileFormat = "FL", fileType = c("GE","PR"))
dats <- get_datasets(drcs)
dat_read <- lapply(dats,readRDS)
barcs <- search_variable_labels(names(dats),"bar|Cluster number|blood sample|Case Identification|^Line number$")
extr  <- extract_dhs(barcs,add_geo = TRUE)

# make the barcode for the individual as taken from different sources often
extr$CDPR50FL$Barcode <- apply(extr$CDPR50FL[,4:5],1,function(x){
  ret <- unique(grep("[A-Z][0-9][A-Z][0-9][A-Z]",x,value=TRUE))
  if(length(ret)==0){
    ret<-NA
  }
  return(ret)
})
extr$CDPR61FL$HIV_Barcode <- apply(extr$CDPR61FL[,5:6],1,function(x){
  ret <- unique(grep("[A-Z][0-9][A-Z][0-9][A-Z]",x,value=TRUE))
  if(length(ret)==0){
    ret<-NA
  }
  return(ret)
})
extr$CDPR61FL$VaccSmear_Barcode <- apply(extr$CDPR61FL[,c(4,7)],1,function(x){
  ret <- unique(grep("[A-Z][0-9][A-Z][0-9][A-Z]",x,value=TRUE))
  if(length(ret)==0){
    ret<-NA
  }
  return(ret)
})
extr$CDPR61FL$Barcode <- apply(extr$CDPR61FL[,c(4:7)],1,function(x){
  ret <- unique(grep("[A-Z][0-9][A-Z][0-9][A-Z]",x,value=TRUE))
  if(length(ret)==0){
    ret<-NA
  }
  return(ret)
})

# because barcodes are nt unique between surveys, find ones where there is cross over and remove
# as the barcode alone is not enough to discern what survey it is from
non_miss_id <- ID_split[-na.omit(match(missing$Barcode,ID_split$Barcode)),]
non_miss_id <- non_miss_id[non_miss_id$DHS,]
hits <- cbind(match(non_miss_id$Barcode,extr$CDPR50FL$Barcode),
              match(non_miss_id$Barcode,extr$CDPR61FL$HIV_Barcode),
              match(non_miss_id$Barcode,extr$CDPR61FL$VaccSmear_Barcode))
difficults <- non_miss_id[which(apply(hits,1,function(x) length(na.omit(unique(x))))>1),]
ID_split <- ID_split[-match(difficults$ID,ID_split$ID),]

# matching function
match_bar_code <- function(ID_row){

  barcode <- ID_row$Barcode
  grab <- c("hhid","hvidx","hv001","LATNUM","LONGNUM","ADM1NAME","DHSREGNA","SurveyId","Barcode")
  dfret <- data.frame("hhid"=NA,"hvidx"=NA,"hv001"=NA,"LATNUM"=NA,
                         "LONGNUM"=NA,"ADM1NAME"=NA,"DHSREGNA"=NA,
                         "SurveyId"=NA,"Barcode"=barcode,"Notes"="",
                      stringsAsFactors = FALSE)

  # is it a missing child
  if(ID_row$MissingChildren){

    mat <- match(barcode,missing$Barcode)
    dfret$SurveyId <- "CD2013DHS"
    dfret$hv001 = missing$hv001[mat]
    dfret[c("LATNUM","LONGNUM","ADM1NAME","DHSREGNA")] <-
      dat_read$CDGE61FL@data[dat_read$CDGE61FL$DHSCLUST==dfret$hv001,
                             c("LATNUM","LONGNUM","ADM1NAME","DHSREGNA")]
    dfret$Notes <- "hv001 (used for lat/long etc) from MissingChildren_Cluster"

  # else is it a DHS sample
  } else if(ID_row$DHS) {
    # is it from 2007
    mat <- match(barcode,extr$CDPR50FL$Barcode)
    if (!is.na(mat)){
      dfret <- extr$CDPR50FL[mat,grab]
    # or is it from 2013
    } else {
      mat <- match(barcode,extr$CDPR61FL$Barcode)
      if(length(mat)>1) stop("huh")
      if(!is.na(mat)){
        dfret <- extr$CDPR61FL[mat,grab]
      } else {
          dfret$Notes <- "Not in DHS nor MissingChildren"
        }
      }
    }

  return(dfret)
}

# make our meta match data table
ID_split$ROW <- factor(seq_len(nrow(ID_split)), levels=seq_len(nrow(ID_split)))
get <- by(ID_split,ID_split$ROW,match_bar_code) %>% data.table::rbindlist()

# add it and remove extra barcode column save
ID_split <- cbind(ID_split,get)
ID_split <- ID_split[,-which(names(ID_split)=="Barcode")[2]]
saveRDS(ID_split, "analysis/data/derived/clean_meta.rds")

#####################
# Required packages #
#####################



require(data.table)
require(dplyr)
require(pbmcapply)
require(chillR)



##############################################################################################################################################
##############################################################################################################################################



################
# Define Paths #
################



# 1. Input
##########


#Climate
EOBS.dir  = "Input/EOBS_data"

#Photoperiod
photo.dir = "Input"

#Phenology
PEP.dir   = "Output/Merged_file"


# 2. Output
###########

PEP_models.dir = "Output_models"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



## Phenology data
#################

PEP.df <- fread(paste(PEP.dir, "pep_drivers_spring_data.csv", sep="/")) %>% 
  # select columns
  dplyr::select(c(s_id, species, timeseries, year, lat, lon, leaf_out,
                  GDD, GDD.chill.lo, GDD.chill.hi, GDD.chill.DL)) %>% 
  # Mean GDD sums for each timeseries
  group_by(s_id, species, timeseries) %>% 
  filter(n() >= 20) %>% 
  mutate(GDD          = mean(GDD),
         GDD.chill.lo = mean(GDD.chill.lo),
         GDD.chill.hi = mean(GDD.chill.hi),
         GDD.chill.DL = mean(GDD.chill.DL)) %>% 
  ungroup()


## Photoperiod
##############

photo.df = fread(paste(photo.dir, "Photoperiod.csv", sep="/"))


## Import daily climatic datasets from GLDAS
############################################

#define climate variables
list.files(EOBS.dir)
vn <- c("tn_ens_mean_0.1deg_reg_v30.0e.csv", # daily minimum temperature
        "tx_ens_mean_0.1deg_reg_v30.0e.csv") # daily maximum temperature

#create empty list
DataList <- replicate(length(vn),data.frame())
#loop through climate variables
for(i in 1:length(vn)) {
  #read data
  data = fread(paste0(EOBS.dir, "/", vn[i]))
  #add table to list
  DataList[[i]] = data }
#add names to list
names(DataList)=vn



##############################################################################################################################################
##############################################################################################################################################



######################
## Helper functions ##
######################



################################
# Chilling and GDD calculation #
################################


# Chilling units
Chill.fun = function(Temp, threshold_low, threshold_up){
  if(between(Temp,threshold_low, threshold_up)) {CU=1} else {CU=0}
  return(CU)
}

# Forcing units
GDD.fun = function(Temp, threshold_low){
  if(Temp >= threshold_low) {GDD=Temp-threshold_low} else {GDD=0}
  return(GDD)
}


#####################
# Photoperiod model #
#####################


# Function for photoperiod sensitivity with different k and p values
photoperiod_response <- function(DL, k, p) {
  DL_adjusted <- pmax(DL - 6, 0)  # Shift DL so that it starts at 6 hours instead of 0
  f_DL <- (log(1 + k * DL_adjusted) / log(1 + k * (20 - 6)))^p  # Normalize to DL = 20
  f_DL <- pmin(f_DL, 1)  # Cap at 1 from DL = 20 to 24
  return(f_DL)
}


##################
# Chilling model #
##################

# Function for chilling dependency scaling factor
chilling_response <- function(Chilling, requirement) {
  factor <- pmin(Chilling / requirement, 1)  # chilling dependency reaches 1 at X chill days
  return(factor)
}
# Chilling = actual chilling unit 
# requirement = total chilling requirement (in days), i.e., after how many chill days is chilling completely fulfilled?


#####################################
# Function to compute adjusted GDDs #
#####################################

adjust_gdd <- function(GDD, Chilling, requirement=150, include_photo = F, DL = NULL, k = 1, p = 1) {
  f_Chill <- chilling_response(Chilling, requirement)
  
  if (include_photo) {
    if (is.null(DL)) stop("DL must be provided when include_photo is TRUE")
    f_Photo <- photoperiod_response(DL, k, p)
    GDD_adjusted <- GDD * f_Chill * f_Photo
  } else {
    GDD_adjusted <- GDD * f_Chill
  }
  
  return(GDD_adjusted)
}


#################################
# Replace NA with neighbor mean #
#################################


replace_na_with_neighbor_mean <- function(x) {
  if (!is.numeric(x)) return(x)  # Ensure only numeric columns are processed
  
  na_indices <- which(is.na(x))  # Identify NA indices
  known_indices <- which(!is.na(x))  # Identify non-NA values
  
  if (length(known_indices) == 0) {
    return(x)  # If all values are NA, return as is
  }
  
  for (i in na_indices) {
    # Find nearest non-NA before and after
    before_idx <- ifelse(any(known_indices < i), max(known_indices[known_indices < i]), NA)
    after_idx <- ifelse(any(known_indices > i), min(known_indices[known_indices > i]), NA)
    
    if (!is.na(before_idx) & !is.na(after_idx)) {
      # Compute mean of the surrounding values
      x[i] <- mean(c(x[before_idx], x[after_idx]), na.rm = TRUE)
    } else if (!is.na(before_idx)) {
      x[i] <- x[before_idx]  # If NA is at the end, take the previous value
    } else if (!is.na(after_idx)) {
      x[i] <- x[after_idx]  # If NA is at the beginning, take the next value
    }
  }
  
  return(x)
}



##############################################################################################################################################
##############################################################################################################################################



##############################
## Merge datasets into list ##
##############################



# Identifier 1 (all site x year combinations)
PEP.df$site_year = paste0(PEP.df$s_id,"_",PEP.df$year)

# Identifier 2 (all timeseries x year combinations)
PEP.df$ts_yr     = paste0(PEP.df$timeseries,"_",PEP.df$year)
timeseries_year  = unique(PEP.df$ts_yr)

# add PEP data (+plant functional type label) and photoperiod to list
DataList[[3]] = photo.df
DataList[[4]] = PEP.df

rm(photo.df, data, PEP.df)
names(DataList)=c(vn,"photoperiod","PEP")
names(DataList)



#############################################################
## Loop through each observations using parallel computing ##
#############################################################



parallelCalc <- function(timeseries_years){ 

  # Subset input data by time-point
  #################################
  
  #phenology data
  pheno.sub  <- DataList[[4]][which(DataList[[4]]$ts_yr==timeseries_years),]
  
  # Define the length of the period of interest in calendar units
  year1     <- as.character(pheno.sub$year-1)
  year2     <- as.character(pheno.sub$year)
  start_doy <- paste(year1,"-01-01", sep="") 
  end_doy   <- paste(year2,"-06-30", sep="")
  days      <- seq(as.Date(start_doy), as.Date(end_doy), by="days")
  
  #get number of days in year 1 and year 2
  days_year1 <- length(seq(as.Date(start_doy), as.Date(paste(year1,"-12-31", sep="")), by="days"))
  days_year2 <- length(seq(as.Date(paste0(year2,"-01-01")), as.Date(end_doy), by="days"))
  
  
  # ------------------------------------------------------------------------
  # Short Name        ; Long Name                               ; Unit [d-1]
  # ------------------------------------------------------------------------
  # TMIN              ; Minimum Air Temperature                 ; degC
  # TMAX              ; Maximum Air Temperature                 ; degC
  # ------------------------------------------------------------------------
  
  
  #Tmin
  T_min.sub2 <- DataList[[1]][which(DataList[[1]]$site_year==pheno.sub$site_year),]%>%
    dplyr::select(as.character(1:days_year2))
  
  #--------------------------------------
  #Skip timeseries for which there is: 
  if (
    # 1) no climate data for current year and/or
    nrow(T_min.sub2)==0 |
    # 2) no climate data for previous year
    nrow(DataList[[1]][which(DataList[[1]]$s_id==pheno.sub$s_id & DataList[[1]]$year==year1),])==0
  ) {} else {
    #--------------------------------------
    
    T_min.sub1 <- DataList[[1]][which(DataList[[1]]$s_id==pheno.sub$s_id & DataList[[1]]$year==year1),]%>% 
      dplyr::select(as.character(1:days_year1))
    T_min.sub  <- cbind(T_min.sub1, T_min.sub2)
    
    #Tmax
    T_max.sub1 <- DataList[[2]][which(DataList[[2]]$s_id==pheno.sub$s_id & DataList[[2]]$year==year1),]%>% 
      dplyr::select(as.character(1:days_year1))
    T_max.sub2 <- DataList[[2]][which(DataList[[2]]$site_year==pheno.sub$site_year),]%>%
      dplyr::select(as.character(1:days_year2))
    T_max.sub  <- cbind(T_max.sub1, T_max.sub2)
    
    #Photo
    photo.sub1 <- DataList[[3]][which(DataList[[3]]$lat==pheno.sub$lat),][1]%>% 
      dplyr::select(as.character(1:days_year1))
    photo.sub2 <- DataList[[3]][which(DataList[[3]]$lat==pheno.sub$lat),][1]%>%
      dplyr::select(as.character(1:days_year2))
    photo.sub  <- cbind(photo.sub1, photo.sub2)
    
    
    ##############################################################################################################################################
    
    
    ###############################
    # Create table of daily climate
    ###############################
    
    
    # Generate sub-dataframe to store results
    factors.sub <- pheno.sub %>% 
      dplyr::select(s_id,species,timeseries,year,leaf_out)
    factors.sub = as.data.frame(factors.sub) 
    
    #create table
    daily_vals <- data.frame(Year     = lubridate::year(as.Date(days,origin=days[1])),
                             Month    = lubridate::month(as.Date(days,origin=days[1])),
                             Day      = lubridate::day(as.Date(days,origin=days[1])),
                             DOY      = as.numeric(names(T_min.sub)),
                             Tmin     = as.numeric(T_min.sub), 
                             Tmax     = as.numeric(T_max.sub), 
                             DL       = as.numeric(photo.sub)) %>% 
      #replace NAs with neighbor mean
      mutate(across(c(Tmin, Tmax, DL), replace_na_with_neighbor_mean))
    
    #shortest day of year
    solstice = which(daily_vals$DL==min(daily_vals$DL))[1] 
    
    
    ##############################################################################################################################################
    
    
    ######################
    # Add chilling and GDD
    ######################
    
    
    #Get hourly values
    hourly_vals = stack_hourly_temps(daily_vals, latitude=pheno.sub$lat)$hourtemps
    
    # Chilling and Forcing hours
    Chill_days = hourly_vals %>%
      #add Chilling and Forcing units
      mutate(Chill = sapply(Temp, Chill.fun, threshold_low=-10, threshold_up=10),
             GDD   = sapply(Temp, GDD.fun,   threshold_low=5)) %>%
      #summarise daytime hours
      group_by(Year, Month, Day) %>%
      summarise(Chill = mean(Chill, na.rm=T),
                GDD   = mean(GDD, na.rm=T))%>%
      ungroup() %>%
      #order
      dplyr::select(Chill, GDD) 
    
    # define November first as start date of chilling accumulation
    start_row <-  which(daily_vals$Month==11 & daily_vals$Day==1)  
    
    # Cbind
    daily_vals = cbind(daily_vals, Chill_days) %>%
      mutate(
        #continous count from Jan 1
        id = row_number(),
        # Sum chilling days from 1 November
        Chill_sum = ifelse(id < start_row, 0, cumsum(ifelse(id >= start_row, Chill, 0))),
        # Add chilling/photoperiod-adjusted GDDs
        GDD.chill.lo = adjust_gdd(GDD=GDD, Chilling=Chill_sum, requirement=75),
        GDD.chill.hi = adjust_gdd(GDD=GDD, Chilling=Chill_sum, requirement=150),
        GDD.chill.DL = adjust_gdd(GDD=GDD, Chilling=Chill_sum, requirement=150, include_photo = T, DL = DL)) %>% 
      #order
      dplyr::select(id, Year, Month, Day, DOY,
                    GDD, GDD.chill.lo, GDD.chill.hi, GDD.chill.DL) %>% 
      #calculate GDD sums from winter solstice
      filter(!id < solstice) %>% 
      mutate(GDD          = cumsum(GDD),
             GDD.chill.lo = cumsum(GDD.chill.lo),
             GDD.chill.hi = cumsum(GDD.chill.hi),
             GDD.chill.DL = cumsum(GDD.chill.DL)) %>% 
      #delete days of previous year
      filter(!Year == year1) 
    
    
    ##############################################################################################################################################
    
    
    ########################
    # Predict leaf-out dates
    ########################
    
    
    factors.sub = factors.sub %>% 
      mutate(leaf_out_GDD          = which(daily_vals$GDD          >= pheno.sub$GDD)[1],
             leaf_out_GDD.chill.lo = which(daily_vals$GDD.chill.lo >= pheno.sub$GDD.chill.lo)[1],
             leaf_out_GDD.chill.hi = which(daily_vals$GDD.chill.hi >= pheno.sub$GDD.chill.hi)[1],
             leaf_out_GDD.chill.DL = which(daily_vals$GDD.chill.DL >= pheno.sub$GDD.chill.DL)[1])
    
    
    ##############################################################################################################################################
    
    
    # Safe the table   
    return(factors.sub)
    
  }
}



##############################################################################################################################################
##############################################################################################################################################



#initialize the loop
outputlist <- pbmclapply(timeseries_year, parallelCalc, mc.cores=48, mc.preschedule=T)
climate.factors.df <- rbindlist(outputlist)



##############################################################################################################################################
##############################################################################################################################################



###################
## Safe the data ##
###################



#Safe table
write.csv(climate.factors.df, paste(PEP_models.dir, "pep_model_predictions.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



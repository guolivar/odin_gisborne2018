#'---
#'
#'---


#########################################################################
# Parsing, cleaning ODIN-SD offline data with intermittent clock updates.
# Developed for the project "Air Quality Measurements in Gisborne - 2018"
# Authors:
#         Gustavo Olivares
#########################################################################

# Prepare libraries
library(librarian) # To more flexibly manage packages
shelf(readr,
      openair)
# Set the data folder ####
data_path <- "~/data/ODIN_SD/Gisborne/WORKING/"
folders_list <- dir(data_path,pattern = '00')
# Define time average for output
tavg <- '1 hour' # This is for TimeAverage
txt_tavg <- '1hr' # This is for file naming
# Load data ####
# Cycle through the folders to work with all the DATA.TXT files
for (i in (1:length(folders_list))){
  folder <- folders_list[i]
  print(folder)
  odin.data <- readr::read_delim(paste0(data_path,folder,"/DATA.TXT"),
                                 delim = ';',
                                 skip = 1,
                                 col_names = c('framelength',
                                               'PM1',
                                               'PM2.5',
                                               'PM10',
                                               'PM1x',
                                               'PM2.5x',
                                               'PM10x',
                                               'GasSN',
                                               'Gasppm',
                                               'GasT',
                                               'Gas2mV',
                                               'Temperature',
                                               'RH',
                                               'ODINid',
                                               'ODINsn',
                                               'RTCdate',
                                               'RTCtime',
                                               'GSMdate',
                                               'GSMtime',
                                               'RTCdate2',
                                               'RTCtime2'))
  # Pad missing GSM date-time with RTC date-time for simplicity
  odin.data$GSMtime[is.na(odin.data$GSMtime)] <- odin.data$RTCtime[is.na(odin.data$GSMtime)]
  odin.data$GSMdate[is.na(odin.data$GSMdate)] <- odin.data$RTCdate[is.na(odin.data$GSMdate)]
  
  # Construct POSIX objects for RTC and GSM timestamps
  odin.data$RTCtimestamp <- as.POSIXct(paste(odin.data$RTCdate,odin.data$RTCtime),tz='UTC')
  odin.data$GSMtimestamp <- as.POSIXct(paste(odin.data$GSMdate,odin.data$GSMtime),tz='UTC')
  # Find the correction from RTC to GSM timestamps
  time_correction.all <- as.numeric(difftime(odin.data$GSMtimestamp,odin.data$RTCtimestamp,units = 'secs'))
  time_correction.all[time_correction.all==0] <- NA
  time_diff <- mean(time_correction.all,na.rm = TRUE)
  # Calculate the "real" timestamp for the records
  odin.data$date <- odin.data$RTCtimestamp + time_diff
  odin.data <- odin.data[,c('date',
                            'PM1',
                            'PM2.5',
                            'PM10',
                            'Temperature',
                            'RH',
                            'ODINsn')]
  # Construct the ALLDATA frame
  if (i == 1){
    # This is the first iteration so we just copy the "odin.data" dataframe
    all.data <- odin.data
    all.data.tavg <- timeAverage(odin.data,avg.time = tavg)
    all.data.tavg$ODINsn <- odin.data$ODINsn[1]
  } else {
    # We already have "all.data" so we need to append the current "odin.data"
    all.data <- rbind(all.data,odin.data)
    tmp1 <- timeAverage(odin.data,avg.time = tavg)
    tmp1$ODINsn <- odin.data$ODINsn[1]
    all.data.tavg <- rbind(all.data.tavg,tmp1)
    # Remove all.data to clean for next iteration
    rm(odin.data)
  }
}

# Get devices locations ####
odin_locations <- readr::read_delim("odin_locations.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
nsites <- length(odin_locations$serialn)
all.data$lon <- 0
all.data$lat <- 0
all.data.tavg$lon <- 0
all.data.tavg$lat <- 0
for (j in (1:nsites)){
  print(odin_locations$serialn[j])
  loc_id <- (substr(all.data$ODINsn,3,6)==substr(odin_locations$serialn[j],7,11))
  all.data$lon[loc_id] <- odin_locations$lon[j]
  all.data$lat[loc_id] <- odin_locations$lat[j]
  loc_id.tavg <- (substr(all.data.tavg$ODINsn,3,6)==substr(odin_locations$serialn[j],7,11))
  all.data.tavg$lon[loc_id.tavg] <- odin_locations$lon[j]
  all.data.tavg$lat[loc_id.tavg] <- odin_locations$lat[j]
}

ggplot(data = all.data.tavg, aes(x=date)) +
  geom_line(aes(y=PM2.5,colour = ODINsn))
# This shows that ODIN-0005 malfunctioned so those data will be removed
remove_idx <- (all.data$ODINsn == "SD0005")
all.data[remove_idx,] <- NA
all.data <- na.exclude(all.data)
remove_idx <- (all.data.tavg$ODINsn == "SD0005")
all.data.tavg[remove_idx,] <- NA
all.data.tavg <- na.exclude(all.data.tavg)
# Save the "all.data" dataframe
save(all.data,file = 'alldata.RData')
save(all.data.tavg,file = paste0('alldata',txt_tavg,'.RData'))

#'---
#'title: "ODIN-SD Summary"
#'author: Gustavo Olivares
#'---

#'This script parses and cleans the ODIN-SD data from their memory cards to then calculate summary statistics 
#'and plot draft maps of PM~10~ and PM~2.5~
#'

#'
#' ## Prepare libraries
#'  
library(librarian) # To more flexibly manage packages
shelf(readr,
      openair,
      automap,
      raster,
      gstat,
      sp,
      rgdal,
      ggmap,
      ggplot2,
      scales)
#' ## Set constants
data_path <- path.expand("~/data/ODIN_SD/Gisborne/WORKING/")
folders_list <- dir(data_path,pattern = '00')
# Define time average for output
tavg <- '1 day'
#' ## Load BAM data
#'
#'  

bam.data <- read_csv(paste0(data_path,"BAM data PM10.csv"),
                     col_types = cols(Day = col_character(),
                                      SiteName = col_skip(),
                                      Time = col_character()))
bam.data$ODINsn <- 'GDC'
# Move the date to GMT for consistency
bam.data$date <- as.POSIXct(paste(bam.data$Day,bam.data$Time),format = "%d/%m/%Y %I:%M:%S %p", tz= "GMT") - 12*3600
bam.data.tavg <- timeAverage(bam.data,avg.time = tavg)
bam_lat <- -(38 + 39/60 + 37.20961/3600)
bam_lon <- 178 + 00/60 + 42.5312/3600

#' ## Load ODIN data
#'
#'  Cycle through the folders to work with all the DATA.TXT files
#'  
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


#'
#'  Get devices locations
#'  
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
# Subset for the campaign
start_date <-as.POSIXct("2018-06-22 12:00:00",tz="UTC")
end_date <-as.POSIXct("2018-08-16 12:00:00",tz="UTC")
nminutes <- difftime(end_date,start_date,units = 'min')
all.data <- subset(all.data, (date >= start_date) & (date <= end_date))
all.data.tavg <- subset(all.data.tavg, (date >= start_date) & (date <= end_date))

#'
#' ## Summary statistics
#'
#'Note that the statistics are calculated on 10 minutes averages and that the campaign
#'went from June 22nd 12:00 UTC until August 16th 12:00 UTC (55 days in total).
#'
#'The units for the parameters are:
#' * PM~2.5~ [$mu$g/m^3^]
#' * PM~10~ [$mu$g/m^3^]
#' * Temperature [celsius]
#' * RH [%]
#' 

# Calculate the summary table for each unit
summary_mean <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data.tavg, FUN = mean)
summary_max <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data.tavg, FUN = max)
summary_min <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data.tavg, FUN = min)
summary_sd <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data.tavg, FUN = sd)
summary_N <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data.tavg, FUN = length)
summary_pct <- summary_N
summary_pct[,2:5] <- format(100 * (summary_N[,2:5] / (as.numeric(nminutes) / 10)),digits = 2)

#' 
#' ## Average concentrations
#' 
print(format(summary_mean,digits = 1))
#' 
#' ## Maximum concentrations
#' 
print(format(summary_max,digits = 1))
#' 
#' ## Minimum concentrations
#' 
print(format(summary_min,digits = 1))
#' 
#' ## Standard deviation
#' 
print(format(summary_sd,digits = 1))
#' 
#' ## Data availability [%]
#' 
print(format(summary_pct,digits = 1))

#'
#' ### Malfunctioning ODIN
#'  The previous summaries show that ODIN-0005 malfunctioned so it will be removed from the analysis.
#'   Also, after the QA process we identified ODIN-0028 as malfunctional so we also removed that unit from the analysis.
remove_idx <- ((all.data$ODINsn == "SD0005") | (all.data$ODINsn == "SN0028"))
all.data[remove_idx,] <- NA
all.data <- na.exclude(all.data)
remove_idx <- ((all.data.tavg$ODINsn == "SD0005") | (all.data.tavg$ODINsn == "SN0028"))
all.data.tavg[remove_idx,] <- NA
all.data.tavg <- na.exclude(all.data.tavg)

# Merge BAM dataset with TAVG
bam.data.tavg <- subset(bam.data.tavg, (date >= start_date) & (date <= end_date))

ggplot(data = all.data.tavg, aes(x=date)) +
  geom_line(data = subset(all.data.tavg,ODINsn == "SD0011"),aes(y=PM10,colour = "G6")) +
  geom_line(data = subset(all.data.tavg,ODINsn == "SD0025"),aes(y=PM10,colour = "G1")) +
  geom_line(data = subset(all.data.tavg,ODINsn == "SD0010"),aes(y=PM10,colour = "G16")) +
  geom_line(data = bam.data.tavg, aes(x=date,y=PM10,colour = "GDC"),size = 2)


#'
#' ## Time series
#' 

#'
#' ### PM~2.5~
#' 
ggplot(data = all.data.tavg, aes(x=date)) +
  geom_line(aes(y=PM2.5,colour = ODINsn))
#'
#' ### PM~10~
#' 
ggplot(data = all.data.tavg, aes(x=date)) +
  geom_line(aes(y=PM10,colour = ODINsn))
#'
#' ### Temperature
#' 
ggplot(data = all.data.tavg, aes(x=date)) +
  geom_line(aes(y=Temperature,colour = ODINsn))
#'
#' ### Relative Humidity
#' 
ggplot(data = all.data.tavg, aes(x=date)) +
  geom_line(aes(y=RH,colour = ODINsn))


#'
#' ## Average Maps
#' 
# Some useful constants
proj4string_NZTM <- CRS('+init=epsg:2193')
proj4string_latlon <- CRS('+init=epsg:4326')
# Assign coordinates to the dataframe
summary_mean_map <- aggregate(cbind(Temperature, RH, PM2.5, PM10, lon, lat) ~ODINsn, all.data.tavg, FUN = mean)
coordinates(summary_mean_map) <- ~ lon + lat
proj4string(summary_mean_map) <- proj4string_latlon

# Get the basemap
centre_lat <- mean(summary_mean_map$lat)
centre_lon <- mean(summary_mean_map$lon)
ca <- get_googlemap(
  c(lon=centre_lon,lat=centre_lat),
  zoom=13,
  scale=2,
  color="bw",
  key = "AIzaSyACi3pNvPQTxZWx5u0nTtke598dPqdgySg")

#' ### PM~2.5~
ggmap(ca) + 
  geom_label(data = odin_locations, aes(x=lon,y=lat,label = siteid), nudge_x = -0.004) +
  geom_point(data=as.data.frame(summary_mean_map),aes(x=lon,y=lat,colour = PM2.5),size = 5) +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map$PM2.5)),
                          name = "PM2.5", oob=squish)

#' ### PM~10~
ggmap(ca) +  
  geom_label(data = odin_locations, aes(x=lon,y=lat,label = siteid), nudge_x = -0.004) +
  geom_label(aes(x=bam_lon,y=bam_lat,label = 'BAM'),size = 5, nudge_y = 0.004) +
  geom_point(data = as.data.frame(summary_mean_map),aes(x=lon,y=lat,colour = PM10),size = 5) +
  geom_point(aes(x=bam_lon,y=bam_lat,colour = mean(bam.data.tavg$PM10)),size = 5) +
  geom_point(aes(x=bam_lon,y=bam_lat),size = 1.3, colour = 'black') +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map$PM10)),
                          name = "PM10", oob=squish)

#' ### PM~coarse~
ggmap(ca) +  
  geom_label(data = odin_locations, aes(x=lon,y=lat,label = siteid), nudge_x = -0.004) +
  geom_point(data=as.data.frame(summary_mean_map),aes(x=lon,y=lat,colour = PM10 - PM2.5),size = 5) +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map$PM10 - summary_mean_map$PM2.5)),
                          name = "PMcoarse", oob=squish)



# Save the "all.data" dataframe
save(all.data,file = 'alldata.RData')
save(all.data.tavg,file = 'alldataTAVG.RData')
data.output <- all.data.tavg[,c('date','PM2.5','PM10','Temperature','RH','ODINsn')]
write_excel_csv(data.output,paste0('./gisborne_',tavg,'.csv'))

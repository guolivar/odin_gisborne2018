#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
# Define server logic required to draw a histogram

shinyServer(function(input, output) {
# Read the secrets
secret_hologram <- read_delim("./secret_hologram.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
# Get the devices ID
base_url <- "https://dashboard.hologram.io/api/1/devices?"
tag <- "gisborne"
built_url <- paste0(base_url,
                    "orgid=",secret_hologram$orgid,"&",
                    "tagname=",tag,"&",
                    "apikey=",secret_hologram$apikey)
req1 <- curl_fetch_memory(built_url)
jreq1 <- fromJSON(rawToChar(req1$content))$data
nsites <- length(jreq1)
curr_data <- data.frame(deviceid = (1:nsites),ODIN = NA)
for (i in (1:nsites)){
  curr_data$deviceid[i] <- jreq1[[i]]$id
  curr_data$ODIN[i] <- jreq1[[i]]$name
}
curr_data$Last_reading <- as.POSIXct('2018-05-01 00:00:00',tz='UTC')
# Get the latest measurements
base_url <- "https://dashboard.hologram.io/api/1/csr/rdm?"

i <- 1
for (i in (1:nsites)){
  built_url <- paste0(base_url,
                      "deviceid=",curr_data$deviceid[i],"&",
                      "limit=1&",
                      "apikey=",secret_hologram$apikey)
  req2 <- curl_fetch_memory(built_url)
  jreq2 <- fromJSON(rawToChar(req2$content))$data
  payload <- fromJSON(rawToChar(base64decode(fromJSON(jreq2[[1]]$data)$data)))
  curr_data$Last_reading[i] <- as.POSIXct(jreq2[[1]]$logged,tz='UTC')
}

curr_data$delay <- floor(difftime(Sys.time(),curr_data$Last_reading, units = 'secs'))

curr_data$mask <- as.numeric(curr_data$delay < 3600)
reboot_odins <- subset(curr_data,mask == 0)
output$table <- DT::renderDataTable({
  DT::datatable({
    dead_odins <- length(reboot_odins$deviceid)
    if (dead_odins > 0){
      out <- reboot_odins[,c('ODIN','Last_reading')]
    } else{
      out <- data.frame(message = "All ODIN units are happy!")
    }
    out
    },
    options = list(pageLength = 18))
  })

# Get devices locations
#proj4string <- CRS('+init=epsg:2193')

odin_locations <- read_delim("./odin_locations.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
curr_data$lat <- NA
curr_data$lon <- NA
for (i in (1:nsites)){
  loc_id <- which(substr(odin_locations$serialn,7,11)==substr(curr_data$ODIN[i],6,9))
  curr_data$lon[i] <- odin_locations$lon[loc_id]
  curr_data$lat[i] <- odin_locations$lat[loc_id]
}
centre_lat <- mean(curr_data$lat)
centre_lon <- mean(curr_data$lon)

# Now get datafor last 12 hours and calculate average for mapping
# Cycle through each deviceID and calculate the 12hr average up to now
curr_data$PM1 <- NA
curr_data$PM2.5 <- NA
curr_data$PM10 <- NA
curr_data$Temperature <- NA
curr_data$RH <- NA
max_nmeas <- 60*12
ndev <- length(curr_data$deviceid)
# t_start is 12 hours before now
t_start <- floor(as.numeric(Sys.time()-12*3600))

base_url <- "https://dashboard.hologram.io/api/1/csr/rdm?"
for (i_dev in (1:ndev)){
  built_url <- paste0(base_url,
                      "deviceid=",curr_data$deviceid[i_dev],"&",
                      "limit=",max_nmeas,"&",
                      "timestart=",t_start,"&",
                      "orgid=",secret_hologram$orgid,"&",
                      "apikey=",secret_hologram$apikey)
  req2 <- curl_fetch_memory(built_url)
  jreq2 <- fromJSON(rawToChar(req2$content))$data
  
  ndata <- length(jreq2)
  c_data <- data.frame(id = (1:ndata))
  c_data$PM1 <- NA
  c_data$PM2.5 <- NA
  c_data$PM10 <- NA
  c_data$PMc <- NA
  c_data$GAS1 <- NA
  c_data$Tgas1 <- NA
  c_data$GAS2 <- NA
  c_data$Temperature <- NA
  c_data$RH <- NA
  c_data$timestamp <- NA
  if (ndata < 1){
    next
  }
  for (i in (1:ndata)){
    payload <- fromJSON(rawToChar(base64decode(fromJSON(jreq2[[i]]$data)$data)))
    # {"PM1":4,"PM2.5":6,"PM10":6,"GAS1":-999,"Tgas1":0,"GAS2":204,"Temperature":7.35,"RH":80.85}
    c_data$PM1[i] <- payload[1]
    c_data$PM2.5[i] <- payload[2]
    c_data$PM10[i] <- payload[3]
    c_data$PMc[i] <- payload[3] - payload[2]
    c_data$GAS1[i] <- payload[4]
    c_data$Tgas1[i] <- payload[5]
    c_data$GAS2[i] <- payload[6]
    c_data$Temperature[i] <- payload[7]
    c_data$RH[i] <- payload[8]
    c_data$timestamp[i] <- jreq2[[i]]$logged
  }
  curr_data$PM1[i_dev] <- mean(c_data$PM1,na.rm = TRUE)
  curr_data$PM2.5[i_dev] <- mean(c_data$PM2.5,na.rm = TRUE)
  curr_data$PM10[i_dev] <- mean(c_data$PM10,na.rm = TRUE)
  curr_data$Temperature[i_dev] <- mean(c_data$Temperature,na.rm = TRUE)
  curr_data$RH[i_dev] <- mean(c_data$RH,na.rm = TRUE)
}

curr_data$mask[is.na(curr_data$PM2.5)] <- 0
curr_data$PM2.5[is.na(curr_data$PM2.5)] <- 0

output$plot1 <- renderPlot({
  cmap <- get_map(c(centre_lon,centre_lat),zoom=13,scale = 2)
  ggmap(cmap) +
    geom_point(data = curr_data,
               aes(x=lon,y=lat,colour=PM2.5),
               alpha=curr_data$mask,
               size=15) +
    geom_text(data=curr_data,aes(x=lon,y=lat,label=substring(ODIN,1,9)),colour = 'red') +
    ggtitle("PM2.5 average for the last 12 hours") +
    scale_colour_gradient(low="white", high="red")
  },width = 1024,height = 1024)

})

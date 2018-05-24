# Sample data ... fetch data and play
library(readr)
library(RJSONIO)
library(curl)
library(base64enc)
library(ggplot2)
library(zoo)

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
devices <- data.frame(deviceid = (1:nsites),ODIN = NA)
for (i in (1:nsites)){
  devices$deviceid[i] <- jreq1[[i]]$id
  devices$ODIN[i] <- jreq1[[i]]$name
}

# Get the last X measurements
nmeas <- 60*12
base_url <- "https://dashboard.hologram.io/api/1/csr/rdm?"
ndev <- length(devices$deviceid)
for (i_dev in (1:ndev)){
    n_steps <- ceiling(nmeas/1000)
    for (step in (1:n_steps)){
      if (step == 1){
        built_url <- paste0(base_url,
                            "deviceid=",devices$deviceid[i_dev],"&",
                            "limit=",nmeas,"&",
                            "timestart=1526249648&",
                            "orgid=",secret_hologram$orgid,"&",
                            "apikey=",secret_hologram$apikey)
        req2 <- curl_fetch_memory(built_url)
        jreq2 <- fromJSON(rawToChar(req2$content))$data
      } else {
        built_url <- paste0(base_url,
                            "deviceid=",devices$deviceid[i_dev],"&",
                            "limit=",nmeas,"&",
                            "timestart=1526249648&",
                            "startat=",startat,"&",
                            "orgid=",secret_hologram$orgid,"&",
                            "apikey=",secret_hologram$apikey)
        req2 <- curl_fetch_memory(built_url)
        jreq2 <- append(jreq2,fromJSON(rawToChar(req2$content))$data)
      }
      last <- length(jreq2)
      startat <- jreq2[[last]]$id
    }
    ndata <- length(jreq2)
  if (ndata > 1){
    c_data <- data.frame(id = (1:ndata))
    c_data$serialn <- devices$ODIN[i_dev]
    c_data$device <- devices$deviceid[i_dev]
    c_data$PM1 <- NA
    c_data$PM2.5 <- NA
    c_data$PM10 <- NA
    c_data$PMc < NA
    c_data$GAS1 <- NA
    c_data$Tgas1 <- NA
    c_data$GAS2 <- NA
    c_data$Temperature <- NA
    c_data$RH <- NA
    c_data$timestamp <- NA
    
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
    if (i_dev == 1){
      all_data <- c_data
    } else {
      all_data <- rbind(all_data,c_data)
    }
    rm(c_data)
  }
}

all_data$date <- as.POSIXct(all_data$timestamp,tz="UTC")
ggplot(data = all_data,aes(x=date)) +
  geom_line(aes(y=rollmean(PM2.5,60,fill = NA),colour=serialn))

select_SN <- devices$ODIN[grep("0019",devices$ODIN)]
plot_data <- subset(all_data,serialn==select_SN)
ggplot(data = plot_data,aes(x=date)) +
  geom_line(aes(y=rollmean(PM1,60,na.pad=TRUE),colour='PM1'))+
  geom_line(aes(y=rollmean(PM2.5,60,na.pad=TRUE),colour='PM2.5'))+
  geom_line(aes(y=rollmean(PM10,60,na.pad=TRUE),colour='PM10')) +
  geom_line(aes(y=rollmean(Temperature * 5,60,na.pad=TRUE),colour='T * 10')) +
  ggtitle(select_SN)

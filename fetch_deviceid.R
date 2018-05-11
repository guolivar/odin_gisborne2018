#' ---
#' title: "Fetch IDs for the ODINs deployed in <TAG>"
#' author: "Gustavo Olivares"
#' ---

library(readr)
library(RJSONIO)
library(curl)

# Read the secrets
secret_hologram <- read_delim("/data/GusData/repositories/odin_gisborne2018/secret_hologram.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
base_url <- "https://dashboard.hologram.io/api/1/devices?"
tag <- "gisborne"
built_url <- paste0(base_url,
                    "orgid=",secret_hologram$orgid,"&",
                    "tagname=",tag,"&",
                    "apikey=",secret_hologram$apikey)
req1 <- curl_fetch_memory(built_url)
jreq <- fromJSON(rawToChar(req1$content))$data
devices <- data.frame(deviceid = (1:length(jreq)))
for (i in (1:length(jreq))){
  devices$deviceid[i] <- jreq[[i]]$id
}


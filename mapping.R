#########################################################################
# Taking cleaned, averaged ODIN data and generate maps and animations.
# Developed for the project "Air Quality Measurements in Gisborne - 2018"
# Authors:
#         Gustavo Olivares
#########################################################################
##### Load relevant packages #####
library(librarian) # To more flexibly manage packages
shelf(readr,
      reshape2,
      automap,
      raster,
      gstat,
      sp,
      rgdal,
      ggmap,
      ggplot2,
      scales,
      gstat,
      RNetCDF,
      RJSONIO,
      curl,
      base64enc,
      zoo,
      openair,
      stringi,
      viridis,
      dplyr,
      RColorBrewer,
      purrr,
      magick)
# Set the data folder ####
data_path <- "~/data/ODIN_SD/Gisborne/WORKING/"
plot_path <- "~/data/ODIN_SD/Gisborne/"

# Load data
load('alldata.RData')
load('alldataTAVG.RData')
# Some useful constants
proj4string_NZTM <- CRS('+init=epsg:2193')
proj4string_latlon <- CRS('+init=epsg:4326')
# Assign coordinates to the dataframes
coordinates(all.data) <- ~ lon + lat
proj4string(all.data) <- proj4string_latlon
coordinates(all.data.tavg) <- ~ lon + lat
proj4string(all.data.tavg) <- proj4string_latlon

centre_lat <- mean(all.data.tavg$lat)
centre_lon <- mean(all.data.tavg$lon)

## Prepare the map to plot animations #####

# Get the basemap
ca <- get_map(
  c(lon=centre_lon,lat=centre_lat),
  zoom=13,crop=T,
  scale="auto",color="bw",source="google",
  maptype="terrain") # can change to terrain

# Re-project to NZTM #####
all.data.tavg <- spTransform(all.data.tavg,proj4string_NZTM)
print("Starting the kriging")

#Setting the  prediction grid properties #####
cellsize <- 100 #pixel size in projection units (NZTM, i.e. metres)
min_x <- all.data.tavg@bbox[1,1] - cellsize#minimun x coordinate
min_y <- all.data.tavg@bbox[2,1] - cellsize #minimun y coordinate
max_x <- all.data.tavg@bbox[1,2] + cellsize #mximum x coordinate
max_y <- all.data.tavg@bbox[2,2] + cellsize #maximum y coordinate

x_length <- max_x - min_x #easting amplitude
y_length <- max_y - min_y #northing amplitude

ncol <- round(x_length/cellsize,0) #number of columns in grid
nrow <- round(y_length/cellsize,0) #number of rows in grid

grid <- GridTopology(cellcentre.offset=c(min_x,min_y),cellsize=c(cellsize,cellsize),cells.dim=c(ncol,nrow))

#Convert GridTopolgy object to SpatialPixelsDataFrame object. #####
grid <- SpatialPixelsDataFrame(grid,
                               data=data.frame(id=1:prod(ncol,nrow)),
                               proj4string=CRS('+init=epsg:2193'))


# Get rid of NA containing rows
all.data.tavg <- subset(all.data.tavg,!is.na(PM2.5))
all_dates <- sort(unique(all.data.tavg$date))
valid_dates <- FALSE * (1:length(all_dates))
# limits for colorscales #####
cmin <- min(all.data.tavg$PM2.5)
cmax <- max(all.data.tavg$PM2.5) * 0.5
## Interpolate and plot #####
ndates <- length(all_dates)
breaks <- as.numeric(quantile((1:ndates),c(0,0.5,1), type = 1))
nbreaks <- length(breaks)
i <- 0
for (d_slice in (1:ndates)){
  c_data <- subset(all.data.tavg,subset = (date==all_dates[d_slice]))
  
  if (length(unique(c_data$ODINsn))<2){
    next
  }
  valid_dates[d_slice] <- TRUE
  #  surf.krig <- autoKrige(pm2.5 ~ 1,data=c_data,new_data = grid, input_data=c_data)
  #  surf.krig$krige_output$timestamp <-d_slice
  #  proj4string(surf.krig$krige_output) <- CRS('+init=epsg:2193')
  
  # surf.idw <- idw(PM2.5 ~ 1,newdata = grid, locations = c_data, idp = 1,na.action = na.omit)
  # surf.idw$timestamp <-d_slice
  # proj4string(surf.idw) <- CRS('+init=epsg:2193')
  
  surf.idw2 <- idw(PM2.5 ~ 1,newdata = grid, locations = c_data, idp = 2)
  surf.idw2$timestamp <-d_slice
  proj4string(surf.idw2) <- CRS('+init=epsg:2193')
  
  if (i==0){
    #    x_data <- surf.krig$krige_output@data
    #    x_bbox <- surf.krig$krige_output@bbox
    #    x_coords <- surf.krig$krige_output@coords
    #    x_coords.nrs <- c(1,2)
    #    to_rast.krig <- surf.krig$krige_output
    #    r0.krig <- rasterFromXYZ(cbind(to_rast.krig@coords,to_rast.krig@data$var1.pred))
    #    crs(r0.krig) <- '+init=epsg:2193'
    #    raster_cat.krig <- r0.krig
    
    # to_rast.idw <- surf.idw
    # r0.idw <- rasterFromXYZ(cbind(surf.idw@coords,surf.idw$var1.pred))
    # crs(r0.idw) <- '+init=epsg:2193'
    # raster_cat.idw<- r0.idw
    
    to_rast.idw2 <- surf.idw2
    r0.idw2 <- rasterFromXYZ(cbind(surf.idw2@coords,surf.idw2$var1.pred))
    crs(r0.idw2) <- '+init=epsg:2193'
    raster_cat.idw2<- r0.idw2
    i <- 1
  }
  else {
    #    x_data <- rbind(x_data,surf.krig$krige_output@data)
    #    x_coords <- rbind(x_coords,surf.krig$krige_output@coords)
    #    to_rast.krig <- surf.krig$krige_output
    #    r0.krig <- rasterFromXYZ(cbind(to_rast.krig@coords,to_rast.krig@data$var1.pred))
    #    crs(r0.krig) <- '+init=epsg:2193'
    #    raster_cat.krig <- addLayer(raster_cat.krig,r0.krig)
    
    # to_rast.idw <- surf.idw
    # r0.idw <- rasterFromXYZ(cbind(surf.idw@coords,surf.idw$var1.pred))
    # names(r0.idw) <- as.character(all_dates[d_slice])
    # crs(r0.idw) <- '+init=epsg:2193'
    # raster_cat.idw<- addLayer(raster_cat.idw,r0.idw)
    
    to_rast.idw2 <- surf.idw2
    r0.idw2 <- rasterFromXYZ(cbind(surf.idw2@coords,surf.idw2$var1.pred))
    names(r0.idw2) <- as.character(all_dates[d_slice])
    crs(r0.idw2) <- '+init=epsg:2193'
    raster_cat.idw2<- addLayer(raster_cat.idw2,r0.idw2)
  }
  # rtp <- rasterToPolygons(projectRaster(r0.idw,crs = "+proj=longlat +datum=WGS84"))
  rtp2 <- rasterToPolygons(projectRaster(r0.idw2,crs = "+proj=longlat +datum=WGS84"))
  points <- data.frame(spTransform(c_data,CRS('+init=epsg:4326')))
  
  # Build the animation
  # map_out <- ggmap(ca) + geom_polygon(data = rtp,aes(x = long, y = lat, group = group, 
  #                                                    fill = rep(rtp[[1]], each = 5)), 
  #                                     size = 0, 
  #                                     alpha = 0.85) +
  #   scale_fill_gradient(low="white", high="red",limits=c(0, cmax), name = "PM2.5", oob=squish) +
  #   geom_point(data=points,aes(x=lon,y=lat),colour = "black") +
  #   ggtitle(paste(as.character(all_dates[d_slice]+12*3600),"NZST"))
  # ggsave(filename=paste0(data_path,'../idw/',format(all_dates[d_slice]+12*3600,format = "%Y-%m-%d %H:%M"),'.png'), plot=map_out, width=6, height=6, units = "in")
  
  map_out <- ggmap(ca) + geom_polygon(data = rtp2,aes(x = long, y = lat, group = group, 
                                                      fill = rep(rtp2[[1]], each = 5)), 
                                      size = 0, 
                                      alpha = 0.8) +
    scale_fill_gradient(low="white", high="red",limits=c(0, cmax), name = "PM2.5", oob=squish) +
    geom_point(data=points,aes(x=lon,y=lat),colour = "black") +
    ggtitle(paste(as.character(all_dates[d_slice]+12*3600),"NZST"))
  ggsave(filename=paste0(data_path,'../idw2/',format(all_dates[d_slice]+12*3600,format = "%Y-%m-%d %H:%M"),'.png'),
         plot=map_out,
         width=6,
         height=6,
         units = "in")
  
}
# save('raster_cat.idw',file = paste0(data_path,'raster_cat.idw.RData'))
save('raster_cat.idw2',file = paste0('raster_cat.idw2.Rdata'))

print("Done with interpolating ...")

# raster_cat_idw_LL <- projectRaster(raster_cat.idw,crs = "+proj=longlat +datum=WGS84")
raster_cat_idw2_LL <- projectRaster(raster_cat.idw2,crs = "+proj=longlat +datum=WGS84")
save(list = c('raster_cat_idw2_LL'),file = paste0("raster_odin_LL_IDW.RData"))

# Plot time series ####
plot_tseries <- ggplot(data.frame(all.data.tavg),aes(x=date)) +
  geom_line(aes(y=PM2.5,colour=ODINsn))
ggsave(filename = paste0(data_path,
                         't_series_',
                         format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
                         format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
                         ".png"),
       plot = plot_tseries,
       width = 12,
       height = 6,
       units = 'in')

# Write NetCDF files ####
# # IDW
# lat_dim <- unique(coordinates(raster_cat_idw_LL)[,2])
# lon_dim <- unique(coordinates(raster_cat_idw_LL)[,1])
# tim_dim <- all_dates[valid_dates]
# nc.idw <- create.nc("odin_idw.nc")
# # Dimensions specifications
# dim.def.nc(nc.idw, "time", unlim=TRUE)
# dim.def.nc(nc.idw, "latitude",length(lat_dim))
# dim.def.nc(nc.idw, "longitude",length(lon_dim))
# # Variable specifications
# var.def.nc(nc.idw,"time","NC_INT","time")
# att.put.nc(nc.idw,"time","units","NC_CHAR","seconds since 1970-01-01 00:00:0.0")
# att.put.nc(nc.idw,"time","long_name","NC_CHAR","time")
# 
# var.def.nc(nc.idw,"latitude","NC_FLOAT","latitude")
# att.put.nc(nc.idw,"latitude","units","NC_CHAR","degrees_north")
# att.put.nc(nc.idw,"latitude","long_name","NC_CHAR","latitude")
# att.put.nc(nc.idw,"latitude","standard_name","NC_CHAR","latitude")
# 
# var.def.nc(nc.idw,"longitude","NC_FLOAT","longitude")
# att.put.nc(nc.idw,"longitude","units","NC_CHAR","degrees_east")
# att.put.nc(nc.idw,"longitude","long_name","NC_CHAR","longitude")
# att.put.nc(nc.idw,"longitude","standard_name","NC_CHAR","longitude")
# 
# var.def.nc(nc.idw,"pm2p5","NC_FLOAT",c("longitude","latitude","time"))
# att.put.nc(nc.idw,"pm2p5","units","NC_CHAR","ug m**-3")
# att.put.nc(nc.idw,"pm2p5","long_name","NC_CHAR","Mass concentration of PM2.5 ambient aerosol particles in air")
# att.put.nc(nc.idw,"pm2p5","standard_name","NC_CHAR","mass_concentration_of_pm2p5_ambient_aerosol_particles_in_air")
# att.put.nc(nc.idw,"pm2p5","cell_methods","NC_CHAR","time: mean (interval: 15 minutes)")
# att.put.nc(nc.idw,"pm2p5","missing_value","NC_FLOAT",-999.9)
# 
# # Global attributes
# att.put.nc(nc.idw,"NC_GLOBAL","title","NC_CHAR","PM2.5 interpolated surface (Inverse Square Distance)")
# att.put.nc(nc.idw,"NC_GLOBAL","Conventions","NC_CHAR","CF-1.7")
# att.put.nc(nc.idw,"NC_GLOBAL","Institution","NC_CHAR","NIWA (National Institute of Water and Atmospheric Research, Auckland, New Zealand)")
# att.put.nc(nc.idw,"NC_GLOBAL","project_id","NC_CHAR","CONA - 2018")
# att.put.nc(nc.idw,"NC_GLOBAL","history","NC_CHAR",paste0(format(max(all.data.tavg$date),format = "%Y%m%d"),
#                                                          " Data generated and formatted"))
# att.put.nc(nc.idw,"NC_GLOBAL","comment","NC_CHAR","Data for visualisation only")
# 
# # Load data
# var.put.nc(nc.idw,"latitude",lat_dim)
# var.put.nc(nc.idw,"longitude",lon_dim)
# var.put.nc(nc.idw,"time",as.numeric(tim_dim))
# rast_data <- getValues(raster_cat_idw_LL)[,(1:length(tim_dim))]
# dim(rast_data) <- c(length(lon_dim),
#                     length(lat_dim),
#                     length(tim_dim))
# var.put.nc(nc.idw,"pm2p5",rast_data)
# 
# # Close the file and save
# close.nc(nc.idw)

# IDW2
lat_dim <- unique(coordinates(raster_cat_idw2_LL)[,2])
lon_dim <- unique(coordinates(raster_cat_idw2_LL)[,1])
tim_dim <- all_dates[valid_dates]
nc.idw2 <- create.nc("odin_idw2.nc")
# Dimensions specifications
dim.def.nc(nc.idw2, "time", unlim=TRUE)
dim.def.nc(nc.idw2, "latitude",length(lat_dim))
dim.def.nc(nc.idw2, "longitude",length(lon_dim))
# Variable specifications
var.def.nc(nc.idw2,"time","NC_INT","time")
att.put.nc(nc.idw2,"time","units","NC_CHAR","seconds since 1970-01-01 00:00:0.0")
att.put.nc(nc.idw2,"time","long_name","NC_CHAR","time")

var.def.nc(nc.idw2,"latitude","NC_FLOAT","latitude")
att.put.nc(nc.idw2,"latitude","units","NC_CHAR","degrees_north")
att.put.nc(nc.idw2,"latitude","long_name","NC_CHAR","latitude")
att.put.nc(nc.idw2,"latitude","standard_name","NC_CHAR","latitude")

var.def.nc(nc.idw2,"longitude","NC_FLOAT","longitude")
att.put.nc(nc.idw2,"longitude","units","NC_CHAR","degrees_east")
att.put.nc(nc.idw2,"longitude","long_name","NC_CHAR","longitude")
att.put.nc(nc.idw2,"longitude","standard_name","NC_CHAR","longitude")

var.def.nc(nc.idw2,"pm2p5","NC_FLOAT",c("longitude","latitude","time"))
att.put.nc(nc.idw2,"pm2p5","units","NC_CHAR","ug m**-3")
att.put.nc(nc.idw2,"pm2p5","long_name","NC_CHAR","Mass concentration of PM2.5 ambient aerosol particles in air")
att.put.nc(nc.idw2,"pm2p5","standard_name","NC_CHAR","mass_concentration_of_pm2p5_ambient_aerosol_particles_in_air")
att.put.nc(nc.idw2,"pm2p5","cell_methods","NC_CHAR","time: mean (interval: 15 minutes)")
att.put.nc(nc.idw2,"pm2p5","missing_value","NC_FLOAT",-999.9)

# Global attributes
att.put.nc(nc.idw2,"NC_GLOBAL","title","NC_CHAR","PM2.5 interpolated surface (Inverse Square Distance)")
att.put.nc(nc.idw2,"NC_GLOBAL","Conventions","NC_CHAR","CF-1.7")
att.put.nc(nc.idw2,"NC_GLOBAL","Institution","NC_CHAR","NIWA (National Institute of Water and Atmospheric Research, Auckland, New Zealand)")
att.put.nc(nc.idw2,"NC_GLOBAL","project_id","NC_CHAR","CONA - 2018")
att.put.nc(nc.idw2,"NC_GLOBAL","history","NC_CHAR",paste0(format(max(all.data.tavg$date),format = "%Y%m%d"),
                                                          " Data generated and formatted"))
att.put.nc(nc.idw2,"NC_GLOBAL","comment","NC_CHAR","Data for visualisation only")

# Load data
var.put.nc(nc.idw2,"latitude",lat_dim)
var.put.nc(nc.idw2,"longitude",lon_dim)
var.put.nc(nc.idw2,"time",as.numeric(tim_dim))
rast_data2 <- getValues(raster_cat_idw2_LL)[,(1:length(tim_dim))]
dim(rast_data2) <- c(length(lon_dim),
                     length(lat_dim),
                     length(tim_dim))
var.put.nc(nc.idw2,"pm2p5",rast_data2)

# Close the file and save
close.nc(nc.idw2)

## Create MP4 video ####
# system(paste0("ffmpeg -f image2 -r 6 -pattern_type glob -i '",
#               data_path,
#               "../idw/",
#               "*.png' ",
#               data_path,
#               "../idw/",
#               format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
#               format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
#               ".mp4"))
system(paste0("ffmpeg -f image2 -r 6 -pattern_type glob -i '",
              "/home/gustavo/data/ODIN_SD/Gisborne/idw2/*.png'",
              "*.png' ",
              plot_path,
              "idw2/",
              format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
              format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
              ".mp4"))



## Upload to youtube ####
system(paste0("youtube-upload --title=\"Gisborne ",
              format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d %H:%M"),
              " to ",
              format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d %H:%M"),
              "\" --client-secrets=client_secrets.json ",
              plot_path,
              "idw2/",
              format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
              format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
              ".mp4 --playlist=\"Gisborne 2018 - ODIN\""))

# Compress TXT files ####
system(paste0("tar -zcvf ",
              data_path,
              'all.data',
              format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
              format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
              ".tgz ",
              data_path,
              'all.data',
              format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
              format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
              ".txt"))

system(paste0("tar -zcvf ",
              data_path,
              'all.data.tavg',
              format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
              format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
              ".tgz ",
              data_path,
              'all.data.tavg',
              format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
              format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
              ".txt"))

## Upload data ####

# RCurl::ftpUpload(paste0(data_path,"odin_idw.nc"),
#                  "ftp://ftp.niwa.co.nz/incoming/GustavoOlivares/odin_alexandra/odin_idw.nc")
RCurl::ftpUpload(paste0(data_path,"odin_idw2.nc"),
                 "ftp://ftp.niwa.co.nz/incoming/GustavoOlivares/odin_alexandra/odin_idw2_gisborne.nc")

RCurl::ftpUpload(paste0(data_path,
                        'all.data',
                        format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
                        format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
                        ".tgz"),
                 paste0("ftp://ftp.niwa.co.nz/incoming/GustavoOlivares/odin_alexandra/",
                        'all.data_gisborne',
                        format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
                        format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
                        ".tgz"))
RCurl::ftpUpload(paste0(data_path,
                        'all.data.tavg',
                        format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
                        format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
                        ".tgz"),
                 paste0("ftp://ftp.niwa.co.nz/incoming/GustavoOlivares/odin_alexandra/",
                        'all.data.tavg_gisborne',
                        format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
                        format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
                        ".tgz"))

## Remove files ####
# system(paste0("rm -rf ",
#               data_path,
#               "../idw/*"))
# system(paste0("rm -rf ",
#               data_path,
#               "../idw2/*"))
# system(paste0('rm -f ',
#               data_path,
#               'all.data',
#               format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
#               format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
#               ".txt"))
# system(paste0('rm -f ',
#               data_path,
#               'all.data.tavg',
#               format(min(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),"_",
#               format(max(all.data.tavg$date) + 12*3600,format = "%Y%m%d"),
#               ".txt"))
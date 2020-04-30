library(rgdal) # useful spatial library


#setwd("C:/Users/mklopez/Documents/Scratch")

form <- read.csv('coords_2020.csv')

str(form)
str(gps.tabl)

#####used for ODK GPS table, when bringing in data from the field and needing to convert it
#####the following attributes in quotes are how the data comes in

gps.tbl <- subset(form, select=c("group_wr94l28.Site_Code", "group_wr94l28.Tree_Number",
                                 "Tree_Location.Latitude", "Tree_Location.Longitude",     
                                 "Tree_Location.Altitude"))

#####I use this for data copied from the Access Database
gps.tabl <- subset(form, select=c("GPSID", "Easting",
                                 "Northing", "Latitude",     
                                 "Longitude"))

####this is some generic stuff that I change from time to time, whenever I need to convert something###
gps.tabl <- subset(form, select=c("Easting", "Northing", "series"))

# define this as a spatial object (SpatialPointsDataFrame)
coordinates(gps.tabl) <- c("Longitude", "Latitude")
coordinates(gps.tabl) <- c("Easting", "Northing")
coordinates(gps.tabl) <- c("X_UTM", "Y_UTM")
# Define the projection system via CRS (this is true for all tablet data

proj4string(gps.tabl) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #for lat long in your csv

proj4string(gps.tabl) <- CRS('+proj=utm +zone=11 +north +units=m +ellps=GRS80')  #for UTM in your csv, make sure correct datum is placed

# Transform the WGS84 to UTM (Jemez is zone 13, so change that)

gps.utm <- spTransform(gps.tabl, CRS("+proj=utm +zone=13 +north +units=m +ellps=GRS80"))

# Transform the UTM to WGS84

gps.latlong <- spTransform(gps.tabl, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# merge these new data frames and include both lat-long and utm in your database table (like a badass) 

gps.merge.tbl <- merge(gps.tabl, gps.latlong, by=c('series'))

write.csv(gps.merge.tbl, 'lomacoords.csv')

??rgdal

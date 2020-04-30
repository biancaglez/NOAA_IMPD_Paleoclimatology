## Query script preparing Jemez data for the The International Multiproxy Paleofire Database (IMPD) 
## at the World Data Center for Paleoclimatology using Jemez data.
## Author: Bianca Gonzalez 
rm(list=ls(all=TRUE))
setwd("C:/Users/bgonzalez/Desktop/NOAA/NOAA_IMPD_Paleoclimatology")

## write note with file path and notes on the script on things that didn't work/worked 
## do this for NOAA and for SFNF and for the storymap -- relfect etc

library(dplR)
library(burnr)
library(vctrs)
library(stringr)
library(dplyr)
library(tidyverse)
library(rgdal)
library(raster)
library(rgdal)
library(sp)
library(tigris)



## read in all fhx files in folder 
fhx_all <- list.files(path = "C:/Users/bgonzalez/Desktop/access_databases/Jemez.database.merge/jemez.final.fhx.files",  
                      pattern = "*fhx", full.names = TRUE) 

# Identify all csv files in folder
fhx<- 
  read_fhx("C:/Users/bgonzalez/Desktop/access_databases/Jemez.database.merge/jemez.final.fhx.files/mpb.FHX")

#### NOTES BELOW from Ellis for metadata sheet 


#### DONE: FHX_filename: IE POT.fhx ####

# initally run on the sites that were already in there
# cleanest data are the ones that are already in there - fhx folder? 
# all.fhx.jemez -- start with metadata that is in access database already -- 
# prior to my era jemez_final_fhx -- these are in fire history jemez final fhx files 
# the february

# queries also in access database 

# read in all fhx files in folder
list_files <- lapply(fhx_all, read_fhx)

filenames<- list()
for(i in seq_along(list_files)){
  
  # sort of hack getting rid of 0-9 nums as well as the - symbol 
  no_num<-sub("([0-9]+).*$", "", list_files[[i]]$series[1])
  filenames[[i]]<-sub("\\-.*", "", no_num)
  
}

#### Needs Ellis attention: contributors ####
# concatenate collectors? 
# for Ellis to deal with 



#### DONE 4/21/2020: Site_code is the FHX filename ####
# IE POT 
fhx_filenames <- unlist(filenames)

# this we have to get from the access database -- might be able to write a query script for each
# fhx file -- ask Ellis! 
# site code is abbreviation -- subsetting the data later 
# write a little script ^
# Access database: sites table -- site code 
#compare fhx filename (chopped off) with the sites table -- abbreivation 

#### DONE: First_year [AD] ####
# pull from access database SQL min and max from database of year 
# min max of unique sample ID per site ---  
# years for all of the trees of this site 
# min/max of years per site 

# DONE: Last_year [AD]
#### DONE: Species_code [of sampled species]: PIPO, PSME #### 
# give me all species for all trees in site X Y Z  - with comma in between 
# right now I only have single instance of that species 

#### DONE: Collection_date: ####
#Year only 2017 - one year or both? (comma in between? possible formatting)follow the template she provided #

#### DONE: Number_samples ####
# -- use & in between years -- collection trip table - count of # of trees

#### DONE: collectors is just concatenating the collectors ####
#### UTM_easting [include datum and zone]: 466386 (NAD 83, zone 13N) ####
#### UTM_northing [include datum and zone] ####
# common name 
#### DONE 4/21/2020: crossdaters -  look in database for this -- Ellis can help clean this up ####
# prefer the table query
# look at existing queries
# concatenate crossdaters? 

#### DONE: Lat/Long ####
# most lat/long should be in the database
# need to average Lat long of each site
# Latitude  [decimal degrees  (WGS 84)] IE: 36.4277 
# Longitude   [decimal degrees  (WGS 84)] IE: -105.375

pts<- (c(3975244.477,	358312.1521))
+zone=13
sputm <- SpatialPoints(pts, proj4string=CRS("+proj=utm +zone=13 +datum=WGS84"))  
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))

#### DONE 04/21/2020: Mean_elevation [meters] IE: 2690 ####

# average the coordinates to a single pt in access and you get a center coordinate 
# raster extract() script after that
                       

###read in FULL jemez 30 m DEM in UTM 13 proj - exported dem from ARc as geotiff 
dem.30m<-raster("C:/Users/bgonzalez/Desktop/NOAA/NOAA_IMPD_Paleoclimatology/data/Jemez30m_DEM_UTM1/Jemez30m_DEM_UTM1.tif")

# let's use and keep the original projection
dem.30m@crs

# leave to given projection -- reproject the pts to the dem projection 
# +proj=utm +zone=13 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs  

#read in csv with northing and easting of site coordinates (exported from database)
jemez.coords<-read.csv(paste0(getwd(), "/","data", "/", "NOAA-eastnorth.csv"))

jemez.pts <-jemez.coords[,1:2] #grab only east/north

### make coords a spatial data points data frame with coordinates() - identify col names that hold spatial info
coordinates(jemez.pts)<- c("AvgOfEasting","AvgOfNorthing")

##set CRSto raster CRS
jemez.pts<- SpatialPoints(jemez.pts, proj4string=CRS("+proj=utm +zone=13 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))  

## make sure they overlap and are in same CRS
plot(dem.30m)
plot(jemez.pts,  col = "purple", add = TRUE)

# they do! 

#extract elevation from DEM for site avg pts 
elev.fs.trees<- raster::extract(dem.30m, jemez.pts@coords,method='bilinear')

#add back to df with site ids
cbind(elev.fs.trees, jemez.coords[3:4])

# -- later to do: run r script to get elevation then fill out elevation in gpsid and then query the elevation
                       
#### DONE:  Country: United States####

country<- rep("USA", length(list_files))
write.csv(country, "filepath")

#### State: New Mexico ####
# Bring in lat/ long points from each site and confirm if they are within the NM boundary 
# autopopulated

##### DONE: Area_sampled [hectares]for one site #### 

#VTN and FWH sites -- no lat long or easting northing provided  
# calculate in R 
# coordinates for all trees in a site - convex hole
# populate back in database in R 

# they have come in as LAT LONG -- the workflow was changed? ODK import database ? 
# they came in as lat long -- in ODK form -- they come in as LAT LONG -- select LAT LONG 
# they were always WGS84? -- yes ODK collects as WGS84? 
# they can be established as WGS84 --  
# use the UTMS instead? take all coords in there and reproject them -- 

convex.coords<-read.csv(paste0(getwd(), "/","data", "/", "NOAA_convex_utms.csv"))

# going to need to group by abbreviation and then pass to lapply fn 
# try to see the reality space and also check out the units of this projection again 
bol<-convex.coords %>% filter(Abbreviation=='BOL') %>% dplyr::select("Easting", "Northing")

# make spatial object first 
coordinates(bol) <- c("Easting", "Northing")

# set CRS using Mannie's for UTMS
proj4string(bol) <- CRS('+proj=utm +zone=13 +north +units=m +ellps=GRS80')

# does CHULL take spatial points?
mm<-chull(bol@coords) # convex hull eq. gives row #s (integar vector) of lat long that are farthest pts 
plot(bol@coords)

convex_coords_bol <- bol@coords[c(mm, mm[1]), ]  # select rows of mm and the last pt to join together
lines(convex_coords_bol, col="red")

# make convex_coords into spatial polygon

#list within list of polygon and polygons
# Polygon/Polygons - creates obj of spatial class Spatial Polygon
# POlygons ID	arg is character vector of length one with identifier

sp_poly<- SpatialPolygons(list(Polygons(list(Polygon(convex_coords_bol)), ID=1))) # sp_poly@proj4string # no CRS yet

#set crs
sp_poly@proj4string<- CRS("+proj=utm +zone=13 +datum=WGS84")

# lets make sure BOL is in the Jemez first - visual check 
plot(dem.30m)
plot(sp_poly, add=TRUE, col="red")
lines(convex_coords_bol, col="red") # the tiniest little blimp in the lower half center

#units are of current projection - squared meters !
# area fn returns spatial object in squared meters https://www.rdocumentation.org/packages/raster/versions/3.0-12/topics/area
sp_poly$area_hectares <- area(sp_poly) / 1000 # make hectares

#### DONE: -- Area in hectares for all sites ####

convex.coords<-read.csv(paste0(getwd(), "/","data", "/", "NOAA_convex_utms.csv"))
#store convex polygons in list 
sp_list<-list()
#store areas in vector
vec<- vector(length=length((unique(convex.coords$Abbreviation))))
options(scipen=999) # set global options - Positive values bias towards fixed -- get rid of scientific notation

for(i in seq_along(unique(convex.coords$Abbreviation))){
  print(i)
      ###### newewst edits
      # make spatial object first 
      temp<- convex.coords %>% filter(Abbreviation==unique(convex.coords$Abbreviation)[i]) %>% dplyr::select("Easting","Northing") # need to change to i
      coordinates(temp) <- c("Easting", "Northing")
      
      # set CRS using Mannie's for UTMS
      proj4string(temp) <- CRS('+proj=utm +zone=13 +north +units=m +ellps=GRS80')
      
      # CHULL needs the matrix of coords
      mm<-chull(temp@coords) # convex hull eq. gives row #s (integar vector) of lat long that are farthest pts 
      convex_coords_temp <- temp@coords[c(mm, mm[1]), ]  # select rows of mm and the last pt to join together
      
      # make convex_coords into spatial polygons, set CRS, Calc area
      sp_list[[i]]<- SpatialPolygons(list(Polygons(list(Polygon(convex_coords_temp)), ID=1)))
      sp_list[[i]]@proj4string<- CRS("+proj=utm +zone=13 +datum=WGS84")
      
      # area fn returns spatial object in squared meters https://www.rdocumentation.org/packages/raster/versions/3.0-12/topics/area
      sp_list[[i]]$area_hec <- area(sp_list[[i]]) / 10000 
      
      vec[i]<- sp_list[[i]]$area_hec    # now let's take all area_hec out of the list and add to a vector 
      
      # this is a normal number now yay! 
      }

# cbind vector with the unique instances of sites
area_sites<- as.data.frame(cbind(vec, unique(convex.coords$Abbreviation)))

#visually check code ^ 
plot(dem.30m)
plot(sp_list[[1]], add=TRUE, col= "red") # it is on the DEM! wahoo! 

# sources: 
# https://gis.stackexchange.com/questions/177106/area-units-for-output-polygon-in-wgs-84-utm-zone-32n-using-qgis
# https://stackoverflow.com/questions/25606512/create-convex-hull-polygon-from-points-and-save-as-shapefile

#### NEEDS REVIEW: County #### 
# center coordinates for the site- and getting the township -- 
# from a polygon/raster layer -- for the NA fire scar  
# right datum etc 

# bring in county layer 
counties <- shapefile(paste0(getwd(), "/", "data", "/", "tl_2019_us_county", "/", "tl_2019_us_county.shp"))
#counties@proj4string<- CRS("+proj=utm +zone=13 +datum=WGS84") # change to local CRS for this NM purpose

#NM counties
t_cty <- tigris::counties(state="35")

#read in csv with lat and long of site coordinates (exported from DB)
jemez.latlong<-read.csv(paste0(getwd(), "/","data", "/", "NOAA-latlong.csv"))
jemez.gcs <-jemez.latlong[,1:2] #grab only lat/long - geographic coordinate system (gcs)

# originally collected in UTMS and then converted to LAT/LONG -- ??? basic lat long projection for NM and test a few
# plot on counties or jemez DEM 
# need to assign the coords their home projection -- but where do I grab this info? -- coordinate system projection col?
# bring in DEM to do a check -- original projection assigned may not be true to what these points were originally projected to? 

### make coords a spatial data points data frame with coordinates() - identify col names that hold spatial info
coordinates(jemez.gcs)<- c("AvgOfLatitude","AvgOfLongitude")

# double check projection 
jemez.gcs<- SpatialPoints(jemez.gcs, proj4string=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))  ##set CRS to relevent UTM

# SpatialPointsDataFrame(jemez.gcs@coords, data = data.frame())
# use extract function here! 
#extract elevation from DEM for site avg pts 

# do these shapes intersect at all? T/F 
rgeos::gIntersects(tigris::states(), jemez.gcs)

#county first
sp::over(t_cty, jemez.gcs)

raster::intersect(t_cty, jemez.gcs)

meep<-SpatialPointsDataFrame(jemez.gcs@coords, proj4string=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0 "), df)
rgeos::gIntersects(meep, townships)


#### NEEDS REVIEW: Township: T25N ####
# center coordinates for the site- and getting the township -- 
# from a polygon/raster layer -- for the NA fire scar  

townships<- shapefile((paste0(getwd(),"/","data", "/", "PLSSTownship", "/", "PLSSTownship.shp")))


#read in csv with lat and long of site coordinates (exported from DB)
jemez.latlong<-read.csv(paste0(getwd(), "/","data", "/", "NOAA-latlong.csv"))
jemez.gcs <-jemez.latlong[,1:2] #grab only lat/long - geographic coordinate system (gcs)

coordinates(jemez.gcs) <- c("AvgOfLatitude", "AvgOfLongitude")

# set CRS to townships CRS
proj4string(jemez.gcs) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')

# row 147 is failing so let's take that row out 
jemez.gcs@coords <- jemez.gcs@coords[-c(147),] 

proj4string(jemez.gcs) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')

# all three over, intersect and gintersect functions dont work... 
sf::st_intersection(jemez.gcs@coords, townships)
rgeos::gIntersects(townships, jemez.gcs)
# townships first 
sp::over()
#### NEEDS REVIEW: Range: R15E #### 
# center coordinates for the site- and getting the township -- 
# from a polygon/raster layer -- for the NA fire scar 


#### NEEDS REVIEW: Lowest_elev [meters] 2658 meters #### 
# Highest_elev [meters]

# bring in pts grouped by site 
# find min elev 
# find max elev 



#### NEEDS REVIEW: EDIT NOAA Species output and concatenate #### 

#### NEEDS REVIEW: national/federal land ownership layer (for forest and park) ####
## can find the national layer here: https://www.usgs.gov/core-science-systems/science-analytics-and-synthesis/gap/science/pad-us-data-download?qt-science_center_objects=0#qt-science_center_objects
## takes too long to download in the park

NM_ownership <- shapefile((paste0(getwd(),"/","data", "/", "PADUS2_0NM_Shapefile", "/", "PADUS2_0Proclamation_NM.shp")))

NM_ownership@proj4string # this is crazy and doesnt plot on top of dem.30

proj4string(NM_ownership) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')

# intesect lat long points with this puppy 

#read in csv with lat and long of site coordinates (exported from DB)
jemez.latlong<-read.csv(paste0(getwd(), "/","data", "/", "NOAA-latlong.csv"))
jemez.gcs <-jemez.latlong[,1:2] #grab only lat/long - geographic coordinate system (gcs)

coordinates(jemez.gcs) <- c("AvgOfLatitude", "AvgOfLongitude")

# set to CRS from ownership layer 
proj4string(jemez.gcs) <- CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m
+no_defs +ellps=GRS80 +towgs84=0,0,0 ')

# now lets overlap them and see what we get 

sf::st_intersection(NM_ownership, jemez.gcs)

rgeos::gIntersects(NM_ownership, jemez.gcs)
# townships first 
sp::over()
raster::intersect()

# going to see what's going on here --- both are not plotting on the DEM
plot(jemez.gcs, add = TRUE)



## 

#### list of additional filenames
# for Ellis 
#### we can sample ranger district, national forest, park monument, township - 30 - 35 grab layers (first NM then nationwide) -- topo quads (line 40) - 
# tree level min/max elevation 


#query for all of the sites in the sites table -- 
# subset the sites and then add that to the query (feature of the query)
# start with all of the sites -- 
# monument canyon -- ignore 
# dont get hung up on legacy sites 
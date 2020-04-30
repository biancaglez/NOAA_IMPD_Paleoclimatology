# script to process Jemez data for NOAA

## Query script preparing Jemez data for the The International Multiproxy Paleofire Database (IMPD)
## at the World Data Center for Paleoclimatology using Jemez data.
## Author: Bianca Gonzalez
import glob
import pandas as pd

#read fhx files in list
path =r'C:\Users\bgonzalez\Desktop\access_databases\Jemez.database.merge\jemez.final.fhx.files'
all_files = glob.glob(path + "/*.fhx")



# FHX_filename. IE POT.fhx

# Site_code IE: POT

# Latitude  [decimal degrees  (WGS 84)] IE: 36.4277

# Longitude   [decimal degrees  (WGS 84)] IE: -105.375

# Mean_elevation [meters] IE: 2690

# Country: United States

# State: New Mexico

# Area_sampled [hectares]

# First_year [AD]

# Last_year [AD]

# Species_code [of sampled species]: PIPO,PSME

# Collection_date: Year only 2017

# Number_samples

# County

# Township: T25N

# Range: R15E

# UTM_easting [include datum and zone]: 466386 (NAD 83, zone 13N)

# UTM_northing [include datum and zone]

# Topographic_map
# Lowest_elev [meters] 2658 meters

# Highest_elev [meters]
# Slope
# Aspect
# Substrate_type:  pumice

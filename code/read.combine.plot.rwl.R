###1-28-2020
###reading rwl files in directory, combine, then plot sample depth and ring width measurements

# often a useful chunk of code to have at the beginning of each scrip is the followinw hich removes all named objects
# from the current R session
rm(list=ls(all=TRUE))
setwd("C:/ellis/research/Taos/BLM_RGDN_monument/climate_reconstruction/data/measurements/olla/all")
library(dplR)

# reads in names of all files with .rwl extension in the folder
#  OR.PSME.rwl.file.names <- list.files(pattern=".rwl$")

# shorter code to read in names of all files with .rwl extension and make a list
file.names<-Sys.glob("*.rwl")

#for loop that reads in the data for all file names in .rwl file name list created above
for(i in file.names) 
{ x <- read.rwl(i, ) #header=TRUE  #only add this code after "i," if there is a header - cuts out first 3? lines
# removes quotes from file names for column names
assign(print(i, quote=FALSE), x)}

# combines first two .rwl files by "get"ing and "paste" ing the 1st [1] and 2nd [2] file names from OR.PSME.rwl.file.names
aa<-combine.rwl(get(paste(file.names[1])), get(paste(file.names[2])))

#for loop that combines remaining files (3rd on) from file names with first two (aa)
for (i in 3:length(file.names)){
  aa<-combine.rwl(aa, get(paste(file.names[i])))
}

summary(all.rwl)


#rename
all.rwl<-aa

str(all.rwl)

#setwd("C:/ellis/research/Ceanothus2014/data/Measurements")

#writes rwl of all individual files joined, long.names allows 8 character sample names
#write.rwl(cea.all, 'cea.RWL', 
#          long.names = TRUE, 
#          prec=0.01)

#write csv file with all measurements
#write.csv(cea.all, 'cea.all.csv')

length(all.rwl) #29 series

#line plot of inner & outer dates of all samples
#set up tiff
tiff('./figs/age.plot.tiff', 
     height = 6, width=8, compression = 'lzw', type = "cairo", antialias = "gray", units="in",res=301)
jpeg('./figs/age.plot.jpeg', 
     height = 6, width=8, type = "cairo", antialias = "gray", units="in",res=301)

plot(all.rwl)
dev.off()


#plot all raw ring widths
#set up tiff
tiff('./figs/all.rwl.tiff', 
     height = 6, width=8, compression = 'lzw', type = "cairo", antialias = "gray", units="in",res=301)
jpeg('./figs/all.rwl.jpeg', 
     height = 6, width=8, type = "cairo", antialias = "gray", units="in",res=301)

spag.plot(all.rwl, zfac = 0.5)

dev.off()

#detrend with spline
rwi<-detrend(all.rwl, method = "Spline")

###plot ring width indices
spag.plot(rwi, zfac = 0.7)

#make default crn
crn<-chron(rwi, "oll")

plot(crn)



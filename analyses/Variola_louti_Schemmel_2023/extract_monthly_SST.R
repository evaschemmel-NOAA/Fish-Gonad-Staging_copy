#Code below is used get mean monthly water temp for Guam and create a location map
#Want month mean water temp from 2010 through 2017 - dont need summarized to year - just need monthly averages

#https://oceanwatch.pifsc.noaa.gov/erddap/search/index.html?page=1&itemsPerPage=1000&searchFor=CRW_sst_v1_0_monthly+
#mapping and erddap packages

library(ncdf4)
library(httr)
library(lubridate)
library(plyr)
library(rerddap)
library(sf)
library(maps)
library(mapdata)
library(maptools)
library(parsedate)
library(readr)
library(rgdal)
library(date)
library(plotdap)
library(raster)
library(sp)
library(rerddapXtracto)

#########


#Extract the regional boundaries from the downloaded shapefile

#shapefile but only need region 1 = Guam"
setwd("~/Documents/github/satellite_course/marianaregionshapefile")
marianashape=st_read("Mariana_Boundary.shp") # shapefile for all of Marianas that Eva created
mariana=st_coordinates(marianashape)


#Let's set the bounds of our region: 
xlim = c(144,147) 
ylim = c(12,21)



ERDDAP_Node="https://oceanwatch.pifsc.noaa.gov/erddap/"
dataInfo <- rerddap::info('CRW_sst_v1_0_monthly', url=ERDDAP_Node)
#Let's look at the variables in our dataset:
dataInfo$variable$variable_name

tcoord <- c("2010-01-01", "2017-12-31")

temp.rangeall=list()
I=which(mariana[,4]==1)  #select the feature # select only Guam region
poly=mariana[I,1:2]
#write.csv(mariana,'mariana.csv') 

poly=data.frame(poly)
names(poly) <- c("lon","lat")

xcoord <- poly$lon 
ycoord <- poly$lat

parameter=dataInfo$variable$variable_name[1]
#extract data
sst <- rxtractogon (dataInfo, parameter=parameter, xcoord=xcoord, ycoord=ycoord, tcoord=tcoord)

mean_sst=as.matrix(apply(sst$analysed_sst,c(1,2),mean,na.rm=TRUE))

sst <- rxtractogon (dataInfo, parameter=parameter, xcoord=xcoord, ycoord=ycoord, tcoord=tcoord)
str(sst)

```{r create averages using a for loop}

dims <- dim(sst$analysed_sst)
clim_mon <- array(dim=c(dims[1],dims[2],12))

dimnames(clim_mon) <- list(Longitude=sst$longitude, 
                           Latitude=sst$latitude, 
                           Month=month.abb[seq(1,12)])
str(clim_mon)#check to see if the dimmensions are correct
for (im in 1:12) {
  mtime <- which(month(sst$time) == im) 
  month1 <- sst$analysed_sst[,,mtime]
  avg <- apply(month1,c(1,2),function(x) mean(x,na.rm=TRUE))
  clim_mon[,,im] <- avg
}

```

```{r reformat data for ggplot}

clim_month_long <- reshape2::melt(clim_mon,value.name="sst")
```

#summarize to get average sst per month
summary<-clim_month_long  %>%
  dplyr::group_by(Month)%>%
  dplyr::summarise(mean_sst=mean(sst, na.rm=TRUE))


```{r Map the data }

# Use the facet feature of ggplot to automatically create 12 maps.
library(ggplot2)
library(ggmap)
library(mapdata)
library(MAP)
coast <- map_data("worldHires")
coast <- map_data("world", ylim = ycoord, xlim = xcoord)
coast <- map_data("world2", ylim = ycoord, xlim = xcoord)
map.world <- map_data(map="world")
#png(file="Monthly_Averages.png")
ggplot(
  data = clim_month_long, 
  aes(x = Longitude, y = Latitude, fill = sst)) +
  geom_tile(na.rm=T) +
  geom_polygon(data = coast, aes(x=long, y = lat, group = group), fill = "grey80") + #where do we get the coast layer?
  theme_bw(base_size = 12) + ylab("Latitude") + xlab("Longitude") +
  coord_fixed(1.3, xlim = c(144,146), ylim = c(12.5,14)) + #set the spatial limits of the plot
  scale_y_continuous(expand=c(0,0)) + 
  scale_x_continuous(expand=c(0,0)) +
  scale_fill_gradientn(colours = rev(rainbow(12)), na.value = NA) + 
  facet_wrap(~ Month)


#dev.off()


#####basic marianas map
install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial") #scale bar and north arrow


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

## [1] "sf"  
## [1] "data.frame"
st_crs(coast)



library(ggplot2)
library(ggmap)
library(mapdata)
library(MAP)
coast <- map_data("worldHires")
coast <- map_data("world", ylim = ycoord, xlim = xcoord)
coast <- map_data("world2", ylim = ycoord, xlim = xcoord)
map.world <- map_data(map="world")

marianas<-ggplot() +
  geom_polygon(data = coast, aes(x=long, y = lat, group = group), fill = "black") + #where do we get the coast layer?
  geom_sf() +
  coord_sf(xlim = c(144, 147), ylim = c(13,21), expand = FALSE)+
   xlab("Longitude") + ylab("Latitude") +
  theme_bw( ) +
  annotate(geom = "text", x = 145.4, y = 13.5, label = "Guam", color = "grey22", size = 4)+
  annotate(geom = "text", x = 146.3, y = 18.1, label = "Pagan", color = "grey22", size = 4)+
  annotate(geom = "text", x = 146.1, y = 19.7, label = "Asuncion", color = "grey22", size = 4)+
  annotate(geom = "text", x = 145.8, y = 20.1, label = "Maug", color = "grey22", size = 4)+
  annotate(geom = "text", x = 145.5, y = 20.5, label = "Uracus", color = "grey22", size = 4)+
  theme(text=element_text(  family="Times ", size=15,colour="black"), axis.text.x= element_text(colour="white"))
  
#didnt work to reset x axis labels
marianas + scale_x_discrete(breaks=c(144.0,144.5,145.0,1345.5,146.0,146.5, 147.0),
                       labels=c("144"," ","145"," ","146","", "147"))  
  
marianas
ggsave("marianas.png", width = 3, height = 6, dpi = 300)


###from Mike Kinney
library(ggplot2)
library(ggspatial)
worldmap <- map_data ("world2")

ggplot()+
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill="darkgrey", colour = "black") +
  labs(x = "Longitude",y="Latitude")+
  coord_fixed(xlim = c(144, 147) ,ylim= c(12,22))+
  annotation_scale(location = "bl", width_hint = 0.5)+
  annotation_north_arrow(location = "bl",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)


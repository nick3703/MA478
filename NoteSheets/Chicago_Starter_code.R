library(tidyverse)
library(sf)


x = as(sf::read_sf("Graded Events/Projects/Chicago_data/chipocsouth.shp"), "sf")

x2<-subset(x, select = !duplicated(names(x)))



test<- read.csv("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/crime.csv")

test <- test %>% dplyr::select(-X)


tot_burg <- apply(test,1,sum)

x2$tot_burg <- tot_burg

ggplot() +
  geom_sf(data = x2, aes(fill= tot_burg))+
  scale_fill_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous"))


time_burg <- apply(test,2,sum)


test<-readMM("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/neighborhood.mtx")

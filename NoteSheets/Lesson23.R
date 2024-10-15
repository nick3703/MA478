library(tidyverse)
library(faraway)
library(INLA)
library("spdep")
library("lme4")
library("spatstat")
library("sp")
library("maptools") #This is depricated
library("lattice")
library(sf)

data(penicillin)
summary(penicillin)


#Prior for Precision
prec.prior <- list(prec = list(param = c(0.001, 0.001)))

penicillin$treat <- relevel(penicillin$treat, "D")


inla.pen <- inla(yield ~ 1 + treat + f(blend, model = "iid",
                                       hyper = prec.prior),
                 data = penicillin, control.predictor = list(compute = TRUE))
summary(inla.pen)


test<-inla.pen$marginals.fixed[[3]]


ggplot(data.frame(inla.smarginal(test)), aes(x, y)) +
  geom_line() +
  theme_bw()


Z <- as(model.matrix(~ 0 + blend, data = penicillin), "Matrix")


data("sleepstudy")

sleepstudy$Reaction <- sleepstudy$Reaction / 1000


inla.sleep <- inla(Reaction ~ 1 + Days + f(Subject, model = "iid"),
                   data = sleepstudy)


#sleepstudy %>% ggplot(aes(x=Days,y=Reaction))+geom_point()+
#  stat_smooth(method="lm",se=F)


#sleepstudy %>% ggplot(aes(x=Days,y=Reaction))+geom_point()+
#  stat_smooth(method="lm",se=F)+facet_wrap(~Subject)


inla.sleep.w <- inla(Reaction ~ 1 + f(Subject, Days, model = "iid"),
                     data = sleepstudy, control.predictor = list(compute = TRUE))
summary(inla.sleep.w)

#Subject is the group, Days is allowed to differ between group



#Predicted Values

post.means <- inla.sleep.w$summary.linear.predictor[,1]

sleepstudy.fitted <- sleepstudy
sleepstudy.fitted$postmeans<-post.means

#Visualizing

sleepstudy.fitted %>% ggplot(aes(x=Days,y=Reaction,col="Fitted Line"))+geom_line(aes(y=post.means,col="Posterior Mean"),lwd=1)+
  stat_smooth(method="lm",se=F)+geom_point(aes(y=Reaction,col="Values"))+facet_wrap(~Subject)




## -----------------------------------------------------------------------------
nyc.stops <- read.table(file = "frisk_with_noise.dat", skip = 6, 
                        header = TRUE)

# Add labels to factors
nyc.stops$eth <- as.factor(nyc.stops$eth)
levels(nyc.stops$eth) <- c("black", "hispanic", "white")
nyc.stops$eth <- relevel(nyc.stops$eth, "white")

nyc.stops$crime <- as.factor(nyc.stops$crime)
levels(nyc.stops$crime) <- c("violent", "weapons", "property", "drug")


## -----------------------------------------------------------------------------
nyc.stops.agg <- aggregate(cbind(stops, past.arrests, pop) ~ precinct + eth,
                           data = nyc.stops, sum)

# Population is summed 4 times
nyc.stops.agg$pop <- nyc.stops.agg$pop / 4


## -----------------------------------------------------------------------------
nyc.inla <- inla(stops ~ eth + f(precinct, model = "iid"),
                 data = nyc.stops.agg, offset = log((15 / 12) * past.arrests),
                 family = "poisson")

summary(nyc.inla)


## -----------------------------------------------------------------------------
# Ethnicity precinct index
nyc.stops.agg$ID <- 1:nrow(nyc.stops.agg)
nyc.inla2 <- inla(stops ~ eth + f(precinct, model = "iid") + 
                    f(ID, model = "iid"), 
                  data = nyc.stops.agg, offset = log((15/12) * past.arrests),
                  family = "poisson")

summary(nyc.inla2)


## ----label = "eff-precinct", fig = TRUE, echo = FALSE, fig.cap = "Posterior means and 95% credible intervals of the random effects of the 75 NYC precincts."----
tab <- nyc.inla2$summary.random$precinct
names(tab)[c(4, 6)] <- c("IClower", "ICupper")

ggplot(data = tab, aes(x = mean, y = ID)) +
  geom_point() +  ylim(1, 75) + xlab(expression(u[p])) + ylab ("Precinct") +
  geom_errorbarh(aes(xmin = IClower, xmax = ICupper) ) +
  geom_vline (xintercept = 0, linetype = 2)


#With DIC

nyc.inla <- inla(stops ~ eth + f(precinct, model = "iid"),
                 data = nyc.stops.agg, offset = log((15 / 12) * past.arrests),
                 family = "poisson",control.compute = list(dic=TRUE,config = TRUE))


nyc.inla2 <- inla(stops ~ eth + f(precinct, model = "iid") + 
                    f(ID, model = "iid"), 
                  data = nyc.stops.agg, offset = log((15/12) * past.arrests),
                  family = "poisson",control.compute=list(dic=TRUE,config = TRUE))


nyc.inla$dic$dic

nyc.inla2$dic$dic

ns=100

#Posterior Samples
xs = inla.posterior.sample(ns,nyc.inla)

xs2 = inla.posterior.sample(ns,nyc.inla2)




## Spatial


library("spatstat")
library("sp")
library("maptools")
library("lattice")
data(bei)

# Create SpatialPoints object
bei.pts <- as(bei, "SpatialPoints")



#Create grid
bei.poly <- as(as.im(bei$window, dimyx=c(5, 10)), "SpatialGridDataFrame")
bei.poly <- as(bei.poly, "SpatialPolygons")




#Number of observations per cell
idx <- over(bei.pts, bei.poly)
tab.idx <- table(idx)
#Add number of trees
d <- data.frame(Ntrees = rep(0, length(bei.poly)))
row.names(d) <- paste0("g", 1:length(bei.poly))
d$Ntrees[as.integer(names(tab.idx))] <- tab.idx



bei.trees <- SpatialPolygonsDataFrame(bei.poly, d)

#Mapping
idx.mapping <- as.vector(t(matrix(1:50, nrow = 10, ncol = 5)))
bei.trees2 <- bei.trees[idx.mapping, ]


#Summary statistics of covariates
covs <- lapply(names(bei.extra), function(X) {
  layer <- bei.extra[[X]]
  res <- lapply(1:length(bei.trees2), function(Y) {
    summary(layer[as.owin(bei.trees2[Y, ])])})
  res <- as.data.frame(do.call(rbind, res))
  
  names(res) <- paste0(X, ".", c("min", "1Q", "2Q", "mean", "3Q", "max"))
  return(res)
})

covs <- do.call(cbind, covs)

#Add to SPDF
bei.trees2@data <- cbind(bei.trees2@data, covs)



bei.adj.q <- poly2nb(bei.trees2)
bei.adj.r <- poly2nb(bei.trees2, queen = FALSE)


spplot(bei.trees2, "Ntrees", col.regions = rev(gray.colors(16)),
       sp.layout = list(list("sp.points", bei.pts))
)



require(gridExtra)

p1 <- spplot(bei.trees2, "elev.mean", col.regions = rev(hcl.colors(16, "RedOr")))
p2 <- spplot(bei.trees2, "grad.mean", col.regions = rev(hcl.colors(16, "OrYel")))

grid.arrange(p1, p2, ncol = 2)






## ---- fig = TRUE, echo = FALSE, label = "adj", fig.cap = "Queen (left) versus rook (right) adjacency."----
par(mfrow = c(1,2), mar = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))

bei.xy <- coordinates(bei.trees2)

plot(bei.trees2)#, main = "Queen adjacency")
plot(bei.adj.q, bei.xy, add = TRUE)

plot(bei.trees2)#, main = "Rook adjacency")
plot(bei.adj.r, bei.xy, add = TRUE)


W.bin <- nb2listw(bei.adj.q, style = "B")
W.rs <- nb2listw(bei.adj.q, style = "W")
W.bin



library("INLA")

#Log-Poisson regression
m0 <- inla(Ntrees ~ elev.mean + grad.mean, family = "poisson",
           data = as.data.frame(bei.trees2) )

#Log-Poisson regression with random effects
bei.trees2$ID <- 1:length(bei.trees2)
m0.re <- inla(Ntrees ~ elev.mean + grad.mean + f(ID), family = "poisson",
              data = as.data.frame(bei.trees2) )



####



london.gen <- rgdal::readOGR("LDNSuicides.shp")

x = as(sf::read_sf("LDNSuicides.shp"), "Spatial")

library(spdep)
temp <- poly2nb(x)
nb2INLA("LDN.graph",temp)
H <- inla.read.graph(filename="LDN.graph")
image(inla.graph2matrix(H))
LDN.adj<- paste(getwd(),"/LDN.graph",sep="")

#inla.matrix2graph()


#--- The data ---#
y=c(75,145,99,168,152,173,152,169,130,117,124,119,134,90,
    98,89,128,145,130,69,246,166,95,135,98,97,202,75,100,
    100,153,194)
E=c(80.7,169.8,123.2,139.5,169.1,107.2,179.8,160.4,147.5,
    116.8,102.8,91.8,119.6,114.8,131.1,136.1,116.6,98.5,
    88.8,79.8,144.9,134.7,98.9,118.6,130.6,96.1,127.1,97.7,
    88.5,121.4,156.8,114)

#Social Depravation (x1)
#Social Fragmentation (x2)

x1=c(0.87,-0.96,-0.84,0.13,-1.19,0.35,-0.84,-0.18,-0.39,0.74,
     1.93,0.24,0.59,-1.15,-0.8,-0.59,-0.12,1.43,-0.04,-1.24,1,
     0.53,-0.75,1.36,-0.93,-1.24,1.68,-1.04,2.38,0.03,-0.2,0.14)
x2=c(-1.02,-0.33,-1.43,-0.1,-0.98,1.77,-0.73,-0.2,-0.96,-0.58,
     0.36,1.48,0.46,-0.93,-1.62,-0.96,-0.48,0.81,2.41,-0.4,
     0.71,-0.05,-0.33,-0.47,-0.92,0.06,0.22,-0.73,0.1,-0.59,
     0.7,2.28)

names<- sort(london.gen$NAME)
data <- data.frame(NAME=names, y=y, E=E, x1=x1, x2=x2)
Nareas <- length(data[,1])

#E is equivalent to saying offset=log(E)



boroughs=london.gen
data.boroughs=attr(boroughs, "data")
order <- match(data.boroughs$NAME,data$NAME)
data <- data[order,]
ID<-seq(1,32)
data <- cbind(ID,data)

print(data)


              
plot(x)

formula1 <- y ~ 1+f(ID,model = "iid",
                   hyper=list(prec=list(prior="loggamma",
                                        param=c(1,0.00001))))
mod1 <- inla(formula1,family="poisson",data=data,E=E)



x@data$county_effect <- mod1$summary.random$ID$mean

x<-as(x,"sf")

ggplot() +
  geom_sf(data = x, aes(fill= county_effect))+
  scale_fill_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous"))



#https://inla.r-inla-download.org/r-inla.org/doc/latent/besagproper.pdf

#https://faculty.washington.edu/jonno/SISMIDmaterial/3-spatial1.pdf

#We now prepare the BYM model and run INLA
formula <- y ~ 1 + f(ID, model="bym", graph=LDN.adj)
mod <- inla(formula,family="poisson",data=data,E=E)


#We calculate zeta=exp(csi) where csi=upsilon + nu
m <- mod$marginals.random$ID[1:Nareas]
zeta <- lapply(m,function(x)inla.emarginal(exp,x))



#We now calculate the probability that the spatial effects zeta are above 1, 
#identifying areas with excess risk of suicides. This is equivalent to 
#calculate the probability that csi is above 0, which is easier to obtain
a=0
inlaprob<-lapply(mod$marginals.random$ID[1:Nareas], function(X){
  1-inla.pmarginal(a, X)
})

x$excess <- unlist(inlaprob)


ggplot() +
  geom_sf(data = x, aes(fill= excess))





#Now we add the covariates (deprivation - x1 and social fragmentation - x2) and repeat the steps
########################
formula.cov <- y ~ 1+ f(ID, model="bym", graph=LDN.adj) + x1 + x2
mod.cov <- inla(formula.cov,family="poisson",data=data,E=E)
mod.cov$summary.fixed
m <- mod.cov$marginals.random$ID[1:Nareas]

zeta.cov <- lapply(m,function(x)inla.emarginal(exp,x))

a=0

inlaprob.cov<-lapply(mod.cov$marginals.random$ID[1:Nareas], function(X){
  1-inla.pmarginal(a, X)
})

x$excess_cov <- unlist(inlaprob.cov)


ggplot() +
  geom_sf(data = x, aes(fill= excess_cov))


m<-mod.cov$marginals.random$ID
mat.marg<-matrix(NA, nrow=Nareas, ncol=1000)
for (i in 1:Nareas){
  u<-m[[Nareas+i]]
  s<-inla.rmarginal(1000, u)
  mat.marg[i,]<-s}
var.RRspatial<-mean(apply(mat.marg, 2, sd))^2
var.RRhet<-inla.emarginal(function(x) 1/x,
                          mod.cov$marginals.hyper$"Precision for ID (iid component)")
var.RRspatial/(var.RRspatial+var.RRhet)

#Finally we build the maps. First we create a dataset with all the relevant quantities 
#and classes of SMR and posterior probabilities. Then transform the continuous SMR and 
#posterior probabilities in factors, Merge the spatial polygon of London boroughs with the 
#data and map the quantities.
Spatial.results<- data.frame(NAME=data$NAME,SMR=unlist(zeta),
                             pp=unlist(inlaprob), SMR.cov = unlist(zeta.cov), pp.cov = unlist(inlaprob.cov))
SMR.cutoff<- c(0.6, 0.9, 1.0, 1.1,  1.8)
pp.cutoff <- c(0,0.2,0.8,1)
#Transform SMR and pp in factors
SMR.DM=cut(Spatial.results$SMR,breaks=SMR.cutoff,include.lowest=TRUE)
pp.DM=cut(Spatial.results$pp,breaks=pp.cutoff,include.lowest=TRUE)
SMR.COV=cut(Spatial.results$SMR.cov,breaks=SMR.cutoff,include.lowest=TRUE)
pp.COV=cut(Spatial.results$pp.cov,breaks=pp.cutoff,include.lowest=TRUE)
maps.SMR.factors <- data.frame(NAME=data$NAME,
                               SMR.DM=SMR.DM,pp.DM=pp.DM,SMR.COV=SMR.COV,pp.COV=pp.COV)
attr(boroughs, "data")=merge(data.boroughs,maps.SMR.factors,by="NAME")
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "SMR.DM", col.regions=gray(3.5:0.5/4),main="")
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "pp.DM", col.regions=gray(2.5:0.5/3),main="")
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "SMR.COV", col.regions=gray(3.5:0.5/4))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "pp.COV", col.regions=gray(2.5:0.5/3))
###########################################################################################
###########################################################################################


#Spatio-Temporal
H <- inla.read.graph(filename="Ohio_data/Ohio.graph")
image(inla.graph2matrix(H))

data_oh <- read.csv("Ohio_data/OhioRespMort.csv")

data_oh %>% head(5)

data_oh$county1=data_oh$county


formula.one <- y ~ 1 + f(county,model="bym",graph=H,constr=TRUE) + 
  f(county1,year, model="iid",constr=TRUE) + year

model.one <- inla(formula.one,family="poisson", data=data_oh,E=E,
                  control.predictor=list(compute=TRUE),
                  control.compute=list(dic=TRUE))


#Posterior Mean of Spatial Effect

m <- model.one$marginals.random[[1]][1:88]
zeta <- unlist(lapply(m,function(x) inla.emarginal(exp,x)))


x = as(sf::read_sf("Ohio_data/tl_2010_39_county00.shp"), "sf")

x$county_effect = zeta



ggplot() +
  geom_sf(data = x, aes(fill= county_effect)) +
  scale_fill_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous"))


data_oh$year2<-data_oh$year

formula.two <- y ~ 1 + f(county,model="bym",graph=H,constr=TRUE) + 
  f(year,model="rw1") + 
  f(year2,model="iid")


model.two <- inla(formula.two,family="poisson", data=data_oh,E=E,
                  control.predictor=list(compute=TRUE),
                  control.compute=list(dic=TRUE))


model.two$dic$dic
model.one$dic$dic

data_oh$county.year<-seq(1,nrow(data_oh))


formula.three <- y~1+ f(county,model="bym",graph=H,constr=TRUE) + 
  f(year,model="rw1") + 
  f(year2,model="iid")+
  f(county.year,model="iid")



model.three <- inla(formula.three,family="poisson", data=data_oh,E=E,
                  control.predictor=list(compute=TRUE),
                  control.compute=list(dic=TRUE))

model.three$dic$dic



formula.four <- y~1+ f(county,model="bym",graph=H,constr=TRUE) + 
  f(year,model="rw2") + 
  f(year2,model="iid")+
  f(county.year,model="iid")

model.four <- inla(formula.four,family="poisson", data=data_oh,E=E,
                    control.predictor=list(compute=TRUE),
                    control.compute=list(dic=TRUE))

model.four$dic$dic

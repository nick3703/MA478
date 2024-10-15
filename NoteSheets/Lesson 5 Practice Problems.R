library(faraway)
library(tidyverse)
data(gala)
gala %>% str() #Examining Variables

gala_mod <- gala %>% select(-Endemics) #Removing Endemics

#4a
Pois_glm <- glm(Species~Area+Elevation+Nearest+Scruz+Adjacent,data=gala_mod,family=poisson())

Pois_glm$coefficients #Coefficients
deviance(Pois_glm)

#4b
#Done on Paper

#4c

y<-gala_mod$Species

mu<-y
eta <- log(mu)
z <- eta + (y-mu)/mu #Adjusted variables
w <- mu
lmod <- lm(z~Area+Elevation+Nearest+Scruz+Adjacent,data=gala_mod,weights=w)
lmod$coefficients
Pois_glm$coefficients
#Sorta close


#4d
eta <- lmod$fitted.values
mu <- exp(eta)
z <- eta + (y-mu)/mu #Adjusted variables
w <- mu
lmod <- lm(z~Area+Elevation+Nearest+Scruz+Adjacent,data=gala_mod,weights=w)
lmod$coefficients
Pois_glm$coefficients
#Getting Closer!

mu_hat <- exp(lmod$fitted.values)

deviance <- 2*sum(y*log(y/mu_hat)-(y-mu_hat))
deviance - deviance(Pois_glm)
#Deviance is close


#4e
eta <- lmod$fitted.values
mu <- exp(eta)
z <- eta + (y-mu)/mu #Adjusted variables
w <- mu
lmod <- lm(z~Area+Elevation+Nearest+Scruz+Adjacent,data=gala_mod,weights=w)
lmod$coefficients
Pois_glm$coefficients
#Even closer...

mu_hat <- exp(lmod$fitted.values)
deviance <- 2*sum(y*log(y/mu_hat)-(y-mu_hat))
deviance - deviance(Pois_glm)
#Deviance is darn close

#4f

for(j in 1:10){
  eta <- lmod$fitted.values
  mu <- exp(eta)
  z <- eta + (y-mu)/mu #Adjusted variables
  w <- mu
  lmod <- lm(z~Area+Elevation+Nearest+Scruz+Adjacent,data=gala_mod,weights=w)

  mu_hat <- exp(lmod$fitted.values)
  deviance <- 2*sum(y*log(y/mu_hat)-(y-mu_hat))
  print(deviance - deviance(Pois_glm))
}


#Now we look at our coefficients
lmod$coefficients
Pois_glm$coefficients


#Exactly the same

#4g
X <- model.matrix(lmod)

Var_B <- solve(t(X)%*%diag(w)%*%X)

sqrt(diag(Var_B))

#Compare to Standard errors in:

summary(Pois_glm)

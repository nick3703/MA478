library(faraway)
library(tidyverse)
cancer_dat <- read.table("https://users.stat.ufl.edu/~aa/glm/data/Cancer.dat",
                         header=T)

model_1 <- glm(count~factor(histology)+factor(stage)+
                 factor(time),family=poisson(),
               data=cancer_dat)


model_2 <- glm(count~factor(histology)+
                 factor(time),family=poisson(),
               data=cancer_dat)

chi_sq_stat <- as.numeric(-2*(logLik(model_2)-logLik(model_1)))

df <- ncol(model.matrix(model_1))-ncol(model.matrix(model_2))

1-pchisq(chi_sq_stat,df)

anova(model_2,model_1,test="Chisq")

chi_sq_stat


#GOF Test

1-pchisq(deviance(model_1),df.residual(model_1))


#Model with Offset

model_3 <- glm(count~factor(histology)+factor(stage)+
                 factor(time)+offset(log(risktime)),family=poisson(),
               data=cancer_dat)

1-pchisq(deviance(model_3),df.residual(model_3))


cancer_dat %>% ggplot(aes(x=risktime,y=count)) + 
  geom_point()




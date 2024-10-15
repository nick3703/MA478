library(faraway)
library(tidyverse)
data(suicide)
suicide %>% glimpse()
st <- xtabs(y~cause+age+sex,data=suicide)

# Test for independence
summary(st)

#Clearly not independent

modified_dat <- suicide %>%
  mutate(agesex = sex:age) %>%
  select(-age,-sex)

modified_dat %>% glimpse()

mod_t <- xtabs(y~cause+agesex,data=modified_dat)

library(FactoMineR)
CA(mod_t)

mosaicplot(mod_t,color=TRUE,main=NULL,las=1)


#Without considering age:

mod2_t <- suicide %>% group_by(sex,cause) %>%
  summarize(tot=sum(y))


st2 <- xtabs(tot~cause+sex, data=mod2_t)
CA(st2)

test <- glm(tot~cause+sex,data=mod2_t)
z<-xtabs(residuals(test,type="pearson")~cause+sex,data=mod2_t)
blah <- svd(z,2,2)

leftsv <- blah$u %*%diag(sqrt(blah$d))
rightsv <- blah$v %*% diag(sqrt(blah$d))
l1 <- 1.1*max(abs(rightsv),abs(leftsv))

plot(rbind(leftsv,rightsv),asp=1,xlim=c(-l1,l1),
     ylim=c(-l1,l1),xlab="SV1",ylab="SV2",type="n")
abline(h=0,v=0)
text(leftsv,dimnames(z)[[1]])
text(rightsv,dimnames(z)[[2]])


prop.table(st2)

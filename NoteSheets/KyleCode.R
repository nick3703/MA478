dat <- data.frame(n=c(1842,2195,1231),suc=c(1400,1624,943),year=c("a","b","c"))

dat2 <- data.frame(response = c(rep(1,1400),rep(0,1842-1400),rep(1,1624),rep(0,2195-1624),
                     rep(1,943),rep(0,1231-943)),
                   year=c(rep("a",1842),rep("b",2195),rep("c",1231)))


test <- glm(response~year,data=dat2,family=binomial)
summary(test)

confint(test)


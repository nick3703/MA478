library(tidyverse)
library(nnet)
library(VGAM)

data <- read.csv("Course End Feedback Data.csv")

clean_dat <- data %>% 
  mutate(Section=ifelse(section =="MA206    H2 1 - CLARK, N","H","I")) %>%
  mutate(Creative=n.this.course..my.ability.to.think.creatively.and.take.intellectual.risks.increased.) %>%
  select(-section,-n.this.course..my.ability.to.think.creatively.and.take.intellectual.risks.increased.,
         -My.fellow.students.contributed.to.my.learning.in.this.course.) %>%
  mutate(Creative = ordered(Creative,levels=c("Disagree","Neutral","Agree","Strongly Agree")))


my_mod <- vglm(Creative~Section,data=clean_dat,
               family=cumulative(parallel = TRUE))

smaller_mod <- vglm(Creative~1,data=clean_dat,
                    family=cumulative(parallel = TRUE))

chi_stat <- -2*(logLik(smaller_mod)-logLik(my_mod))
1-pchisq(chi_stat,1)

confint(my_mod)
exp(confint(my_mod)[4,])

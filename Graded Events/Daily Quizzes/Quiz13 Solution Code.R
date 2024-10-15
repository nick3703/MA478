library(tidyverse)
library(lme4)
library(nlme)

Chi_dat <- read.csv("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/crime.csv")


Chi_dat <- Chi_dat %>% dplyr::select(-X)

monthly_dat <- apply(Chi_dat,2,sum)


monthly_df <- data.frame(obs=rep(c("A","B","C","D","E","F"),each=12),val=monthly_dat)

plot(monthly_df$val)
monthly_df_mod <- monthly_df %>% mutate(obs_2=seq(1,nrow(monthly_df)))

model_1 <- lme(val~1,
               random=~1|obs,          
               correlation=corCompSymm(form=~1|obs),
                         data=monthly_df,
                         method="ML")

model_2  <- lme(val~1,
                          random=~1|obs,
                          correlation=corAR1(form=~1|obs),
                          data=monthly_df,
                          method="ML")


model_3  <- lmer(val~1,
                random=~1|obs_2,
                data=monthly_df_mod)


getVarCov(model_2,obs=1,type="marginal")
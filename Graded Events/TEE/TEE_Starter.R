### MA478 Generalized Linear Models

### Term End Analysis Competition (TEE) Starter File

## OBJECTIVE 

# A charitable organization wishes to develop generalized linear models (GLM)
# to improve the cost-effectiveness of their direct marketing campaigns
# to previous donors.

# 1) Develop a GLM using data from the most recent campaign that
# can effectively capture likely donors so that the expected net profit is maximized.

# 2) Develop a GLM to predict donation amounts for donors - the data
# for this will consist of the records for donors only.


## DATA EXPLORATION

# Load the data
charity <- read.csv(file.choose()) # load the "charity_full.csv" file
# conduct additional exploratory data analysis


## DATA PREPARATION

# Predictor transformation example
charity.t <- charity
charity.t$avhv <- log(charity.t$avhv)
# add further transformations and data preparations, where appropriate
# for example, some GLM methods can struggle when predictors are highly skewed

# Setup training, validation and testing data subsets for analysis
data.train <- charity.t[charity$part=="train",]
x.train <- data.train[,2:21]
c.train <- data.train[,22] # donr
n.train.c <- length(c.train) # 3984
y.train <- data.train[c.train==1,23] # damt for observations with donr=1
n.train.y <- length(y.train) # 1995

data.valid <- charity.t[charity$part=="valid",]
x.valid <- data.valid[,2:21]
c.valid <- data.valid[,22] # donr
n.valid.c <- length(c.valid) # 2018
y.valid <- data.valid[c.valid==1,23] # damt for observations with donr=1
n.valid.y <- length(y.valid) # 999

data.test <- charity.t[charity$part=="test",]
n.test <- dim(data.test)[1] # 2007
x.test <- data.test[,2:21]

# Standardize the predictor variables
x.train.mean <- apply(x.train, 2, mean)
x.train.sd <- apply(x.train, 2, sd)
x.train.std <- t((t(x.train)-x.train.mean)/x.train.sd) # standardize to have zero mean and unit sd
apply(x.train.std, 2, mean) # check zero mean
apply(x.train.std, 2, sd) # check unit sd
data.train.std.c <- data.frame(x.train.std, donr=c.train) # to classify donr
data.train.std.y <- data.frame(x.train.std[c.train==1,], damt=y.train) # to predict damt when donr=1

x.valid.std <- t((t(x.valid)-x.train.mean)/x.train.sd) # standardize using training mean and sd
data.valid.std.c <- data.frame(x.valid.std, donr=c.valid) # to classify donr
data.valid.std.y <- data.frame(x.valid.std[c.valid==1,], damt=y.valid) # to predict damt when donr=1

x.test.std <- t((t(x.test)-x.train.mean)/x.train.sd) # standardize using training mean and sd
data.test.std <- data.frame(x.test.std)


## GLM MODELING FOR CLASSIFICATION

# EXAMPLE: Binary logistic regression
model.log1 <- glm(donr ~ reg1 + reg2 + reg3 + reg4 + home + chld + hinc + I(hinc^2) + genf + wrat + 
                    avhv + incm + inca + plow + npro + tgif + lgif + rgif + tdon + tlag + agif, 
                  data.train.std.c, family=binomial("logit"))

post.valid.log1 <- predict(model.log1, data.valid.std.c, type="response") # n.valid post probs

# Calculate ordered profit function using average donation = $14.50 and mailing cost = $2
profit.log1 <- cumsum(14.5*c.valid[order(post.valid.log1, decreasing=T)]-2)
plot(profit.log1) # see how profits change as more mailings are made
n.mail.valid <- which.max(profit.log1) # number of mailings that maximizes profits
c(n.mail.valid, max(profit.log1)) # report number of mailings and maximum profit
# 1291.0 11642.5

cutoff.log1 <- sort(post.valid.log1, decreasing=T)[n.mail.valid+1] # set cutoff based on n.mail.valid
chat.valid.log1 <- ifelse(post.valid.log1>cutoff.log1, 1, 0) # mail to everyone above the cutoff
table(chat.valid.log1, c.valid) # classification table
#               c.valid
#chat.valid.log1   0   1
#              0 709  18
#              1 310 981
# check n.mail.valid = 310+981 = 1291
# check profit = 14.5*981-2*1291 = 11642.5

# Results

# n.mail Profit  Model
# 1291   11642.5 Log1

# Select best model (ex: model.log1) with the maximum profit in the validation sample
post.test <- predict(model.log1, data.test.std, type="response") # post probs for test data

# Oversampling adjustment for calculating number of mailings for test set
n.mail.valid <- which.max(profit.log1)
tr.rate <- .1 # typical response rate is .1
vr.rate <- .5 # whereas validation response rate is .5
adj.test.1 <- (n.mail.valid/n.valid.c)/(vr.rate/tr.rate) # adjustment for mail yes
adj.test.0 <- ((n.valid.c-n.mail.valid)/n.valid.c)/((1-vr.rate)/(1-tr.rate)) # adjustment for mail no
adj.test <- adj.test.1/(adj.test.1+adj.test.0) # scale into a proportion
n.mail.test <- round(n.test*adj.test, 0) # calculate number of mailings for test set

cutoff.test <- sort(post.test, decreasing=T)[n.mail.test+1] # set cutoff based on n.mail.test
chat.test <- ifelse(post.test>cutoff.test, 1, 0) # mail to everyone above the cutoff
table(chat.test)
#    0    1 
# 1676  331
# based on this model we'll mail to the 331 highest posterior probabilities

# See below for saving chat.test into a file for submission


## GLM MODELING FOR PREDICTION

# EXAMPLE: Linear regression
model.ls1 <- lm(damt ~ reg1 + reg2 + reg3 + reg4 + home + chld + hinc + genf + wrat + 
                  avhv + incm + inca + plow + npro + tgif + lgif + rgif + tdon + tlag + agif, 
                data.train.std.y)

pred.valid.ls1 <- predict(model.ls1, newdata = data.valid.std.y) # validation predictions
mean((y.valid - pred.valid.ls1)^2) # mean squared error
# 1.867523
sd((y.valid - pred.valid.ls1)^2)/sqrt(n.valid.y) # std error
# 0.1696615

# Results

# MSE       Model
# 1.867523  LS1

# Select best model (ex: model.ls1) since it has minimum MSE in the validation sample
yhat.test <- predict(model.ls1, newdata = data.test.std) # test predictions


## FINAL RESULTS

# Save final results for both best GLM models for classification and regression
length(chat.test) # check length = 2007
length(yhat.test) # check length = 2007
chat.test[1:10] # check this consists of 0s and 1s
yhat.test[1:10] # check this consists of plausible predictions of damt

# Create a data frame with three variables: ID, DONR and DAMT
ip <- data.frame(ID=data.test$ID, donr=chat.test, damt=yhat.test)

# Where donr = 0 make the value of damt = 0
ip$damt[ip$donr == 0] <- 0

# Compute estimated profits to test set mailings to the most likely donors
idx.donors <- ip[(ip$donr > 0), ]
num_donors <- nrow(idx.donors)
est.profits <- sum(idx.donors$damt) - 2*num_donors
round(est.profits,2)

submit <- data.frame(ID=ip$ID, damt=ip$damt) # data frame with two variables: ID and DAMT
write.csv(submit, file="ABC_submission.csv", row.names=FALSE) # use your initials_submission for the file name

# Submit the csv file on Kaggle for evaluation
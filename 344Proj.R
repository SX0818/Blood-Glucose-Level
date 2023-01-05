library(tidyr)
library(dplyr)
library(formattable)
library(knitr)

# data <- healthcare.dataset.stroke.data
data <- read.csv("/Users/yuqianxie/Desktop/healthcare-dataset-stroke-data.csv",header=T)
# data cleansing (remove row with NA value)
data[data == "N/A"] <- NA
data[data == "Unknown"] <- NA
data <- data |> drop_na()
data <- data[-c(1:3)]
data <- data[-c(3:5)]

# Set up
set.seed(20)
data$bmi <- as.numeric(data$bmi)
N <- length(data$avg_glucose_level) # population size = 3426

n <- 346 # sample size = 346 (based on study planning)

# SRS
# Continuous (Vanilla estimator)
ybar.pop <- mean(data$avg_glucose_level) # true value
SRS.indices <- sample.int(N,n,replace = F)
SRS.sample <- data[SRS.indices, ]
# attach(SRS.sample)
vanilla.estimator <- mean(SRS.sample$avg_glucose_level)

se.function  <-  function(sample.value, estimated.value) {
  res <- sample.value - estimated.value  # residual
  temp <- sum(res^2)/(n-1)
  se <- sqrt((1-n/N) *(temp/n))
  return (se)
}

vanilla.se  <- se.function(SRS.sample$avg_glucose_level, vanilla.estimator)
# vanilla.se <- sqrt((1 - n / N) * var(SRS.sample$avg_glucose_level) / n)
vanilla.srs <- c(vanilla.estimator, vanilla.se)
vanilla.lower.limit.CI <- vanilla.estimator - 1.96*vanilla.se
vanilla.upper.limit.CI <- vanilla.estimator + 1.96*vanilla.se

# Continuous (Ratio estimator)
# Note: we use bmi as our auxiliary variable
xbar.pop <- mean(data$bmi)

ratio.estimator <- (mean(SRS.sample$avg_glucose_level) /mean(SRS.sample$bmi))*xbar.pop
ratio.se  <- se.function(SRS.sample$avg_glucose_level, (mean(SRS.sample$avg_glucose_level) /mean(SRS.sample$bmi))*SRS.sample$bmi)
ratio.srs <- c(ratio.estimator,ratio.se)
ratio.lower.limit.CI <- ratio.estimator - 1.96*ratio.se
ratio.upper.limit.CI <- ratio.estimator + 1.96*ratio.se
cor(bmi,avg_glucose_level)


# Table for SRS Continuous Data
tab <- matrix(c(vanilla.estimator,ratio.estimator,vanilla.se,ratio.se,
vanilla.lower.limit.CI,ratio.lower.limit.CI,vanilla.upper.limit.CI, ratio.upper.limit.CI)
,nrow=2,ncol=4)
colnames(tab) <- c('estimates','standard erros','lower confidence intervals','upper confidence intervals')
rownames(tab) <- c('Vanilla','Ratio')
tab <- as.table(tab)
formattable(tab) %>% kable(caption = "SRS_Continuous_Data")

# Some comments:
# We know that the standard error for the vanilla estimator is 1.836228;
# the standard error for the ratio estimator is 2.026795.
# Check the correlation value is 0.1618729; we learned that the ratio beats the vanilla only if
# there is a strong correlation between X and Y, so our output is reasonable.

# SRS Binary
# Note: we are interested in patients' glucose_level whose below 100
count <- 0
SRS.sample$avg_glucose_level <- as.array(SRS.sample$avg_glucose_level)
for (i in 1:n) {
  if (SRS.sample$avg_glucose_level[i] < 100)
    count <- count + 1
}
sample.value <- count
prop.estimator <-  sample.value / n # sample proportion (estimator)
prop.se <- sqrt((prop.estimator*(1-prop.estimator))/n)
prop.lower.limit.CI <- prop.estimator - 1.96*prop.se
prop.upper.limit.CI <- prop.estimator + 1.96*prop.se
tab3 <- matrix(c(prop.estimator,prop.se,prop.lower.limit.CI,prop.upper.limit.CI)
              ,nrow=1,ncol=4)
colnames(tab3) <- c('estimates','standard erros','lower confidence intervals','upper confidence intervals')
rownames(tab3) <- c('SRS (binary)')
tab3 <- as.table(tab3)
formattable(tab3) %>% kable(caption = "SRS_Binary_Data")
# detach(SRS.sample)

# STR continuous
# Note: we are interested in estimating the average glucose level based on smoking status.
N.h <- tapply(data$avg_glucose_level, data$smoking_status, length)   # population size for peoples' smoking status
smoke <- names(N.h) # smoking status

# Estimate the population mean using STR with proportional allocation
STR.sample.prop <- NULL
n.h.prop <- round( (N.h/N) * n)
for (i in 1: length(smoke)) {
row.indices <- which(data$smoking_status == smoke[i])
sample.indices <- sample(row.indices, n.h.prop[i], replace = F)
STR.sample.prop <- rbind(STR.sample.prop, data[sample.indices, ])
}
ybar.h.prop <- tapply(STR.sample.prop$avg_glucose_level, STR.sample.prop$smoking_status, mean)
var.h.prop <- tapply(STR.sample.prop$avg_glucose_level, STR.sample.prop$smoking_status, var)
se.h.prop <- sqrt((1 - n.h.prop / N.h) * var.h.prop / n.h.prop)
# rbind(ybar.h.prop, se.h.prop)
ybar.str.prop <- sum(N.h / N * ybar.h.prop)
se.str.prop <- sqrt(sum((N.h / N)^2 * se.h.prop^2))
str.prop <- c(ybar.str.prop, se.str.prop)
STR.CI <- c(ybar.str.prop - 1.96*se.str.prop , ybar.str.prop + 1.96*se.str.prop)
STR.lower.limit.CI <- ybar.str.prop - 1.96*se.str.prop
STR.upper.limit.CI <- ybar.str.prop + 1.96*se.str.prop

# Table for STR Continuous Data
tab2 <- matrix(c(ybar.str.prop,se.str.prop,STR.lower.limit.CI,STR.upper.limit.CI)
              ,nrow=1,ncol=4)
colnames(tab2) <- c('estimates','standard erros','lower confidence intervals','upper confidence intervals')
rownames(tab2) <- c('STR (continuous)')
tab2 <- as.table(tab2)
formattable(tab2) %>% kable(caption = "STR_Continuous_Data")

# STR binary
ybar.prop <- tapply(as.numeric(STR.sample.prop$avg_glucose_level < 100), STR.sample.prop$smoking_status, mean)
var.prop <- ybar.prop * ( 1- ybar.prop)
se.prop <- sqrt((1 - n.h.prop / N.h) * var.prop / n.h.prop)
ybar.prop <- sum(N.h / N * ybar.prop)
se.new.prop <- sqrt(sum((N.h / N)^2 * se.prop^2))
str.new.prop <- c(ybar.prop, se.new.prop)
STR.new.lower.limit.CI <- ybar.prop - 1.96*se.new.prop
STR.new.upper.limit.CI <- ybar.prop + 1.96*se.new.prop

# Table for STR Binary Data
tab4 <- matrix(c(ybar.prop,se.new.prop,STR.new.lower.limit.CI,STR.new.upper.limit.CI)
               ,nrow=1,ncol=4)
colnames(tab4) <- c('estimates','standard erros','lower confidence intervals','upper confidence intervals')
rownames(tab4) <- c('STR (binary)')
tab4 <- as.table(tab4)
formattable(tab4) %>% kable(caption = "STR_Binary_Data")








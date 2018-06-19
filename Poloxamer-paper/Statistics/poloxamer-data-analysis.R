# library(lme4)

library(nlme)

PBS_data = read.csv("/Users/michelle/Documents/2018_CUMC/Research/PBS.csv")

head(PBS_data)

table(PBS_data$experiment, exclude=NULL)
table(PBS_data$time, exclude=NULL)
table(PBS_data$numHoles, exclude=NULL)
table(PBS_data$holeDiameter, exclude=NULL)
table(PBS_data$holeArea, exclude=NULL)

with(PBS_data, table(holeDiameter, numHoles, exclude=NULL))

# study2 indicator of 1 means different number of holes, same diameter of holes
PBS_data$study2 = as.numeric(PBS_data$holeDiameter == 0 | PBS_data$holeDiameter == 100)

with(PBS_data, table(study2, holeDiameter, exclude=NULL))

fit2 = lme(Concentration ~ factor(numHoles) + time + factor(numHoles):time, random = ~1|experiment, data = PBS_data[which(PBS_data$study2 == 1), ])

summary(fit2)

# global likelihood ratio test
fit2_1 = lme(Concentration ~ factor(numHoles) + time + factor(numHoles):time, random = ~1|experiment, data = PBS_data[which(PBS_data$study2 == 1), ], method="ML")
fit2_0 = lme(Concentration ~ factor(numHoles) + time, random = ~1|experiment, data = PBS_data[which(PBS_data$study2 == 1), ], method="ML")
#compare the two models
anova(fit2_1, fit2_0)

# fit linear mixed model to account for correlation of observations from the same experiment.  outcome was Concentration.  The fixed effects were number of Holes, treated as categorical, time, and interactions between number of holes and time.  A random intercept was included for experiment.  A global test for whether there was an interaction between number of holes and time was examined by using a likelihood ratio test, comparing the previously mentioned model to one without interaction terms.

# Study 1 - PBS data only
# study1 indicator of 1 means different number of holes, different diameter of holes, same area of combined holes.
PBS_data$study1 = as.numeric(PBS_data$numHoles == 0 | PBS_data$holeDiameter == 200 | PBS_data$numHoles == 4)
with(PBS_data, table(study1, holeDiameter, exclude=NULL))
with(PBS_data, table(study1, numHoles, exclude=NULL))

fit1 = lme(Concentration ~ factor(numHoles) + time + factor(numHoles):time, random = ~1|experiment, data = PBS_data[which(PBS_data$study1 == 1), ])

summary(fit1)

# global likelihood ratio test
fit1_1 = lme(Concentration ~ factor(numHoles) + time + factor(numHoles):time, random = ~1|experiment, data = PBS_data[which(PBS_data$study1 == 1), ], method="ML")
fit1_0 = lme(Concentration ~ factor(numHoles) + time, random = ~1|experiment, data = PBS_data[which(PBS_data$study1 == 1), ], method="ML")
#compare the two models
anova(fit1_1, fit1_0)



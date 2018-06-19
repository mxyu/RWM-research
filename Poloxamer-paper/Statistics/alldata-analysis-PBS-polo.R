library(gdata)
library(nlme)
library(ggplot2)
# read data
data_description = read.xls("/Users/michelle/Documents/2018_CUMC/Research/data_PBSandPOLOXAMER.xlsx", sheet=1)
all_data = read.xls("/Users/michelle/Documents/2018_CUMC/Research/data_PBSandPOLOXAMER.xlsx", sheet=2)
all_data$polo_g1 <- as.numeric(all_data$ExpType=="H")
all_data$polo_g4 <- as.numeric(all_data$ExpType=="G")
all_data$time_sq <- all_data$Time^2
all_data$timexpolo_g1 <- all_data$Time * all_data$polo_g1
all_data$timexpolo_g4 <- all_data$Time * all_data$polo_g4
all_data$time_sqxpolo_g1 <- all_data$time_sq * all_data$polo_g1
all_data$time_sqxpolo_g4 <- all_data$time_sq * all_data$polo_g4

polo_data <- subset(all_data, ExpType %in% c('F','G','H'))
fit_output = lme(ConcAdjustedMm2 ~ 1 + polo_g1 + polo_g4 + Time + time_sq + timexpolo_g1 + timexpolo_g4 + time_sqxpolo_g1 + time_sqxpolo_g4, random = ~1|ExperimentID, data = polo_data)
show(summary(fit_output))
# global likelihood ratio test
fit_1 = lme(ConcAdjustedMm2 ~ 1 + polo_g1 + polo_g4 + Time + time_sq + timexpolo_g1 + timexpolo_g4 + time_sqxpolo_g1 + time_sqxpolo_g4, random = ~1|ExperimentID, data = polo_data, method="ML")
fit_0 = lme(ConcAdjustedMm2 ~ 1 + polo_g1 + polo_g4 + Time + time_sq + timexpolo_g1 + time_sqxpolo_g1, random = ~1|ExperimentID, data = polo_data, method="ML")
#compare the two models
anova(fit_1, fit_0)

# Function that will take argument of a list of Experiment Types (A-H) and print a summary of the lme model and global likelihood ratio test comparing those experiments
run_stats <- function(x){
  all_data$studyTmp = as.numeric(all_data$ExpType %in% x)
  fit_output = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time + factor(NumHoles):Time, random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ])
  show(summary(fit_output))
  # global likelihood ratio test
  fit_1 = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time + factor(NumHoles):Time, random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ], method="ML")
  fit_0 = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time, random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ], method="ML")
  #compare the two models
  anova(fit_1, fit_0)
}

run_stats_sq <- function(x){
  all_data$studyTmp = as.numeric(all_data$ExpType %in% x)
  fit_output = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time + I(Time^2) + factor(NumHoles):Time + factor(NumHoles):I(Time^2), random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ])
  show(summary(fit_output))
  # global likelihood ratio test
  fit_1 = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time + I(Time^2) + factor(NumHoles):Time + factor(NumHoles):I(Time^2), random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ], method="ML")
  fit_0 = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time + I(Time^2) + factor(NumHoles):Time, random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ], method="ML")
  #compare the two models
  anova(fit_1, fit_0)
}

run_stats_traj <- function(x){
  all_data$studyTmp = as.numeric(all_data$ExpType %in% x)
  fit_output = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time + I(Time^2) + factor(NumHoles):Time + factor(NumHoles):I(Time^2), random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ])
  show(summary(fit_output))
  # global likelihood ratio test
  fit_1 = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time + I(Time^2) + factor(NumHoles):Time + factor(NumHoles):I(Time^2), random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ], method="ML")
  fit_0 = lme(ConcAdjustedMm2 ~ factor(NumHoles) + Time + I(Time^2), random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ], method="ML")
  #compare the two models
  anova(fit_1, fit_0)
}

# Using experiments typed A-D, comparing different numbers of holes w/ same diameter in PBS.
run_stats(c('A','B','C','D'))

# Experiments A, D, and E, comparing different number of holes with same total hole area (different diameters) in PBS.
run_stats(c('A', 'D', 'E'))

# Experiments F, G, H; comparing different number of holes with same total hole area (different diameters) in 18% Poloxamer
run_stats(c('F','G','H'))
run_stats_sq(c('F','G','H'))
run_stats_traj(c('F','G'))

# Function that will take argument of a list of Experiment Types (A-H) and print a summary of the lme model and global likelihood ratio test comparing those experiments. Model factors categorically by experiment type rather than numHoles, which is how run_stats function works. Allows for comparison of experiments where numHoles is the same between different experimental groups.
run_stats_by_expType <- function(x){
  all_data$studyTmp = as.numeric(all_data$ExpType %in% x)
  fit_output = lme(ConcAdjustedMm2 ~ factor(ExpType) + Time + factor(ExpType):Time, random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ])
  show(summary(fit_output))
  # global likelihood ratio test
  fit_1 = lme(ConcAdjustedMm2 ~ factor(ExpType) + Time + factor(ExpType):Time, random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ], method="ML")
  fit_0 = lme(ConcAdjustedMm2 ~ factor(ExpType) + Time, random = ~1|ExperimentID, data = all_data[which(all_data$studyTmp == 1), ], method="ML")
  #compare the two models
  anova(fit_1, fit_0)
  # return lme
  return(fit_output)
}

# Compare 0 holes PBS and 0 holes Poloxamer
lme_fit = run_stats_by_expType(c('F','A'))
lme_fit = run_stats_by_expType(c('A', 'B', 'C', 'D','E'))
lme_fit = run_stats_by_expType(c('E', 'H'))

#ggplot(all_data[which(all_data$ExpType %in% c('A','B','C','D')),], aes(x=Time,y=ConcAdjM2,group=ExpType,col=ExpType)) + 
#  geom_smooth(method="lm", se=T) #+
  #geom_point(alpha = 0.3)
  # geom_hline(yintercept=0, linetype="dashed") +
  # theme_bw()

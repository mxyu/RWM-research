# read data

all_data <- read.xls("/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/SSI_stats_reformatted.xlsx", sheet=1)
head(all_data)
names(all_data)
all_data$no_ssi <- with(all_data, procedure_count - All_SSI_Model_Infection_Count)
all_data$prop_ssi <- with(all_data, All_SSI_Model_Infection_Count / procedure_count)

# Get properties of the data - break down by site, procedure, year
# How many per site over each year. Both total # procedures and # SSI
usite_all <- unique(all_data$Site)
mat_site_proc_count <- matrix(NA, length(usite_all), 3)
mat_site_SSI_count <- matrix(NA, length(usite_all), 3)
for (i in 1:length(usite_all)) {
  num_procedures_2015 <- sum(with(all_data, procedure_count[Site == usite_all[i] & summary_year == 2015]), na.rm=TRUE)
  num_procedures_2016 <- sum(with(all_data, procedure_count[Site == usite_all[i] & summary_year == 2016]), na.rm=TRUE)
  num_procedures_2017 <- sum(with(all_data, procedure_count[Site == usite_all[i] & summary_year == 2017]), na.rm=TRUE)
  num_SSI_2015 <-sum(with(all_data, All_SSI_Model_Infection_Count[Site == usite_all[i] & summary_year == 2015]), na.rm=TRUE)
  num_SSI_2016 <-sum(with(all_data, All_SSI_Model_Infection_Count[Site == usite_all[i] & summary_year == 2016]), na.rm=TRUE)
  num_SSI_2017 <-sum(with(all_data, All_SSI_Model_Infection_Count[Site == usite_all[i] & summary_year == 2017]), na.rm=TRUE)
  mat_site_proc_count[i,1] <- num_procedures_2015
  mat_site_proc_count[i,2] <- num_procedures_2016
  mat_site_proc_count[i,3] <- num_procedures_2017
  mat_site_SSI_count[i,1] <- num_SSI_2015
  mat_site_SSI_count[i,2] <- num_SSI_2016
  mat_site_SSI_count[i,3] <- num_SSI_2017
}
colnames(mat_site_proc_count) <- c("2015", "2016", "2017")
rownames(mat_site_proc_count) <- usite_all
mat_site_proc_prop <- prop.table(mat_site_proc_count, 2)
mat_site_proc_prop_percent <- 100*round(mat_site_proc_prop, digits=3)
colnames(mat_site_SSI_count) <- c("2015", "2016", "2017")
rownames(mat_site_SSI_count) <- usite_all
mat_site_SSI_prop <- prop.table(mat_site_SSI_count, 2)
mat_site_SSI_prop_percent <- 100*round(mat_site_SSI_prop, digits=3)


# Div by procedure type over each year. Both total # procedures and # SSI
uproc_all <- unique(all_data$Procedure_type)
mat_ptype_proc_count <- matrix(NA, length(uproc_all), 3)
mat_ptype_SSI_count <- matrix(NA, length(uproc_all), 3)
for (i in 1:length(uproc_all)) {
  num_procedures_2015 <- sum(with(all_data, procedure_count[Procedure_type == uproc_all[i] & summary_year == 2015]), na.rm=TRUE)
  num_procedures_2016 <- sum(with(all_data, procedure_count[Procedure_type == uproc_all[i] & summary_year == 2016]), na.rm=TRUE)
  num_procedures_2017 <- sum(with(all_data, procedure_count[Procedure_type == uproc_all[i] & summary_year == 2017]), na.rm=TRUE)
  num_SSI_2015 <-sum(with(all_data, All_SSI_Model_Infection_Count[Procedure_type == uproc_all[i] & summary_year == 2015]), na.rm=TRUE)
  num_SSI_2016 <-sum(with(all_data, All_SSI_Model_Infection_Count[Procedure_type == uproc_all[i] & summary_year == 2016]), na.rm=TRUE)
  num_SSI_2017 <-sum(with(all_data, All_SSI_Model_Infection_Count[Procedure_type == uproc_all[i] & summary_year == 2017]), na.rm=TRUE)
  mat_ptype_proc_count[i,1] <- num_procedures_2015
  mat_ptype_proc_count[i,2] <- num_procedures_2016
  mat_ptype_proc_count[i,3] <- num_procedures_2017
  mat_ptype_SSI_count[i,1] <- num_SSI_2015
  mat_ptype_SSI_count[i,2] <- num_SSI_2016
  mat_ptype_SSI_count[i,3] <- num_SSI_2017
}
colnames(mat_ptype_proc_count) <- c("2015", "2016", "2017")
rownames(mat_ptype_proc_count) <- uproc_all
mat_ptype_proc_prop <- prop.table(mat_ptype_proc_count, 2)
mat_ptype_proc_prop_percent <- 100*round(mat_ptype_proc_prop, digits=3)
colnames(mat_ptype_SSI_count) <- c("2015", "2016", "2017")
rownames(mat_ptype_SSI_count) <- uproc_all
mat_ptype_SSI_prop <- prop.table(mat_ptype_SSI_count, 2)
mat_ptype_SSI_prop_percent <- 100*round(mat_ptype_SSI_prop, digits=3)

# Print to file
write.csv(mat_site_proc_count, "/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/proc_ct_by_site_year.csv")
write.csv(mat_site_SSI_count, "/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/SSI_ct_by_site_year.csv")
write.csv(mat_site_proc_prop_percent, "/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/proc_prop_by_site_year.csv")
write.csv(mat_site_SSI_prop_percent, "/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/SSI_prop_by_site_year.csv")

write.csv(mat_ptype_proc_count, "/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/proc_ct_by_proctype_year.csv")
write.csv(mat_ptype_SSI_count, "/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/SSI_ct_by_proctype_year.csv")
write.csv(mat_ptype_proc_prop_percent, "/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/proc_prop_by_proctype_year.csv")
write.csv(mat_ptype_SSI_prop_percent, "/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/SSI_prop_by_proctype_year.csv")



# Exclude 2016 data
data1517 <- subset(all_data, summary_year != 2016)

with(data1517, procedure_count[Site == "AH"])

usite <- unique(data1517$Site)
mat_site = matrix(NA, length(usite), 2)

for (i in 1:length(usite)) {
  num_procedures_2015 <- sum(with(data1517, procedure_count[Site == usite[i] & summary_year == 2015]))
  num_procedures_2017 <- sum(with(data1517, procedure_count[Site == usite[i] & summary_year == 2017]))
  # num_SSI <-sum(with(data1517, All_SSI_Model_Infection_Count[Site == usite[i]]))
  # num_no_SSI <- num_procedures - num_SSI;
  mat_site[i,1] <- num_procedures_2015
  mat_site[i,2] <- num_procedures_2017
}
colnames(mat_site) <- c("2015", "2017")
rownames(mat_site) <- usite

mat_site
prop.table(mat_site, 2)
chisq.test(mat_site)
fisher.test(mat_site)



# Procedure type
uproc <- unique(data1517$Procedure_type)
mat_proc = matrix(NA, length(uproc), 2)

for (i in 1:length(uproc)) {
  num_procedures_2015 <- sum(with(data1517, procedure_count[Procedure_type == uproc[i] & summary_year == 2015]))
  num_procedures_2017 <- sum(with(data1517, procedure_count[Procedure_type == uproc[i] & summary_year == 2017]))
  # num_SSI <-sum(with(data1517, All_SSI_Model_Infection_Count[Site == usite[i]]))
  # num_no_SSI <- num_procedures - num_SSI;
  mat_proc[i,1] <- num_procedures_2015
  mat_proc[i,2] <- num_procedures_2017
}
colnames(mat_proc) <- c("2015", "2017")
rownames(mat_proc) <- uproc

mat_proc
prop.table(mat_proc, 2)
chisq.test(mat_proc)
fisher.test(mat_proc)

# data_without_cardiac <- read.xls("/Users/michelle/Documents/2018_CUMC/Research/SSI_rates/SSI_stats_reformatted.xlsx", sheet=1)
# num_procedures_2015 <- sum(as.numeric(data_without_cardiac$procedure_count[data_without_cardiac$summary_year==2015]), na.rm=TRUE)
# num_procedures_2017 <- sum(as.numeric(data_without_cardiac$procedure_count[data_without_cardiac$summary_year==2017]), na.rm=TRUE)
# total_procedures <- num_procedures_2015+num_procedures_2017


with(all_data, table(summary_year, exclude=NULL))
num_procedures_2015 <- sum(as.numeric(all_data$procedure_count[all_data$summary_year==2015]), na.rm=TRUE)
num_procedures_2016 <- sum(as.numeric(all_data$procedure_count[all_data$summary_year==2016]), na.rm=TRUE)
num_procedures_2017 <- sum(as.numeric(all_data$procedure_count[all_data$summary_year==2017]), na.rm=TRUE)
total_procedures <- num_procedures_2015+num_procedures_2017

#pwr.chisq.test(N = total_procedures, df = 1, sig.level=0.05, power=0.8)

SSI_2015 <- sum(as.numeric(all_data$All_SSI_Model_Infection_Count[all_data$summary_year==2015]), na.rm=TRUE)
SSI_2016 <- sum(as.numeric(all_data$All_SSI_Model_Infection_Count[all_data$summary_year==2016]), na.rm=TRUE)
SSI_2017 <- sum(as.numeric(all_data$All_SSI_Model_Infection_Count[all_data$summary_year==2017]), na.rm=TRUE)
procedure_no_SSI_2015 <- num_procedures_2015-SSI_2015
procedure_no_SSI_2017 <- num_procedures_2017-SSI_2017

mat <- matrix(c(SSI_2015, procedure_no_SSI_2015, SSI_2017, procedure_no_SSI_2017), nrow=2, ncol=2)
# data <- as.data.frame(mat, row.names=c("SSI","no_SSI"))
# colnames(data) <- c("2015","2017")                      
# chisq.test(data)
# fisher.test(data)
# prop.table(as.matrix(data), 2)

colnames(mat) <- c("2015","2017")                      

mat
prop.table(mat, 2)
chisq.test(mat)
fisher.test(mat)


# How much improvement in SSI can we detect given our sample size? 

power.prop.test(n=6614, p2=0.01, power=0.80)
power.prop.test(n=5779, p2=0.01, power=0.80)


# Would like to determine what other factors may significantly affect SSI rate
# e.g. procedure type, site

# logistic regression
# glm(, family=binomial(link = "logit"))


fit1 <- glm(cbind(All_SSI_Model_Infection_Count, no_ssi) ~ factor(summary_year), family=binomial(link = "logit"), data=data1517)

# does NOT work
# fit2 <- glm(All_SSI_Model_Infection_Count ~ factor(summary_year), weights=procedure_count, family=binomial(link = "logit"), data=data1517)

fit2 <- glm(prop_ssi ~ factor(summary_year), weights=procedure_count, family=binomial(link = "logit"), data=data1517)

summary(fit1) # coef estimate for year variable is log(odds ratio)
summary(fit2)

full_fit <- glm(cbind(All_SSI_Model_Infection_Count, no_ssi) ~ factor(summary_year) + Procedure_type + Site, family=binomial(link = "logit"), data=data1517)

summary(full_fit)

exp(coef(fit1)) # unadjusted odds ratio
exp(coef(full_fit)) # odds ratio adjusted for site and procedure type


# How to use SIR
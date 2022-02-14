# import packages
if (!require(c("tidyverse", "ggalt","gdata","gridExtra","rgenoud","MatchIt","plm","Hmisc","lmtest","sandwich","multiwayvcov")))
  install.packages(c("tidyverse","ggalt","gdata","gridExtra","rgenoud","MatchIt","plm","Hmisc","lmtest","sandwich","multiwayvcov"))

library(tidyverse)
library(ggalt)
library(gdata)
library(gridExtra)
library(rgenoud)
library(MatchIt)
library(plm)
library(Hmisc)
library(lmtest)
library(sandwich)
library(multiwayvcov)

## data import

ELC <- read.csv(...) #5,925,488 - residential electric usage 2004-2019
ELC <- ELC %>% 
  mutate(yearELC = ..., #factor variable
         monthELC = ...) #factor variable
PropertyStats <- read.csv(...) #16,680 - property data

## fixed effects model

# Treatment variable encompass both the fact of being treated and the time of treatment => 'Treatment' equals 1 when the treated unit gets treatment
ELC_property <- merge(ELC, PropertyStats, by = 'AddressIndex') #2,931,700 - merge two datasets
ELC_property <- ELC_property %>%
  mutate(NormConsumption = Consumption/size, #standardize consumption
         Treatment = ifelse(InitialPeriod<=Period,1,0)) #create treatment variable
ELC_property <- ELC_property %>% mutate(Treatment = ifelse(is.na(Treatment),0,Treatment))
ELC_property <- ELC_property %>% filter(NormConsumption <= quantile(NormConsumption, c(0.9999), na.rm = TRUE)) #remove outliers

# regression adjustments
fe_reg <- plm(log(NormConsumption) ~ Treatment + yearELC + monthELC + CoolingDays + HeatingDays, data = ELC_property, model='within', index = c('ID','Period')) #fixed effects to remove unobserved heterogeneity between properties and periods
summary(fe_reg) 
# Treatment - a dummy variable indicating program participation, 
# yearELC & monthELC - time effects,
# CoolingDays & HeatingDays - monthly number of either days retrieved from National Oceanic and Atmospheric Administration (NOAA), 
# model 'within' - fixed effects model
summary(fe_reg)  #call regression output

## PSM - Propensity Score Matching

# propensity scores
psm_match <- matchit(Group ~ ... [your covariates] ...,
                     method='nearest', data=PropertyStats, replace = TRUE, ratio=21)
# method 'nearest' indicates propensity score matching method (nearest neighbor matching), 
# ratio - number of control matches per treatment unit,
# replace = TRUE - multiple time usage of controls
summary(psm_match) # call regression output

# extract matched data
psm_matched_data <- match.data(psm_match) #6,767
psm_matched_data$Index <- 1:nrow(psm_matched_data)

# merge data and create factors
psmPanel <- merge(ELC, psm_matched_data, by = 'AddressIndex') #1,170,765
psmPanel <- psmPanel %>% 
  mutate(NormConsumption = Consumption/size,
         Treatment = ifelse(InitialPeriod<=Period,1,0)) # Treatment variable encompass both the fact of being treated and the time of treatment - turns out 1 when the treated unit gets treatment
psmPanel <- psmPanel %>% mutate(Treatment = ifelse(is.na(Treatment),0,Treatment))
table(psmPanel$Treatment)
psmPanel <- psmPanel %>% filter(NormConsumption <= quantile(NormConsumption, c(0.9999), na.rm = TRUE)) 

# by group
psm_reg <- plm(log(NormConsumption) ~ Treatment + yearELC + monthELC + CoolingDays + HeatingDays, data = psmPanel, model = 'within', index = c('ID','Period')) # no weights included as the results are almost identical (-0.0565 vs.-0.0562) but no clustering allowed for a weighted panel regression with fixed effects 
summary(psm_reg) 

## visualization for bias reduction in standardized percent bias
psm_mean_treated_bef <- summary(psm_match, data = PropertyStats)$sum.all[2:14,1]

psm_mean_control_bef <- summary(psm_match, data = PropertyStats)$sum.all[2:14,2]

psm_avg_var_bef <- c(sqrt((var(PropertyStats$BaselineConsumption[PropertyStats$Group==1]) + var(PropertyStats$BaselineConsumption[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$size[PropertyStats$Group==1]) + var(PropertyStats$size[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$beds[PropertyStats$Group==1]) + var(PropertyStats$beds[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$baths[PropertyStats$Group==1]) + var(PropertyStats$baths[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$PropertyAge[PropertyStats$Group==1]) + var(PropertyStats$PropertyAge[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$market[PropertyStats$Group==1]) + var(PropertyStats$market[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$assessment[PropertyStats$Group==1]) + var(PropertyStats$assessment[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$MedianIncome[PropertyStats$Group==1]) + var(PropertyStats$MedianIncome[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$PovertyBelow[PropertyStats$Group==1]) + var(PropertyStats$PovertyBelow[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$FemaleHouseholder[PropertyStats$Group==1]) + var(PropertyStats$FemaleHouseholder[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$Black[PropertyStats$Group==1]) + var(PropertyStats$Black[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$RentAsIncome35[PropertyStats$Group==1]) + var(PropertyStats$RentAsIncome35[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$SNAP[PropertyStats$Group==1]) + var(PropertyStats$SNAP[PropertyStats$Group==0]))/2))

psm_std_mean_dif_bef <- 100*(psm_mean_treated_bef-psm_mean_control_bef)/psm_avg_var_bef


psm_mean_treated_aft <- summary(psm_match, data = PropertyStats)$sum.matched[2:14,1]
psm_mean_control_aft <- summary(psm_match, data = PropertyStats)$sum.matched[2:14,2]

psm_avg_var_aft <- c(sqrt((var(psm_matched_data$BaselineConsumption[psm_matched_data$Group==1]) + var(psm_matched_data$BaselineConsumption[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$size[psm_matched_data$Group==1]) + var(psm_matched_data$size[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$beds[psm_matched_data$Group==1]) + var(psm_matched_data$beds[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$baths[psm_matched_data$Group==1]) + var(psm_matched_data$baths[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$PropertyAge[psm_matched_data$Group==1]) + var(psm_matched_data$PropertyAge[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$market[psm_matched_data$Group==1]) + var(psm_matched_data$market[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$assessment[psm_matched_data$Group==1]) + var(psm_matched_data$assessment[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$MedianIncome[psm_matched_data$Group==1]) + var(psm_matched_data$MedianIncome[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$PovertyBelow[psm_matched_data$Group==1]) + var(psm_matched_data$PovertyBelow[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$FemaleHouseholder[psm_matched_data$Group==1]) + var(psm_matched_data$FemaleHouseholder[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$Black[psm_matched_data$Group==1]) + var(psm_matched_data$Black[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$RentAsIncome35[psm_matched_data$Group==1]) + var(psm_matched_data$RentAsIncome35[psm_matched_data$Group==0]))/2),
                     sqrt((var(psm_matched_data$SNAP[psm_matched_data$Group==1]) + var(psm_matched_data$SNAP[psm_matched_data$Group==0]))/2))

psm_std_mean_dif_aft <- 100*(psm_mean_treated_aft-psm_mean_control_aft)/psm_avg_var_aft

psm_reduction <- round(summary(psm_match, data = PropertyStats)$reduction[2:14,1],2)                             

names <- c("Average Baseline Consumption", "Property Size", "No. Beds", "No. Baths", "Property Age", "Market Property Value", "Assessment Property Value")

Covariates <- factor(names, ordered = TRUE, levels = rev(names))

psm_bias_df <- data.frame('Covariates' = Covariates,
                          'Before' = psm_std_mean_dif_bef,
                          'After' = psm_std_mean_dif_aft, 
                          'Reduction' = psm_reduction)

gg_psm <- ggplot(psm_bias_df, aes(x=After, xend=Before, y=Covariates, group=Covariates)) +
  geom_vline(colour='#787878',xintercept = 0) +
  geom_dumbbell(color='#8ab7db', 
                size=2, 
                colour_x = '#0e668b',
                colour_xend ='#d2e8fa',
                dot_guide=TRUE, dot_guide_size=0.25) + 
  labs(x="", y="", title="") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.background=element_rect(fill="#ffffff"),
        panel.background=element_rect(fill="#ffffff"),
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(linetype=2),
        panel.grid.major.x=element_blank(),
        axis.ticks=element_blank(),
        legend.position="top",
        panel.border=element_blank())
plot(gg_psm)


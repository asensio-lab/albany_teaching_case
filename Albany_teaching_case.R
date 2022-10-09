### libraries needed

if (!require(c("tidyverse", "ggalt","gdata","gridExtra","rgenoud","MatchIt","plm","Hmisc")))
  install.packages(c("tidyverse","ggalt","gdata","gridExtra","rgenoud","MatchIt","plm","Hmisc"))

library(tidyverse)
library(ggalt)
library(gdata)
library(gridExtra)
library(rgenoud)
library(MatchIt)
library(plm)
library(Hmisc)

### data import

####### TODO: import energy consumption dataset #######
ELC <- read.csv(...) #5,925,488 - residential electric usage 2004-2019
ELC <- ELC %>% 
  ####### TODO: include factor variables for months and years #######
  mutate(yearELC = ..., #factor variable
         monthELC = ...) #factor variable
####### TODO: import property characteristics dataset #######
PropertyStats <- read.csv(...) #16,680 - property data

### fixed effects model

####### TODO: merge datasets by "AddressIndex" #######
ELC_property <- merge(...) #2,931,700 - merge two datasets
# Treatment variable encompass both the fact of being treated and the time of treatment => 'Treatment' equals 1 when the treated unit gets treatment
ELC_property <- ELC_property %>%
  mutate(NormConsumption = Consumption/size, #standardize consumption
         Treatment = ifelse(InitialPeriod<=Period,1,0)) #create treatment variable
ELC_property <- ELC_property %>% mutate(Treatment = ifelse(is.na(Treatment),0,Treatment))
ELC_property <- ELC_property %>% filter(NormConsumption <= quantile(NormConsumption, c(0.9999), na.rm = TRUE)) #remove outliers

# reference model - fixed effects to remove unobserved heterogeneity between properties and time periods
fe_reg <- plm(log(NormConsumption) ~ Treatment + yearELC + monthELC + CoolingDays + HeatingDays, data = ELC_property, model='within', index = c('ID','Period')) 
# Treatment - a dummy variable indicating program participation, 
# yearELC & monthELC - time effects,
# CoolingDays & HeatingDays - monthly number of either days retrieved from National Oceanic and Atmospheric Administration (NOAA), 
# model 'within' - fixed effects model
summary(fe_reg)  #call regression output

### PSM - Propensity Score Matching

####### TODO: add covariates of your choice #######
psm_match <- matchit(Group ~ ...your covariates...,
                     #method='nearest', data=PropertyStats, replace = TRUE, ratio=21)
# method 'nearest' indicates propensity score matching method (nearest neighbor matching), 
# ratio - number of control matches per treatment unit,
# replace = TRUE - multiple time usage of controls
summary(psm_match) # call matching output
                     
# extract matched data
psm_matched_data <- match.data(psm_match) #6,767
psm_matched_data$Index <- 1:nrow(psm_matched_data)
                                      
### visualization for bias reduction in standardized percent bias    

# get covariate means of treated and control properties before matching (N indicates the number of variables chosen)
psm_mean_treated_bef <- summary(psm_match, data = PropertyStats)$sum.all[2:N,1]
psm_mean_control_bef <- summary(psm_match, data = PropertyStats)$sum.all[2:N,2]

# calculate sum of covariate variances of treated and control properties before matching (include VARS of your choice)
psm_avg_var_bef <- c(sqrt((var(PropertyStats$VAR1[PropertyStats$Group==1]) + var(PropertyStats$VAR1[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$VAR2[PropertyStats$Group==1]) + var(PropertyStats$VAR2[PropertyStats$Group==0]))/2),
                     sqrt((var(PropertyStats$VAR3[PropertyStats$Group==1]) + var(PropertyStats$VAR3[PropertyStats$Group==0]))/2),
                     ...)

# calculate standardized percent bias before matching
psm_std_mean_dif_bef <- 100*(psm_mean_treated_bef-psm_mean_control_bef)/psm_avg_var_bef

####### TODO: calculate covariate means of treated and control properties after matching #######
psm_mean_treated_aft <- ...
psm_mean_control_aft <- ...

####### TODO: calculate sum of covariate variances of treated and control properties after matching #######
psm_avg_var_aft <- ...

####### TODO: calculate standardized percent bias after matching #######
psm_std_mean_dif_aft <- ...                  

# legend and date compilation
names <- c("...your covariates in order of appearance in the formula...")
Covariates <- factor(names, ordered = TRUE, levels = rev(names))
psm_bias_df <- data.frame('Covariates' = Covariates,
                          'Before' = psm_std_mean_dif_bef,
                          'After' = psm_std_mean_dif_aft)

# plot bias reduction
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
                     
### regression adjustments after matching

####### TODO: merge energy consumption and post-PSM datasets and create factors #######
psmPanel <- merge(...) #1,170,765 - merge two datasets
psmPanel <- psmPanel %>% 
  mutate(NormConsumption = ..., #standardize consumption
         Treatment = ...) #create treatment variable
psmPanel <- psmPanel %>% mutate(Treatment = ifelse(is.na(Treatment),0,Treatment))
psmPanel <- psmPanel %>% filter(NormConsumption <= quantile(NormConsumption, c(0.9999), na.rm = TRUE)) #remove outliers

####### TODO: add fixed-effects formula  #######
psm_reg <- plm(..., data = psmPanel, model='within', weights=weights, index = c('ID','Period'))
summary(psm_reg)  #call regression output

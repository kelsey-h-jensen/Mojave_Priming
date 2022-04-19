
##
library(tidyverse)

## Cumulative CO2 cals

resp_long <- read.csv("../Mojave_Priming/data/LATR_respiration_long.csv", header = TRUE)
resp_long$ring <- as.factor(resp_long$ring)
resp_long$time_fct <- as.factor(resp_long$time_fct)
resp_long$rep <- as.factor(resp_long$rep)

resp_long_calcs <- resp_long %>% 
  unite(ID, c("amendment", "ring", "rep"), remove = "FALSE") %>% # create unique identifier
  unite(tmt_ID, c("amendment", "treatment"), remove = "FALSE") %>%
  mutate(soil_ugC = soil_mgC*1000) %>% 
  mutate(abs_ppm = CO2_ppm - background_ppm) %>% # calculate ppm CO2 relative to background
  mutate(abs_ppm = replace(abs_ppm, abs_ppm < 0, 0)) %>% # make all negative rates 0
  mutate(CO2_rate = abs_ppm/hours) %>% # calculate rate by dividing by incubation time
  mutate(CO2_rate = replace(CO2_rate, CO2_rate == "NaN", 0)) %>% # make all negative rates 0
  mutate(ug_co2_c = (abs_ppm*0.260)) %>%   # Ideal gas law to calculate ug C in CO2
  mutate(ug_CO2_c_hr = ug_co2_c/hours) %>% 
  mutate(mineralization = ug_co2_c / soil_mgC,
         mineralizability = ((mineralization / (soil_mgC/soil_g)))) %>% 
  filter(!grepl("6_1", ID)) %>%  # Sample 6_1 has no respiration
  filter(!grepl("1_2", ID)) %>% # Very high control backgrounds
  filter(!grepl("5_2", ID)) # Very high control backgrounds

save(resp_long_calcs, file = "../Mojave_Priming/data/resp_long_calcs_outrm.RData")

resp_long_cumu <- resp_long_calcs %>% 
  group_by(rep, ring, amendment) %>% arrange(time_num) %>% 
  mutate(cumu_ppm = cumsum(abs_ppm)) %>% 
  mutate(cumu_ug_co2_c = cumsum(ug_co2_c)) %>% 
  mutate(cumu_mineralization = cumsum(mineralization),
         cumu_mineralizability = cumsum(mineralizability))
# generates cumulative resp over time (ie retains each time point and adds the subsequent ppm)

save(resp_long_cumu, file = "../Mojave_Priming/data/resp_long_cumu.RData")

### Mean data
mean_data <- resp_long_cumu %>% 
  group_by(treatment, amendment) %>%        
  summarize(mean_cumu_ug_co2_c = mean(cumu_ug_co2_c),
            mean_cumu_mineralization = mean(cumu_mineralization),
            mean_cumu_mineralizability = mean(cumu_mineralizability),
            sd_cumu_ug_co2_c = sd(cumu_ug_co2_c),
            se_cumu_ug_co2_c = sd_cumu_ug_co2_c/sqrt(3),
            sd_cumu_mineralization = sd(cumu_mineralization),
            se_cumu_mineralization = sd_cumu_mineralization/sqrt(3),
            sd_cumu_mineralizability = sd(cumu_mineralizability),
            se_cumu_mineralizability = sd_cumu_mineralizability/sqrt(3))

save(mean_data, file= "../Mojave_Priming/data/mean_cumu_resp.RData")


### Priming data

priming <- read.csv("../Mojave_Priming/data/priming_data.csv")
priming$time_fct <- as.factor(priming$time_fct)

priming <- priming %>% unite(ID, c("ring", "rep"), remove = "FALSE") %>% 
  filter(!grepl("6_1", ID)) %>%  # Sample 6_1 has no respiration
  filter(!grepl("1_2", ID)) %>% # Very high control backgrounds
  filter(!grepl("5_2", ID)) # Very high control backgrounds

save(priming, file = "../Mojave_Priming/data/priming.RData")


# cumulative priming

# Need to deicde which option to save! (w or w/o negative priming)
load(file = "../Mojave_Priming/data/priming.RData")
cumu_priming <- priming %>% 
  group_by(ID, time_num) %>% 
  mutate(cumu_ctrl_CO2_ugC = cumsum(ctrl_soil_CO2_ugC)) %>%
  mutate(cumu_soil_CO2_ugC = cumsum(glu_CO2_ugC)) %>% 
  mutate(cumu_glu_CO2_ugC = cumsum(glu_CO2_ugC)) %>% 
  mutate(cumu_primed_ugC = cumsum(primed_ugC)) # Should we add only positive priming? Could change all negative values to 0

save(cumu_priming, file= "../Mojave_Priming/data/cumu_priming.RData")

# Cumulative priming with no negative priming
cumu_priming <- priming %>% select(!contains("hr")) %>% 
  group_by(rep, ring) %>% arrange(time_num) %>% 
  mutate(primed_ugC = replace(primed_ugC, primed_ugC < 0, 0)) %>% 
  mutate(cumu_ctrl_CO2_ugC = cumsum(ctrl_soil_CO2_ugC)) %>%
  mutate(cumu_soil_CO2_ugC = cumsum(glu_soil_ugC)) %>% 
  mutate(cumu_glu_CO2_ugC = cumsum(glu_CO2_ugC)) %>% 
  mutate(cumu_primed_ugC = cumsum(primed_ugC))

save(cumu_priming, file= "../Mojave_Priming/data/cumu_priming.RData")



mean_data2 <- cumu_priming %>% 
  group_by(treatment, time_fct) %>%        
  summarize(mean_cumu_MAOM_CO2_ugC = mean(cumu_MAOM_CO2_ugC),
            mean_cumu_glu_CO2_ugC = mean(cumu_glu_CO2_ugC),
            mean_cumu_primed_ugC = mean(cumu_primed_ugC),
            sd_cumu_MAOM_CO2_ugC = sd(cumu_MAOM_CO2_ugC),
            se_cumu_MAOM_CO2_ugC = sd_cumu_MAOM_CO2_ugC/sqrt(9), # need to adjust this to "count" since I removed some data. Look at line plots from fractions
            sd_cumu_glu_CO2_ugC = sd(cumu_glu_CO2_ugC),
            se_cumu_glu_CO2_ugC = sd_cumu_glu_CO2_ugC/sqrt(9),
            sd_cumu_primed_ugC = sd(cumu_primed_ugC),
            se_cumu_primed_ugC = sd_cumu_primed_ugC/sqrt(9))

save(mean_data2, file= "../Mojave_Priming/data/mean_cumu_priming_time.RData")


mean_data3 <- cumu_priming %>% filter(time_num == "48") %>% 
  group_by(treatment) %>%        
  summarize(mean_cumu_MAOM_CO2_ugC = mean(cumu_MAOM_CO2_ugC),
            mean_cumu_glu_CO2_ugC = mean(cumu_glu_CO2_ugC),
            mean_cumu_primed_ugC = mean(cumu_primed_ugC),
            sd_cumu_MAOM_CO2_ugC = sd(cumu_MAOM_CO2_ugC),
            se_cumu_MAOM_CO2_ugC = sd_cumu_MAOM_CO2_ugC/sqrt(9), # need to adjust this to "count" since I removed some data. Look at line plots from fractions
            sd_cumu_glu_CO2_ugC = sd(cumu_glu_CO2_ugC),
            se_cumu_glu_CO2_ugC = sd_cumu_glu_CO2_ugC/sqrt(9),
            sd_cumu_primed_ugC = sd(cumu_primed_ugC),
            se_cumu_primed_ugC = sd_cumu_primed_ugC/sqrt(9))

save(mean_data3, file= "../Mojave_Priming/data/mean_cumu_primed.RData")

### Priming long
# 1_2, 5_2, and 6_1 were removed
priming_long <- read.csv(file= "../Mojave_Priming/data/priming_data_long.csv", header= TRUE)

cumu_priming_long <- priming_long %>% select(-contains("hr")) %>% 
  mutate(primed_ugC = replace(primed_ugC, primed_ugC < 0, 0)) %>% # probably need a different solution to this eventually
  arrange(rep, ring, amendment, time_num) %>%
  group_by(rep, ring, amendment) %>% 
  mutate(cumu_MAOM_CO2_ugC = cumsum(MAOM_CO2_ugC)) %>%
  mutate(cumu_glu_CO2_ugC = cumsum(glu_CO2_ugC)) %>% 
  mutate(cumu_primed_ugC = cumsum(primed_ugC))

save(cumu_priming_long, file= "../Mojave_Priming/data/cumu_primed_long.RData")


#### Barplot data
barplot_data <- cumu_priming_long %>% 
  group_by(amendment, treatment) %>%        
  summarize(mean_cumu_MAOM_CO2_ugC = mean(cumu_MAOM_CO2_ugC),
            mean_cumu_glu_CO2_ugC = mean(cumu_glu_CO2_ugC),
            mean_cumu_primed_ugC = mean(cumu_primed_ugC),
            sd_cumu_MAOM_CO2_ugC = sd(cumu_MAOM_CO2_ugC),
            se_cumu_MAOM_CO2_ugC = sd_cumu_MAOM_CO2_ugC/sqrt(length(cumu_MAOM_CO2_ugC)), 
            sd_cumu_glu_CO2_ugC = sd(cumu_glu_CO2_ugC),
            se_cumu_glu_CO2_ugC = sd_cumu_glu_CO2_ugC/sqrt(length(cumu_glu_CO2_ugC)),
            sd_cumu_primed_ugC = sd(cumu_primed_ugC),
            se_cumu_primed_ugC = sd_cumu_primed_ugC/sqrt(length(cumu_primed_ugC)))

write.csv(barplot_data, file= "../Mojave_Priming/data/barplot_priming.csv")
# formatted into *extra* long in excel (C_source, mean, sd, se)

barplot_data48 <- cumu_priming_long %>% 
  filter(time_fct == "48")  %>% 
  group_by(amendment, treatment) %>%        
  summarize(mean_cumu_MAOM_CO2_ugC = mean(cumu_MAOM_CO2_ugC),
            mean_cumu_glu_CO2_ugC = mean(cumu_glu_CO2_ugC),
            mean_cumu_primed_ugC = mean(cumu_primed_ugC),
            sd_cumu_MAOM_CO2_ugC = sd(cumu_MAOM_CO2_ugC),
            se_cumu_MAOM_CO2_ugC = sd_cumu_MAOM_CO2_ugC/sqrt(length(cumu_MAOM_CO2_ugC)), 
            sd_cumu_glu_CO2_ugC = sd(cumu_glu_CO2_ugC),
            se_cumu_glu_CO2_ugC = sd_cumu_glu_CO2_ugC/sqrt(length(cumu_glu_CO2_ugC)),
            sd_cumu_primed_ugC = sd(cumu_primed_ugC),
            se_cumu_primed_ugC = sd_cumu_primed_ugC/sqrt(length(cumu_primed_ugC)))

write.csv(barplot_data48, file= "../Mojave_Priming/data/barplot_priming48.csv")


#### Area plot
area_plot <- cumu_priming_long %>% 
  group_by(amendment, treatment, time_fct) %>% 
  summarize(mean_cumu_MAOM_CO2_ugC = mean(cumu_MAOM_CO2_ugC),
            mean_cumu_glu_CO2_ugC = mean(cumu_glu_CO2_ugC),
            mean_cumu_primed_ugC = mean(cumu_primed_ugC),
            sd_cumu_MAOM_CO2_ugC = sd(cumu_MAOM_CO2_ugC),
            se_cumu_MAOM_CO2_ugC = sd_cumu_MAOM_CO2_ugC/sqrt(length(cumu_MAOM_CO2_ugC)), 
            sd_cumu_glu_CO2_ugC = sd(cumu_glu_CO2_ugC),
            se_cumu_glu_CO2_ugC = sd_cumu_glu_CO2_ugC/sqrt(length(cumu_glu_CO2_ugC)),
            sd_cumu_primed_ugC = sd(cumu_primed_ugC),
            se_cumu_primed_ugC = sd_cumu_primed_ugC/sqrt(length(cumu_primed_ugC)))
write.csv(area_plot, file = "../Mojave_Priming/data/areaplot_calcs.csv")  

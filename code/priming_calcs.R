# Data calculations for incubation experiment

##
library(tidyverse)


resp <- read.csv(file = "../Mojave_Priming/data/alldata_raw.csv", header = TRUE)

resp$ring <- as.factor(resp$ring)
resp$treatment <- as.factor(resp$treatment)
resp$time_fct <- as.factor(resp$time_fct)
resp$rep <- as.factor(resp$rep)

# priming calcs
priming_calcs <- resp %>% 
  unite(ID, c("plant", "ring", "rep"), remove = "FALSE") %>%
  mutate(ctrl_soil.CO2 = ctrl_raw.CO2 - ctrl_bkg.CO2) %>% 
  mutate(ctrl_net.CO2 = ctrl_raw.CO2 - ctrl_bkg.CO2) %>% # Same as soil CO2
  mutate(ctrl_soil.CO2 = replace(ctrl_soil.CO2, ctrl_soil.CO2 < 0, 0)) %>% # make all negative soil resp 0
  mutate(ctrl_net.CO2 = replace(ctrl_net.CO2, ctrl_net.CO2 < 0, 0)) %>% # make all negatives net resp 0
  mutate(ctrl_soil.ugC = ctrl_soil.CO2*0.230) %>%
  mutate(ctrl_net.ugC = ctrl_soil.CO2*0.230) %>% # same as soil derived but name is equivalent to glu_net
  mutate(ctrl_soil.ugC.hr = ctrl_soil.ugC / ctrl_inc.hours) %>% 
  mutate(ctrl_net.ugC.hr = ctrl_net.ugC / ctrl_inc.hours) %>% # same as ctrl_soil.ugC.hr
  mutate(ctrl_soil.ugC.hr = replace(ctrl_soil.ugC.hr, ctrl_soil.ugC.hr == "NaN", 0)) %>% # make NAs 0
  mutate(ctrl_net.ugC.hr = replace(ctrl_net.ugC.hr, ctrl_net.ugC.hr == "NaN", 0)) %>% 
  #Calculate d13C of respired C from soil
  mutate(ctrl_frac.bkg = ctrl_bkg.CO2/ctrl_raw.CO2) %>% 
  mutate(ctrl_frac.soil = 1-ctrl_frac.bkg) %>% 
  mutate(ctrl_soil.isoC = (ctrl_raw.isoC-(ctrl_bkg.isoC*ctrl_frac.bkg))/ctrl_frac.soil) %>%
  mutate(glu_net.CO2 = glu_raw.CO2 - glu_bkg.CO2) %>% 
  mutate(glu_net.CO2 = replace(glu_net.CO2, glu_net.CO2 < 0, 0)) %>% # make all negative rates 0
  mutate(glu_net.ugC = glu_net.CO2*0.230) %>% 
  mutate(glu_net.ugC.hr = glu_net.ugC / glu_inc.hours) %>% 
  mutate(glu_net.ugC.hr = replace(glu_net.ugC.hr, glu_net.ugC.hr == "NaN", 0)) %>% # make NAs 0
  mutate(glu_frac.bkg = glu_bkg.CO2/glu_raw.CO2) %>% 
  # fraction of CO2 from soil + amendment
  mutate(glu_frac.soil.amend = 1-glu_frac.bkg) %>% 
  # d13C from soil + amendment
  mutate(glu_soil.amend.isoC = -((glu_raw.isoC-(glu_bkg.isoC*glu_frac.bkg))/(1-glu_frac.bkg))) %>% 
  # fraction CO2 from soil only
  mutate(glu_frac.amendment = -(glu_soil.amend.isoC-soc_isoC)/(glu_amendment.isoC-soc_isoC)) %>% 
  # fraction CO2 from amendment only
  mutate(glu_frac.soil = 1- glu_frac.amendment) %>% 
  mutate(glu_soil.ugC = glu_frac.soil*glu_net.ugC) %>% 
  mutate(glu_soil.ugC = replace(glu_soil.ugC, glu_soil.ugC < 0, 0)) %>% # make all ugC 0
  mutate(glu_soil.ugC.hr = glu_soil.ugC/glu_inc.hours) %>% 
  mutate(glu_amendment.ugC = glu_net.ugC*glu_frac.amendment) %>%
  mutate(glu_amendment.ugC.hr = glu_amendment.ugC/glu_inc.hours) %>%
  mutate(glu_primed.ugC = glu_soil.ugC - ctrl_soil.ugC ) %>%
  mutate(glu_rel.primed.ugC = glu_primed.ugC / doc_ugC) %>%
  mutate(glu_rel.primed.ugC.hr = glu_primed.ugC / glu_inc.hours) %>%
  mutate(glu_primed.ugC.hr = glu_primed.ugC / glu_inc.hours ) %>% 
  mutate(glu_primed.isoC = (glu_soil.amend.isoC-(glu_amendment.isoC*glu_frac.amendment))/glu_frac.soil) %>% 
  mutate(glu_priming.effect = glu_primed.ugC.hr / ctrl_soil.ugC.hr) %>%
  mutate(ctrl_mineralization.C = ctrl_soil.ugC / soc_mgC) %>% 
  mutate(glu_mineralization.C = glu_soil.ugC / soc_mgC) %>% 
  mutate(ctrl_mineralization.hr = ctrl_mineralization.C / ctrl_inc.hours) %>%
  mutate(glu_mineralization.hr = glu_mineralization.C / glu_inc.hours) %>% 
  mutate(ctrl_rel.doc.C = ctrl_soil.ugC / doc_ugC) %>% 
  mutate(glu_rel.doc.C = glu_soil.ugC / doc_ugC) %>% 
  mutate(ctrl_rel.doc.hr = ctrl_rel.doc.C / ctrl_inc.hours) %>%
  mutate(glu_rel.doc.hr = glu_rel.doc.C / glu_inc.hours) %>%
  filter(!grepl("LATR_6_1", ID)) %>%  # Sample 6_1 has no respiration
  filter(!grepl("LATR_1_6", ID)) %>% # negative priming at -1500
  filter(!grepl("LATR_1_2", ID)) %>% # Very high control backgrounds, check these data
  filter(!grepl("LATR_5_2", ID)) %>%  # Very high control backgrounds, check these data
  filter(!grepl("INSP_5_4", ID)) # Zero respiration for controls

save(priming_calcs, file= "../Mojave_priming/data/priming_calcs.RData")


# Cumulative resp

resp_cumu <- priming_calcs %>% 
  mutate(glu_pos.primed.ugC = replace(glu_primed.ugC, glu_primed.ugC < 0, 0)) %>%
  mutate(glu_neg.primed.ugC = replace(glu_primed.ugC, glu_primed.ugC > 0, 0)) %>%
  group_by(plant, ring, rep) %>% arrange(plant, ring, rep, time_num) %>%
  mutate(ctrl_cumu.soil.ugC = cumsum(ctrl_soil.ugC)) %>% 
  mutate(glu_cumu.soil.ugC = cumsum(glu_soil.ugC)) %>% 
  mutate(glu_cumu.amend.ugC = cumsum(glu_amendment.ugC)) %>%
  mutate(glu_cumu.primed.ugC = cumsum(glu_primed.ugC)) %>% 
  mutate(glu_cumu.rel.primed.ugC = cumsum(glu_rel.primed.ugC)) %>% 
  mutate(ctrl_cumu.mineralization.C = cumsum(ctrl_mineralization.C)) %>% 
  mutate(glu_cumu.mineralization.C = cumsum(glu_mineralization.C)) %>% 
  mutate(ctrl_cumu.rel.doc.C = cumsum(ctrl_mineralization.C)) %>% 
  mutate(glu_cumu.rel.doc.C = cumsum(glu_mineralization.C)) %>% 
  mutate(glu_cumu.pos.primed.ugC = cumsum(glu_pos.primed.ugC)) %>%
  mutate(glu_cumu.neg.primed.ugC = cumsum(glu_neg.primed.ugC)) %>%
  mutate(glu_cumu.net.primed.ugC = cumsum(glu_primed.ugC)) %>% 
  mutate(glu_cumu.rel.primed.ugC = cumsum(glu_rel.primed.ugC)) %>% 
  mutate(glu_cumu.priming.effect = glu_cumu.primed.ugC/ ctrl_cumu.soil.ugC) 

# Pivot long, remove raw data
cumu_resp_long <- resp_cumu %>% 
  pivot_longer(contains('.'), 
               names_to = c("amendment", "measurement"),
               values_to = "values", names_sep = '_') %>% 
  filter(!grepl("frac|isoC|raw|bkg", measurement)) %>% 
  pivot_wider(names_from = "measurement", values_from= "values") %>% 
  replace(is.na(.), 0)

save(cumu_resp_long, file= "../Mojave_Priming/data/cumu_resp_long.RData")

### Cumulative positive, negative, and net priming

cumu_priming <- priming_calcs %>% 
  mutate(glu_pos.primed.ugC = replace(glu_primed.ugC, glu_primed.ugC < 0, 0)) %>%
  mutate(glu_neg.primed.ugC = replace(glu_primed.ugC, glu_primed.ugC > 0, 0)) %>%
  group_by(plant, ring, rep) %>% arrange(plant, ring, rep, time_num) %>%
  mutate(glu_cumu.pos.primed.ugC = cumsum(glu_pos.primed.ugC)) %>%
  mutate(glu_cumu.neg.primed.ugC = cumsum(glu_neg.primed.ugC)) %>%
  mutate(glu_cumu.net.primed.ugC = cumsum(glu_primed.ugC)) %>% 
  mutate(glu_cumu.rel.primed.ugC = cumsum(glu_rel.primed.ugC)) %>%
  pivot_longer(contains('.'), 
               names_to = c("amendment", "measurement"),
               values_to = "values", names_sep = '_') %>% 
  filter(grepl("primed", measurement)) %>% 
  filter(!grepl("ctrl", measurement)) %>% 
  filter(!grepl("isoC", measurement)) %>% 
  pivot_wider(names_from = "measurement", values_from= "values") 


save(cumu_priming, file= "../Mojave_Priming/data/cumu_priming.RData")


### Lineplot means ####

mean_data <- cumu_resp_long %>% 
  group_by(plant, amendment, treatment, time_num, time_fct) %>%        
  summarize(mean_soil.ugC.hr = mean(soil.ugC.hr),
            sd_soil.ugC.hr = sd(soil.ugC.hr),
            se_soil.ugC.hr = sd_soil.ugC.hr/sqrt(length(soil.ugC.hr)),
            mean_cumu.soil.ugC = mean(cumu.soil.ugC),
            sd_cumu.soil.ugC = sd(cumu.soil.ugC),
            se_cumu.soil.ugC = sd_cumu.soil.ugC/sqrt(length(cumu.soil.ugC)), 
            mean_cumu.amend.ugC = mean(cumu.amend.ugC),
            sd_cumu.amend.ugC = sd(cumu.amend.ugC),
            se_cumu.amend.ugC = sd_cumu.amend.ugC/sqrt(length(cumu.amend.ugC)),
            mean_cumu.net.primed.ugC = mean(cumu.net.primed.ugC),
            sd_cumu.net.primed.ugC = sd(cumu.net.primed.ugC),
            se_cumu.net.primed.ugC = sd_cumu.net.primed.ugC/sqrt(length(cumu.net.primed.ugC)),
            mean_cumu.pos.primed.ugC = mean(cumu.pos.primed.ugC),
            sd_cumu.pos.primed.ugC = sd(cumu.pos.primed.ugC),
            se_cumu.pos.primed.ugC = sd_cumu.pos.primed.ugC/sqrt(length(cumu.pos.primed.ugC)), 
            mean_cumu.neg.primed.ugC = mean(cumu.neg.primed.ugC),
            sd_cumu.neg.primed.ugC = sd(cumu.neg.primed.ugC),
            se_cumu.neg.primed.ugC = sd_cumu.neg.primed.ugC/sqrt(length(cumu.neg.primed.ugC)),
            mean_net.ugC.hr = mean(net.ugC.hr),
            sd_net.ugC.hr = sd(net.ugC.hr),
            se_net.ugC.hr = sd_net.ugC.hr/sqrt(length(net.ugC.hr)),
            mean_primed.ugC.hr = mean(primed.ugC.hr),
            sd_primed.ugC.hr = sd(primed.ugC.hr),
            se_primed.ugC.hr = sd_primed.ugC.hr/sqrt(length(primed.ugC.hr)),
            mean_mineralization.hr = mean(mineralization.hr),
            sd_mineralization.hr = sd(mineralization.hr),
            se_mineralization.hr = sd_mineralization.hr/sqrt(length(mineralization.hr))) %>% 
  pivot_longer(contains('.'), 
               names_to = c("stat", "Csource"),
               values_to = "values", names_sep = '_') %>% 
  pivot_wider(names_from = "stat", values_from= "values")

save(mean_data, file = "../Mojave_Priming/data/mean_data.Rdata")



#### I think all of the data below this can be deleted. Need to determine if mean mineralization is useful
baseC_data <- abs_cumu_resp_long %>% 
  mutate(cumu.base.ugC = cumu.soil.ugC - cumu.primed.ugC) %>%
  group_by(plant, amendment, treatment, time_fct) %>% 
  pivot_longer(contains('cumu'), 
               names_to = c("Csource"),
               values_to = "values") %>% 
  filter(!grepl("ctrl", amendment)) %>% 
    filter(!grepl("soil", Csource))

save(baseC_data, file = "../Mojave_Priming/data/baseC_data.RData")

mean_baseC_data <- abs_cumu_resp_long %>% 
  filter(!grepl("ctrl", amendment)) %>%
  mutate(cumu.base.ugC = cumu.soil.ugC - cumu.primed.ugC) %>%
  group_by(plant, treatment, time_fct) %>% 
  summarize(mean_cumu.base.ugC = mean(cumu.base.ugC),
          mean_cumu.primed.ugC = mean(cumu.primed.ugC),
          mean_cumu.amend.ugC = mean(cumu.amend.ugC),
          sd_cumu.base.ugC = sd(cumu.base.ugC),
          se_cumu.base.ugC = sd_cumu.base.ugC/sqrt(length(cumu.base.ugC)), 
          sd_cumu.primed.ugC = sd(cumu.primed.ugC),
          se_cumu.primed.ugC = sd_cumu.primed.ugC/sqrt(length(cumu.primed.ugC)), 
          sd_cumu.amend.ugC = sd(cumu.amend.ugC),
          se_cumu.amend.ugC = sd_cumu.amend.ugC/sqrt(length(cumu.amend.ugC))) %>% 
  pivot_longer(contains('.'), 
               names_to = c("stat", "Csource"),
               values_to = "values", names_sep = '_') %>% 
  filter(!grepl("soil", Csource)) %>% 
  pivot_wider(names_from = "stat", values_from= "values")

save(mean_baseC_data, file = "../Mojave_Priming/data/mean_baseC_data.RData")

mean_mineralization <- resp_long %>% 
  group_by(treatment, plant, amendment, time_num) %>% 
  summarize(mean_mineralization = mean(mineralization.C), 
            sd_mineralization = sd(mineralization.C),
            se_mineralization = sd_mineralization/sqrt(length(mineralization.C)),
            mean_mineralization.hr = mean(mineralization.hr),
            sd_mineralization.hr = sd(mineralization.hr),
            se_mineralization.hr = sd_mineralization.hr/sqrt(length(mineralization.hr))) %>% 
  pivot_longer(contains('min'), 
               names_to = c("stat", "response"),
               values_to = "values", names_sep = '_') %>% 
  pivot_wider(names_from = "stat", values_from= "values")


mean_mineralization <- resp_long %>% 
  group_by(treatment, plant, amendment, time_num) %>% 
  summarize(mean_mineralization.C = mean(mineralization.C), 
            sd_mineralization.C = sd(mineralization.C),
            se_mineralization.C = sd_mineralization.C/sqrt(length(mineralization.C)),
            mean_mineralization.hr = mean(mineralization.hr),
            sd_mineralization.hr = sd(mineralization.hr),
            se_mineralization.hr = sd_mineralization.hr/sqrt(length(mineralization.hr)),
            mean_resp.C = mean(soil.ugC),
            sd_resp.C = sd(soil.ugC),
            se_resp.C = sd_resp.C/sqrt(length(soil.ugC)),
            mean_resp.hr = mean(soil.ugC.hr),
            sd_resp.hr = sd(soil.ugC.hr),
            se_resp.hr = sd_resp.hr/sqrt(length(soil.ugC.hr))) %>% 
  pivot_longer(contains('.'), 
               names_to = c("stat", "response"),
               values_to = "values", names_sep = '_') %>% 
  pivot_wider(names_from = "stat", values_from= "values")

save(mean_mineralization, file= "../Mojave_Priming/data/mean_mineralization.RData")




### Lineplot data

# Group by plant, treatment, amendment, time
# calc means and se
# Same as barplot data?

lineplot_rate <- abs_cumu_resp_long %>% 
  group_by(plant, amendment, treatment, time_fct) %>%        
  summarize(mean_soil.ugC = mean(soil.ugC),
            mean_primed.ugC = mean(primed.ugC),
            sd_soil.ugC = sd(soil.ugC),
            se_soil.ugC = sd_soil.ugC/sqrt(length(soil.ugC)),
            sd_primed.ugC = sd(primed.ugC),
            se_primed.ugC = sd_primed.ugC/sqrt(length(primed.ugC))) %>% 
  pivot_longer(contains('.'), 
               names_to = c("stat", "Csource"),
               values_to = "values", names_sep = '_') %>% 
  pivot_wider(names_from = "stat", values_from= "values")

save(barplot_data, file = "../Mojave_Priming/data/lineplot_data.Rdata")

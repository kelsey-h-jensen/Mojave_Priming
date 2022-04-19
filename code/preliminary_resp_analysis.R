library(tidyverse)

respiration <- read.csv("20211028_LATR_respiration.csv", header = TRUE)
respiration$ring <- as.factor(respiration$ring)
respiration$time_fct <- as.factor(respiration$time_fct)
respiration$rep <- as.factor(respiration$rep)

# Cleaned up data including adjusted baselines for the control samples

resp <- respiration %>% 
  mutate(control_abs = control_CO2 - control_bkg_CO2) %>% 
  mutate(glucose_abs = glucose_CO2 - glucose_bkg_CO2) %>% 
  mutate(control_abs = replace(control_abs, control_abs < 0, 0)) %>% 
  mutate(glucose_abs = replace(glucose_abs, glucose_abs < 0, 0)) %>% 
  mutate(control_rate = control_abs/control_hours) %>% 
  mutate(glucose_rate = glucose_abs/glucose_hours) %>% 
  mutate(control_rate = replace(control_rate, control_rate == "NaN", 0)) %>% 
  mutate(glucose_rate = replace(glucose_rate, glucose_rate == "NaN", 0)) %>%
  mutate(net_resp = glucose_rate - control_rate)

resp %>% group_by(treatment, time_num) %>% 
  summarise_if(is.numeric, mean) %>% 
  unite(ID, c("treatment"), remove = "FALSE") %>% 
  ggplot(aes(x= time_num, y = net_resp, group = ID, color= ID)) +
  geom_point() + geom_path( aes(x = time_num, y = net_resp, group = ID, color= ID))

resp  %>% 
  ggplot( aes(x = time, y = control_rate, group= treatment, color= treatment)) + 
  geom_point() + geom_path(data = resp, aes(x = time, y = control_rate, group= treatment, color= treatment))

resp  %>% 
  ggplot( aes(x = time, y = glucose_rate, group= treatment, color= treatment)) + 
  geom_point() + geom_path(data = resp, aes(x = time, y = glucose_rate, group= treatment, color= treatment))


resp_long <- read.csv("20211028_LATR_respiration_long.csv", header = TRUE)
resp_long$ring <- as.factor(resp_long$ring)
resp_long$time_fct <- as.factor(resp_long$time_fct)
resp_long$rep <- as.factor(resp_long$rep)

resp_long <- resp_long %>% 
  mutate(abs_ppm = CO2_ppm - background_ppm) %>% 
  mutate(CO2_rate = abs_ppm/hours) %>% 
  mutate(abs_ppm = replace(abs_ppm, abs_ppm < 0, 0)) %>% 
  mutate(CO2_rate = replace(CO2_rate, CO2_rate == "NaN", 0))


resp_means <- resp_long %>% 
  group_by(amendment, treatment, time_fct) %>% 
  summarise_if(is.numeric, mean) %>% 
  unite(ID, c("amendment", "treatment"), remove = "FALSE") 

# Rates over time
resp_means %>% 
  ggplot( aes(x = time_fct, y = CO2_rate, group = ID, color= amendment)) +
  geom_point(aes(shape = treatment), size = 3) + geom_line(aes(linetype= treatment), stat= "identity")





# Meeting with Erika

library(lme4)
library(lmerTest)
library(emmeans)


# Ask Erika about how to transform some of these 0 heavy data


load("../Mojave_Priming/data/cumu_resp_long.RData")

cumu_resp_long %>%  ggplot(aes(x=time_num, y= cumu.soil.ugC, color = treatment, group = ID)) + geom_line()+
  facet_grid(amendment ~plant) + scale_y_continuous(trans = "sqrt")

lm1 <- lmer(sqrt(cumu.soil.ugC.hr) ~ (plant * treatment * amendment * time_fct) + 
              (1|ring) + (1|ID) + (1|ID:amendment), data = cumu_resp_long) 

# Degrees of fredom for LATR are much lower than INSP

plot(predict(lm1), residuals(lm1)) 
qqnorm(residuals(lm1)); qqline(residuals(lm1)) 
hist(residuals(lm1)) 

summary(lm1)
anova(lm1)

paircomp <- summary(emmeans(lm1, pairwise ~ treatment|plant|amendment|time_fct)$contrasts)


### Cumulative resp ###

load("../Mojave_Priming/data/abs_cumu_resp_long.RData")

cumu.base <- abs_cumu_resp_long %>% 
  filter(amendment == "ctrl") 


lm2 <- lmer(sqrt(cumu.soil.ugC) ~ (plant * treatment * time_fct) +
             (1|ring), data = cumu.base)
plot(predict(lm2), residuals(lm2)) 
qqnorm(residuals(lm2)); qqline(residuals(lm2))
hist(residuals(lm2))

summary(lm2)
anova(lm2)

paircomp.2 <- summary(emmeans(lm2, pairwise ~ treatment|plant|time_fct)$contrasts)


#### DOC ####

cumu.C.4 <- cumu_resp_long %>% filter(time_fct == "4") %>% filter(amendment == "glu")

lm4 <- lm(sqrt(cumu.soil.ugC) ~ (plant + treatment + amendment + doc_ugC + soc_mgC), data = cumu_resp_long)
car::vif(lm4)

lm3 <- lmer(sqrt(cumu.soil.ugC) ~ 
              (plant + treatment + amendment + doc_ugC + soc_mgC)^5 + (1|ring) + (1|ID), data = cumu_resp_long)

plot(predict(lm3), residuals(lm3)) 
qqnorm(residuals(lm3)); qqline(residuals(lm3))
hist(residuals(lm3))

summary(lm3)
anova(lm3)

lm4 <- lmer(sqrt(cumu.soil.ugC) ~ (plant + treatment + amendment)^3 + (1|ring) + (1|ID), data = cumu.C.4)
anova(lm4)

cumu.C.4 %>% ggplot(aes(x= plant, y= cumu.soil.ugC, color = treatment)) +
  geom_boxplot() + facet_grid(~amendment) + scale_y_continuous(trans = "sqrt")


lm5 <- lmer(sqrt(cumu.soil.ugC) ~ (treatment + amendment + doc_ugC)^3 + (1|ring) + (1|ID), data = cumu.C.4)
anova(lm5)

cumu.C.4 %>% 
  ggplot(aes(x= doc_ugC, y= cumu.soil.ugC, color = treatment)) +
  geom_point() + geom_smooth(aes(group= interaction(treatment, amendment)), method = lm)+
  facet_grid(~amendment) + scale_y_continuous(trans = "sqrt")

lm6 <- lmer(sqrt(cumu.soil.ugC) ~ (treatment + amendment + soc_mgC)^3 + (1|ring) + (1|ID), data = cumu.C.4)
anova(lm6)

cumu.C.4 %>% 
  ggplot(aes(x= soc_mgC, y= cumu.soil.ugC, color = treatment)) +
  geom_point() + geom_smooth(aes(group= interaction(treatment, amendment)), method = lm)+
  facet_grid(~amendment) + scale_y_continuous(trans = "sqrt")


lm7 <- lmer(sqrt(cumu.soil.ugC) ~ (treatment + ugC_gSoil)^2 + (1|ring), data = cumu.C.4)
anova(lm7)

resp_cumu %>% filter(time_fct == "168") %>% 
  ggplot(aes(x= whc, y= ctrl_cumu.soil.ugC, color = treatment)) +
  geom_point() + geom_smooth(aes(group= interaction(treatment)), method = lm)+
   scale_y_continuous(trans = "sqrt") +
  stat_cor(method = "spearman")

resp_cumu %>% filter(time_fct == "168") %>% 
  ggplot(aes(x= whc, y= glu_cumu.soil.ugC, color = treatment)) +
  geom_point() + geom_smooth(aes(group= interaction(treatment)), method = lm)+
  scale_y_continuous(trans = "sqrt") +
  stat_cor(method = "spearman")


########

# Figures for Erika 

ggplot(corr.data.ctrl, aes(x = doc_ugC, y = soc_mgC, fill = plant)) + 
  geom_point() + geom_smooth(method = "lm") + 
  stat_cor(method="spearman", label.y = c(4,5), 
           size= 3, r.accuracy = 0.01, p.accuracy =0.01)

ggplot(cumu_resp_long, aes(x= doc_ugC, y = cumu.soil.ugC, fill= treatment, color = treatment)) + 
  geom_point(aes(shape = plant)) + geom_smooth(method = "lm", aes(linetype = plant), se= FALSE) + 
  facet_grid(amendment~time_fct)

cumu_resp_long %>% filter(time_fct == "4") %>% 
  ggplot(aes(x = doc_ugC, y = cumu.soil.ugC, fill = plant)) + 
  geom_point() + geom_smooth(method = "lm") + facet_grid(~treatment)+
  stat_cor(method="spearman",
           size= 3, r.accuracy = 0.01, p.accuracy =0.01)

cumu_resp_long %>% filter(time_fct == "4") %>% filter(plant == "INSP") %>% 
  ggplot(aes(x = doc_ugC, y = cumu.primed.ugC, fill = treatment, color = treatment)) + 
  geom_point() + geom_smooth(method = "lm") + facet_grid(~amendment)+
  stat_cor(method="spearman", 
           size= 3, r.accuracy = 0.01, p.accuracy =0.01)

cumu_resp_long %>% filter(time_fct == "4") %>% filter(plant == "INSP") %>% 
  ggplot(aes(x = soc_mgC, y = cumu.primed.ugC, fill = treatment, color = treatment)) + 
  geom_point() + geom_smooth(method = "lm") + facet_grid(~amendment)+
  stat_cor(method="spearman", 
           size= 3, r.accuracy = 0.01, p.accuracy =0.01)


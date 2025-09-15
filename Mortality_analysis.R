#morality analysis 
#OR table
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(quotidieR)
library(survival)
library(survminer)
library(gt)
library(gtsummary)


asv <- read.csv("ndi_asv.rds")
asv_sick <- asv %>% ps_filter(month == 0)%>% 
  ps_filter(sick == TRUE)
  
plot_data <- asv_sick %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Genus") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()

death_bact <- plot_data %>% 
  select(death, "Escherichia","Shigella","Enterobacter","Bacteroides","Klebsiella","Enterococcus","Parabacteroides",
         "Haemophilus","Akkermansia") %>%
  tbl_uvregression(
    method = glm,
    y = death,
    method.args = list(family = binomial),
    exponentiate = TRUE,
    pvalue_fun = ~ style_pvalue(.x, digits = 2)
  ) %>%
  add_global_p() %>% # add global p-value
  add_nevent() %>% # add number of events of the outcome
  add_q() 
gtsave(as_gt(death_bact), "overall_death_bacteria_genus.pdf")
death_bact <- death_bact %>% as.data.frame()

death_bact <- death_bact %>% separate("**95% CI**", c('llim', 'ulim'), sep=",")

fpp <- death_bact %>% mutate(bug = fct_reorder(`**Characteristic**`, `**p-value**`, .desc = TRUE)) %>%
  ggplot(aes(x=bug, y= as.numeric(`**OR**`), ymin=as.numeric(llim), ymax=as.numeric(ulim))) +
  geom_pointrange(aes(color = if_else(`**p-value**` < 0.05, TRUE, FALSE))) +
  geom_hline(yintercept=1, lty=2) +
  geom_text(aes(label = paste0("p=", round(as.numeric(`**p-value**`), 2))), 
            hjust = -0.35, 
            vjust = -0.5, 
            color = "black") +
  coord_flip() +
  labs(y="Odds ratio (95% Confidence Interval)",x="Genus",title="Association with All Mortality") +
  scale_color_calc()+
  theme_Publication() +
  theme(plot.title = element_text(margin = margin(b = 20)))
ggsave("bug_OR.pdf", fpp, width = 6, height = 5, dpi = 600, device = "pdf")

plot_data <- asv_sick %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Family") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()
#now for anemia. with hgb < 7
anemia <- asv_sick %>% ps_filter(hgbbase <= 7)
plot_data <- anemia %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Species") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()
death_df <- plot_data %>% drop_na(sm_death12m) %>%
  mutate(days_to_living = ifelse(is.na(deathdaysto12m), 365, deathdaysto12m),
         doe_date = dmy(doe),
         admission_month = month(doe_date))
death_bact <- plot_data %>% 
  select(death, "Escherichia coli","Bifidobacterium pseudocatenulatum","Streptococcus gallolyticus") %>%
  tbl_uvregression(
    method = glm,
    y = death,
    method.args = list(family = binomial),
    exponentiate = TRUE,
    pvalue_fun = ~ style_pvalue(.x, digits = 2)
  ) %>%
  add_global_p() %>% # add global p-value
  add_nevent() %>% # add number of events of the outcome
  add_q() 
gtsave(as_gt(death_bact), "anemia_pd_death_bacteria_species.pdf")
death_bact <- death_bact %>% as.data.frame()

death_bact <- death_bact %>% separate("**95% CI**", c('llim', 'ulim'), sep=",")

e_cox <- coxph(Surv(days_to_living, sm_death12m) ~ `Escherichia coli`, data = death_df)
e<-summary(e_cox)
e_coef <- e$coefficients
e_cont <- e$conf.int
ee <- cbind(e_coef, e_cont)
b_cox <- coxph(Surv(days_to_living, sm_death12m) ~ `Bifidobacterium pseudocatenulatum`, data = death_df)
b<-summary(b_cox)
b_coef <- b$coefficients
b_cont <- b$conf.int 
bb <- cbind(b_coef, b_cont)
s_cox <- coxph(Surv(days_to_living, sm_death12m) ~ `Streptococcus gallolyticus`, data = death_df)
s<-summary(s_cox)
s_coef <- s$coefficients
s_cont <- s$conf.int
ss<-cbind(s_coef, s_cont)
ebs <- rbind(ee, bb, ss)
ebs <- ebs %>% as_tibble(rownames = NA) %>% rownames_to_column()
ebs_bact<- ebs %>% filter(rowname %in% c("`Escherichia coli`", "`Bifidobacterium pseudocatenulatum`","`Streptococcus gallolyticus`" ))

library(ggthemes)
library(ggpubr)
library(quotidieR)
ebs_bact$rowname <- c("Escherichia coli", "Bifidobacterium pseudocatenulatum","Streptococcus gallolyticus" )
ebs_bact <- ebs_bact %>% mutate(rowname = fct_reorder(rowname, coef)) %>%
  mutate(sig = if_else(`Pr(>|z|)` < 0.05, "yes", "no")) 
fpp <- ebs_bact %>% ggplot(aes(x=rowname, y= `exp(coef)`, ymin=`lower .95`, ymax=`upper .95`)) +
  geom_pointrange(aes(color = sig)) +
  geom_hline(yintercept=1, lty=2) +
  geom_text(aes(label = paste0("p=", round(`Pr(>|z|)`, 2))), 
            hjust = -0.3, 
            vjust = -0.5, 
            color = "black") +
  coord_flip() +
  labs(y="Hazard ratio (95% Confidence Interval)") +
  scale_color_calc()+
  theme_Publication()
ggsave("species_anemia_post_discharge_death.png", fpp, width = 7, height = 2.5, dpi = 1200, device = "png")

#inpatient deaths
univariate_inpat <- new %>% mutate(doe_date = dmy(doe),
                                   admission_month = month(doe_date),
                                   Sex = case_when(sex ==1 ~ "Male",
                                                   sex ==2 ~ "Female"),
                                   Site = case_when(site == 1 ~ "Mulago",
                                                    site == 2 ~ "Jinja"),
                                   hbsgt_bin = case_when(hbsgt %in% c(1,2) ~ "AA/AS",
                                                         hbsgt == 3 ~ "SS"),
                                   age_in_months = round(age * 12),
                                   prior_abx_oral = case_when(amoxyl == 1 ~ 1,
                                                              cipro == 1 ~ 1,
                                                              cotrimoxazole == 1 ~ 1,
                                                              flagyl == 1 ~ 1,
                                                              othantibi == "Ampiclox" ~ 1,
                                                              othantibi == "Azithromycin syrup" ~ 1,
                                                              othantibi == "Cefladriox" ~ 1,
                                                              othantibi == "Cephalexin" ~ 1,
                                                              othantibi == "Doxycline" ~ 1,
                                                              othantibi == "Erythromycin" ~ 1,
                                                              .default = 0),
                                   date_to_stool = (as.Date(stool_dt0, format = "%d/%m/%Y") - as.Date(doe, format="%d/%m/%Y")),
                                   stool_day = case_when(date_to_stool == 0 ~ "day 0",
                                                         date_to_stool == 1 ~ "day 1",
                                                         date_to_stool == 2 ~ "day 2",
                                                         date_to_stool == 3 ~ "day 3",
                                                         date_to_stool == 4 ~ "day 4",
                                                         date_to_stool == 5 ~ "day 5",
                                                         date_to_stool == 6 ~ "day 6",
                                                         date_to_stool == 7 ~ "day 7",
                                                         date_to_stool > 7 ~ "past one week",
                                                         .default = "unknown"),
                                   sick = case_when(group %in% c(1,2,3,4,5) ~ "SM",
                                                    group == 6 ~ "CC"),
                                   sma = case_when(enr_sma == 0 ~ "No",
                                                   enr_sma == 1 ~ "Yes",
                                                   .default = NA),
                                   h_uric = case_when(uric_se1 > 7 ~ "Yes",
                                                      uric_se1 <= 7 ~ "No",
                                                      .default = NA),
                                   aki = case_when(haki == 1 ~ "Yes",
                                                   haki == 0 ~ "No"),
                                   coma = case_when(enr_cm == 0 ~ "No",
                                                    enr_cm == 1 ~ "Yes",
                                                    .default = NA),
                                   tff3 = case_when(tff3_pl1 < 4.08 ~ "No",
                                                    tff3_pl1 >= 4.08 ~ "Yes",
                                                    .default = NA),
                                   acid_breaths = case_when(acidotic == 0 ~ "No",
                                                            acidotic == 1 ~ "Yes",
                                                            .default = NA)
)  %>%
  filter(sick == "SM") %>%
  select(outcome, coma, aki, acid_breaths, tff3, h_uric, sma, Sex, Site, age_in_months, hbsgt_bin, admission_month) %>%
  tbl_uvregression(
    method = glm,
    y = outcome,
    method.args = list(family = binomial),
    exponentiate = TRUE,
    pvalue_fun = ~ style_pvalue(.x, digits = 2)
  ) %>%
  add_global_p() %>% # add global p-value
  add_nevent() %>% # add number of events of the outcome
  add_q() %>% # adjusts global p-values for multiple testing
  bold_p() %>% # bold p-values under a given threshold (default 0.05)
  bold_p(t = 0.10, q = TRUE) %>% # now bold q-values under the threshold of 0.10
  bold_labels()
gtsave(as_gt(univariate_inpat), "univariate_inpatient_death.pdf")

#Kaplan Meier deaths
##kaplan curve 
library(survminer)
death_df <- new %>% filter(! group == 6) %>% drop_na(sm_death12m) %>%
  mutate(days_to_living = ifelse(is.na(deathdaysto12m), 365, deathdaysto12m))

km_fit <- survfit(Surv(days_to_living, sm_death12m) ~ 1, data=death_df)
summary(km_fit, times = c(1,30,60,90*(1:10)))
ggsurvplot(km_trt_fit, conf.int = T)
ggsave("pd_death.png", autoplot(km_fit), width = 6, height = 4, dpi = 1200, device = "png") 
#anemia
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ enr_sma, data=death_df)
ane <- ggsurvplot(km_trt_fit, conf.int = T)
surv_pvalue(km_trt_fit)
survdiff(Surv(days_to_living, sm_death12m) ~ enr_sma, data=death_df)
ggsave("pd_death_enr_sma.png", ane$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#coma
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ enr_cm, data=death_df)
coma<-ggsurvplot(km_trt_fit, conf.int = T)
survdiff(Surv(days_to_living, sm_death12m) ~ enr_cm, data=death_df)
surv_pvalue(km_trt_fit)
ggsave("pd_death_enr_cm.png", coma$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#acidotic
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ acidotic, data=death_df)
acid<-ggsurvplot(km_trt_fit, conf.int = T)
surv_pvalue(km_trt_fit)
survdiff(Surv(days_to_living, sm_death12m) ~ acidotic, data=death_df)
ggsave("pd_death_acidotic.png", acid$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#intestinal injury
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ tff3_inj, data=death_df)
int<-ggsurvplot(km_trt_fit, conf.int = T)
surv_pvalue(km_trt_fit)
survdiff(Surv(days_to_living, sm_death12m) ~ tff3_inj, data=death_df)
ggsave("pd_death_tff3_inj.png", int$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#uric_acid
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ hyper_uri, data=death_df)
uri<-ggsurvplot(km_trt_fit, conf.int = T)
surv_pvalue(km_trt_fit)
survdiff(Surv(days_to_living, sm_death12m) ~ hyper_uri, data=death_df)
ggsave("pd_death_hyper_uri.png", uri$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#lacidosis
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ lacidosis, data=death_df)
lacid<-ggsurvplot(km_trt_fit, conf.int = T)
surv_pvalue(km_trt_fit)
survdiff(Surv(days_to_living, sm_death12m) ~ lacidosis, data=death_df)
ggsave("pd_death_lacidosis.png",lacid$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#haki
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ haki, data=death_df)
aki<-ggsurvplot(km_trt_fit, conf.int = T)
surv_pvalue(km_trt_fit)
survdiff(Surv(days_to_living, sm_death12m) ~ haki, data=death_df)
ggsave("pd_death_haki.png", aki$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#sex
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ sex, data=death_df)
sex<-ggsurvplot(km_trt_fit, conf.int = T)
surv_pvalue(km_trt_fit)
survdiff(Surv(days_to_living, sm_death12m) ~ sex, data=death_df)
ggsave("pd_death_sex.png", sex$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#site
km_trt_fit <- survfit(Surv(days_to_living, sm_death12m) ~ site, data=death_df)
site<-ggsurvplot(km_trt_fit, conf.int = T)
surv_pvalue(km_trt_fit)
survdiff(Surv(days_to_living, sm_death12m) ~ site, data=death_df)
ggsave("pd_death_site.png", site$plot, width = 6, height = 4, dpi = 1200, device = "png") 
#Enterob
plot_data <- asv_sick %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Family") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()
death_df <- plot_data %>% drop_na(sm_death12m) %>%
  mutate(days_to_living = ifelse(is.na(deathdaysto12m), 365, deathdaysto12m),
         doe_date = dmy(doe),
         admission_month = month(doe_date),
         high_enterob = case_when(Enterobacteriaceae >= 7.3 ~ 1,
                                  Enterobacteriaceae < 7.3 ~ 0,
                                  .default = NA)) %>%
  drop_na(high_enterob)
km_enterob<- survfit(Surv(days_to_living, sm_death12m) ~ high_enterob, data=death_df)
g1<-ggsurvplot(km_enterob, conf.int = T, palette = c("darkgrey", "orange"), ylim = c(0.9, 1), xlim = c(0,365))
surv_pvalue(km_enterob)
ggsave("pd_death_Enterob_median.pdf", g1$plot, width = 5, height = 4, dpi = 600, device = "pdf") 

##cox regression
#cox regression with Enterobacteriaceae


plot_data <- asv_sick %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Family") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()

death_df <- plot_data %>% drop_na(sm_death12m) %>%
  mutate(days_to_living = ifelse(is.na(deathdaysto12m), 365, deathdaysto12m),
         doe_date = dmy(doe),
         admission_month = month(doe_date),
         Sex = case_when(sex ==1 ~ "Male",
                         sex ==2 ~ "Female"),
         Site = case_when(site == 1 ~ "Mulago",
                          site == 2 ~ "Jinja"),
         hbsgt_bin = case_when(hbsgt %in% c(1,2) ~ "AA/AS",
                               hbsgt == 3 ~ "SS"),
         age_in_months = round(age * 12),
         sma = case_when(enr_sma == 0 ~ "No",
                         enr_sma == 1 ~ "Yes",
                         .default = NA),
         h_uric = case_when(hyper_uri == 1 ~ "Yes",
                            hyper_uri == 0 ~ "No",
                            .default = NA),
         aki = case_when(haki == 1 ~ "Yes",
                         haki == 0 ~ "No"),
         hbsgt_sma = case_when(hbsgt_bin == "SS" & sma == "Yes" ~ "HbSS_SMA",
                               hbsgt_bin == "AA/AS" & sma == "Yes" ~ "HbAA/AS_SMA",
                               sma == "No" ~ "absent_SMA"),
         high_enterob = case_when(Enterobacteriaceae < 7.3 ~ "< median",
                                  Enterobacteriaceae >= 7.3 ~ ">= median",
                                  .default = NA)) 
#models
entero <- coxph(Surv(days_to_living, sm_death12m) ~ high_enterob, data = death_df)
entero %>% tbl_regression(exponentiate = TRUE) %>%  
  bold_p() %>% # bold p-values under a given threshold (default 0.05)
  bold_labels() %>%
  add_nevent(location="level") %>%
  add_n(location="level") %>%
  as_gt() %>%
  gtsave("high_entero_cox.docx")

entero_hbsgt <- coxph(Surv(days_to_living, sm_death12m) ~ high_enterob+hbsgt_bin, data = death_df)
entero_hbsgt %>% tbl_regression(exponentiate = TRUE) %>%  
  bold_p() %>% # bold p-values under a given threshold (default 0.05)
  bold_labels() %>%
  add_nevent(location="level") %>%
  add_n(location="level") %>%
  as_gt()%>%
  gtsave("high_entero_hbsgt_cox.docx") 

sma_entero <- coxph(Surv(days_to_living, sm_death12m) ~ high_enterob+sma, data = death_df)
sma_entero %>% tbl_regression(exponentiate = TRUE) %>%  
  bold_p() %>% # bold p-values under a given threshold (default 0.05)
  bold_labels() %>%
  add_nevent(location="level") %>%
  add_n(location="level") %>%
  as_gt() %>%
  gtsave("sma_noSS_cox.docx")

uric_entero <- coxph(Surv(days_to_living, sm_death12m) ~ high_enterob+h_uric, data = death_df)
uric_entero %>% tbl_regression(exponentiate = TRUE) %>%  
  bold_p() %>% # bold p-values under a given threshold (default 0.05)
  bold_labels() %>%
  add_nevent(location="level") %>%
  add_n(location="level") %>%
  as_gt()%>%
  gtsave("huric_noSS_cox.docx")

aki_entero <- coxph(Surv(days_to_living, sm_death12m) ~ high_enterob+aki, data = death_df)
aki_entero %>% tbl_regression(exponentiate = TRUE) %>%  
  bold_p() %>% # bold p-values under a given threshold (default 0.05)
  bold_labels() %>%
  add_nevent(location="level") %>%
  add_n(location="level") %>%
  as_gt() %>%
  gtsave("aki_noSS_cox.docx")


#1: testing for separation
#create table to see if any of the levels of categorical variables are fully censored
table(death_df$sm_death12m, death_df$high_enterob)
table(death_df$sm_death12m, death_df$hbsgt_bin)
table(death_df$sm_death12m, death_df$sma)
table(death_df$sm_death12m, death_df$h_uric)
table(death_df$sm_death12m, death_df$aki)
table(death_df$high_enterob, death_df$sma)
table(death_df$high_enterob, death_df$aki)
table(death_df$high_enterob, death_df$h_uric)
table(death_df$high_enterob, death_df$hbsgt_bin)
#2: test for colinearity 
#use the variance inflation factor to diagnose collinearity
#too large magnitude is not defined but values > 2.5 are concerning and value > 5 or 10 indicate collinearity issues
#cox regression does not have VIF built in, so convert to linear regression and then test.
pretend.lm <- lm(sm_death12m ~high_enterob+sma+hbsgt_bin,
                 data = death_df)
car::vif(pretend.lm)
#fine
pretend.lm <- lm(sm_death12m ~ high_enterob+h_uric+hbsgt_bin,
                 data = death_df)
car::vif(pretend.lm)
#fine
pretend.lm <- lm(sm_death12m ~ high_enterob+aki+hbsgt_bin,
                 data = death_df)
car::vif(pretend.lm)
#fine


#3:test proportional hazard ratio
#cox regression assumes the hazard ratio of a predictor is constant over time-- hazard itself may change over time, but the ratio of hazards comparing any two individuals should be constant
#cox.zph() works by comparing PH model for each predictor to a model where the predictor's regression coeff vary smoothly with time
entero.zph<-cox.zph(entero)
print(entero.zph)
par(mfrow=c(1,2))
plot(entero.zph)
#globally fine 

entero_hbsgt.zph<-cox.zph(entero_hbsgt)
print(entero_hbsgt.zph)
par(mfrow=c(1,2))
plot(entero_hbsgt.zph)
#globally fine 

sma_entero_zph <- cox.zph(sma_entero)
print(sma_entero_zph)
par(mfrow=c(2,2))
plot(sma_entero_zph)
#globally fine

uric_entero_zph <- cox.zph(uric_entero)
print(uric_entero_zph)
par(mfrow=c(2,2))
plot(uric_entero_zph)
#globally fine

aki_entero_zph <- cox.zph(aki_entero)
print(aki_entero_zph)
par(mfrow=c(2,2))
plot(aki_entero_zph)
#globally fine

#check residulas for linearity
#display more dianostics with ggcoxdiagnostics 
ggcoxdiagnostics(entero, type = "dfbeta")
ggcoxdiagnostics(entero, type = "schoenfeld")
ggcoxdiagnostics(entero, type = "martingale")

#validate findings
library(rms)
dd <- datadist(death_df[,c("high_enterob")])
options(datadist='dd')
time.inc <- 365
S <- rms::Surv(death_df$days_to_living, death_df$sm_death12m)
f <- cph(S ~ high_enterob + hbsgt_bin, x=TRUE, y=TRUE, surv=TRUE, time.inc = time.inc,
         data = death_df)
summary(f)
#calibrate
set.seed(1234)
cal <- calibrate(f, cmethod = "hare", method="boot", m=100, u=time.inc, B=300)
par(mfrow=c(1,1))
plot <- plot(cal, ylim=c(.92, 1), subtitles=FALSE)
calkm <- calibrate(f, u=time.inc, m=360,  cmethod='KM', B=300)
plot <- plot(calkm, add=TRUE)  


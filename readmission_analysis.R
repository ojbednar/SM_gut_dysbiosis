library(tidyverse)
library(ggthemes)
library(ggpubr)
library(gt)
library(gtsummary)

new <- read.csv("full_study.csv")
new <- new %>% mutate(sick = case_when(group != 6 ~ "SM",
                                       group == 6 ~ "CC")) 
inpat <- new %>%
  drop_na(inpatmal12m) %>%
  count(sick, inpatmal12m) %>%       
  group_by(sick) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatmal12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(sick = factor(sick, levels = c("SM", "CC"))) %>%  # Reorder factor levels
  ggplot(aes(sick, pct, fill = sick)) +
  geom_bar(stat = "identity") +
  ylab("% of children with inpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("red", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) + 
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatmal12m_sick_new.pdf", inpat, width = 3, height = 6, dpi = 600, device = "pdf")


inpat_sep <- new %>%
  drop_na(inpatsepsis12m) %>%
  count(sick, inpatsepsis12m) %>%       
  group_by(sick) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatsepsis12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(sick = factor(sick, levels = c("SM", "CC"))) %>%  # Reorder factor levels
  ggplot(aes(sick, pct, fill = sick)) +
  geom_bar(stat = "identity") +
  ylab("% of children with inpatient sepsis over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("red", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatsepsis12m_sick_new.pdf", inpat_sep, width = 3, height = 6, dpi = 600, device = "pdf")

outpat <- new %>%
  drop_na(outpatmal12m) %>%
  count(sick, outpatmal12m) %>%       
  group_by(sick) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(outpatmal12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(sick = factor(sick, levels = c("SM", "CC"))) %>%  # Reorder factor levels
  ggplot(aes(sick, pct, fill = sick)) +
  geom_bar(stat = "identity") +
  ylab("% of children with outpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("red", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("outpatmal12m_sick_new.pdf", outpat, width = 3, height = 6, dpi = 1200, device = "pdf")

inpat_any <- new %>%
  drop_na(inpatany12m) %>%
  count(sick, inpatany12m) %>%       
  group_by(sick) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatany12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(sick = factor(sick, levels = c("SM", "CC"))) %>%  # Reorder factor levels
  ggplot(aes(sick, pct, fill = sick)) +
  geom_bar(stat = "identity") +
  ylab("% of children with any admission over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("red", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatanyl12m_sick_new.pdf", inpat_any,width = 3, height = 6, dpi = 600, device = "pdf")


inpat <- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatmal12m) %>%
  count(Pf, inpatmal12m) %>%       
  group_by(Pf) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatmal12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(Pf = factor(Pf, levels = c("PfPos", "PfNeg"))) %>%  # Reorder factor levels
  ggplot(aes(Pf, pct, fill = Pf)) +
  geom_bar(stat = "identity") +
  ylab("% of children with inpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("lightskyblue2", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatmal12m_cc_new.pdf", inpat, width = 3, height = 6, dpi = 600, device = "pdf")


outpat <- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(outpatmal12m) %>%
  count(Pf, outpatmal12m) %>%       
  group_by(Pf) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(outpatmal12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(Pf = factor(Pf, levels = c("PfPos", "PfNeg"))) %>%  # Reorder factor levels
  ggplot(aes(Pf, pct, fill = Pf)) +
  geom_bar(stat = "identity") +
  ylab("% of children with outpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("lightskyblue2", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("outpatmal12m_cc_new.pdf", outpat, width = 3, height = 6, dpi = 600, device = "pdf")


inpat_sep <- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatsepsis12m) %>%
  count(Pf, inpatsepsis12m) %>%       
  group_by(Pf) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatsepsis12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(Pf = factor(Pf, levels = c("PfPos", "PfNeg"))) %>%  # Reorder factor levels
  ggplot(aes(Pf, pct, fill = Pf)) +
  geom_bar(stat = "identity") +
  ylab("% of children with sepsis over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 2.5),
            size = 7,
            fontface = "bold", 
            color = "black") +
  scale_fill_manual(values = c("lightskyblue2", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatsepsis12m_cc_new.pdf", inpat_sep, width = 3, height = 6, dpi = 600, device = "pdf")

inpat_any <- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatany12m) %>%
  count(Pf, inpatany12m) %>%       
  group_by(Pf) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatany12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(Pf = factor(Pf, levels = c("PfPos", "PfNeg"))) %>%  # Reorder factor levels
  ggplot(aes(Pf, pct, fill = Pf)) +
  geom_bar(stat = "identity") +
  ylab("% of children with any admission over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("lightskyblue2", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatanyl12m_cc_new.pdf", inpat_any, width = 3, height = 6, dpi = 600, device = "pdf")




#for adjusted odds ratios
# function to calculate OR from 2x2 table
or_fun <- function(data, variable, by, ...) {
  table(data[[by]], data[[variable]]) %>%
    fisher.test() %>%
    broom::tidy()
}

#test
or_fun(new, "sick", "inpatsepsis12m")

##for fischer's exact test, two tailed- SM vs CC
summary_sick <- new %>%
  select(sick, outpatmal12m,inpatmal12m, inpatsepsis12m, inpatany12m) %>%
  tbl_summary(by = sick, missing = "no", statistic = all_categorical() ~ "{n} / {N} ({p}%)") %>%
  add_difference(test = everything() ~ or_fun, estimate_fun = everything() ~ style_ratio)%>%
  modify_header(estimate ~ "**Odds Ratio**")
summary_sick %>% add_q(method = "holm") %>% 
  modify_header(q.value ~ "**Adjusted p-value**") %>%
  as_gt() %>%
  gt::gtsave(filename = "OR_cc_v_sm_hospitalization.docx") 

##for fischer's exact test, two tailed- Pfpos vs pfneg
cc <- new %>% filter(group == 6) %>% mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                                                           pcr_fal0 == 0 ~ "PfNeg"))
summary_cc <- cc %>%
  select(Pf, outpatmal12m,inpatmal12m, inpatsepsis12m, inpatany12m) %>%
  tbl_summary(by = Pf, missing = "no", statistic = all_categorical() ~ "{n} / {N} ({p}%)") %>%
  add_difference(test = everything() ~ or_fun, estimate_fun = everything() ~ style_ratio) %>%
  modify_header(estimate ~ "**Odds Ratio**")
summary_cc %>% add_q(method = "holm") %>% 
  modify_header(q.value ~ "**Adjusted p-value**") %>%
  as_gt() %>%
  gt::gtsave(filename = "OR_pfpos_v_pfneg_hospitalization.docx") 
  



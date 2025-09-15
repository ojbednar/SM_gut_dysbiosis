library(microViz)
library(phyloseq)
library(ggthemes)
#alpha diversity 
asv <- readRDS("ndi_asv.rds")
asv_0 <- asv %>% ps_filter(month == 0)
asv_0 <- asv_0 %>% ps_mutate(sick = case_when(group %in% c(1,2,3,4,5) ~ "SM",
                                              group == 6 ~ "CC"),
                             doe_date = dmy(doe),
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
                                                        .default = 0))

#asv_0
alpha_stat <-  asv_0 %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_export <- alpha_stat %>% select(sick, shannon_Species, observed_Species)
write_csv(alpha_export, "asv_0_alpha.csv")
#pfposvpfneg
alpha_stat_cc <-  asv_cc %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_export_cc <- alpha_stat_cc %>% select(Group, shannon_Species, observed_Species)
write_csv(alpha_export_cc, "asv_0_cc_alpha.csv")

#16S SM v CC
beta_1 <- asv_0 %>% 
  ps_mutate(group_order = case_when(group %in% c(1,2,3,4,5) ~ "A",
                                    group == 6 ~ "B")) %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "sick", alpha = 0.6, size = 2) +
  theme_classic(base_size = 15, base_line_size = 0.5) +
  coord_fixed(0.5) + scale_color_calc()

asv_0 %>% tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("sick", "Site", "age_in_months", "prior_abx_oral", "stool_day"), n_perms = 999, seed = 123) %>%
  perm_get()
ggsave("ccvsm_beta_aitchison_asv_0_sm_v_cc.png", beta_1, width = 5, height = 4, dpi = 1200, device = "png")


#16S PfPos vs PfNeg
beta_cc <- asv_0 %>% ps_filter(sick == FALSE) %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Group", alpha = 0.85, size = 3) +
  theme_classic(12) +
  coord_fixed(0.7) + scale_color_brewer(palette = "Paired")
asv_0 %>% ps_filter(sick == FALSE) %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("Group", "site", "sex", "age", "stool_day"), n_perms = 999, seed = 123) %>%
  perm_get()


##alpha
#asv_0
alpha_stat <-  asv_0 %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_export <- alpha_stat %>% select(sick, shannon_Species, observed_Species)
write_csv(alpha_export, "asv_0_alpha.csv")
#pfposvpfneg
alpha_stat_cc <-  asv_cc %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_export_cc <- alpha_stat_cc %>% select(Group, shannon_Species, observed_Species)
write_csv(alpha_export_cc, "asv_0_cc_alpha.csv")


#for whole genome sequencing
library(phyloseq)
library(microViz)
library(corncob)
library(ggthemes)

krak_0 <- readRDS("krak_0.RDS")
metaphlan_0 <-readRDS("metaphlan_0.RDS")

##beta diversity
beta_2 <- metaphlan_0 %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "sick", alpha = 0.6, size = 2) +
  theme_classic(12) +
  coord_fixed(0.7) + #scale_color_brewer(palette = "Set1") +
  labs(title = "Metaphlan") + scale_color_calc()
metaphlan_0 %>% tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("sick", "site", "age", "stool_day"), n_perms = 999, seed = 123) %>%
  perm_get()
ggsave("ccvsm_beta_aitchison_metaphlan_sick.png", beta_2, width = 5, height = 7, dpi = 1200, device = "png")
beta_3 <- krak_0 %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "sick", alpha = 0.6, size = 2) +
  theme_classic(base_size = 15, base_line_size = 0.5) +
  coord_fixed(0.7) + scale_color_calc()+labs(title = "Kraken")
krak_0 %>% tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("sick", "site", "age", "stool_day"), n_perms = 999, seed = 123) %>%
  perm_get()
ggsave("ccvsm_beta_aitchison_kraken_sick.png", beta_3, width = 5, height = 5, dpi = 1200, device = "png")


##alpha diversity table
alpha_stat_krak <-  krak_0 %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_krak_export <- alpha_stat_krak %>% select(sick, shannon_Species, observed_Species) %>% group_by(sick)
write_csv(alpha_krak_export, "kraken_alpha.csv")

alpha_stat_met <-  metaphlan_0 %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_met_export <- alpha_stat_met %>% select(shannon_Species, observed_Species, sick, Group) %>% group_by(sick)
write_csv(alpha_met_export, "met_alpha.csv")

#output tables for species level
met_t <- metaphlan_0 %>% ps_select(c("studyid", "sick"))
plot_data_species_met <- met_t %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Species") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>%
  samdat_tbl()
write_csv(plot_data_species_met, "species_from_met_0.csv")

krak_t <- krak_0 %>% ps_select(c("studyid", "sick"))
plot_data_species_krak <- krak_t %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Species") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% # adds Parabacteroides as sample data!
  samdat_tbl()
write_csv(plot_data_species_krak, "species_from_krak_0.csv")


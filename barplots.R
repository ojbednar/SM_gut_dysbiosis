
##barplots###
#family
asv <- read.csv("ndi_asv.rds")
asv_0 <- asv %>% ps_filter(month == 0)
krak_0 <- readRDS(krak_0, "krak_0.RDS")
metaphlan_0 <-readRDS(metaphlan_0, "metaphlan_0.RDS")

ccvsm_barplot <- asv_0 %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 8, tax_order = c("Enterobacteriaceae","Enterococcaceae",  
                                                                    "Bacteroidaceae", "Tannerellaceae","Lachnospiraceae",
                                                                    "Prevotellaceae","Bifidobacteriaceae","Streptococcaceae")
               ) +
  labs(x= NULL, y = NULL, title = "sick")
ggsave("ccvsm_family_16S.png", ccvsm_barplot, width = 8.3, height = 10, dpi = 1200, device = "png") 

k_barplot <- krak_0 %>%
  ps_mutate(Outcome = cut(outcome, breaks=2, labels=c("survived", "died"))) %>%
  ps_select(Outcome, Group, sick) %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 8, tax_order = c("Enterobacteriaceae","Enterococcaceae",  
                                                                                "Bacteroidaceae", "Tannerellaceae","Lachnospiraceae",
                                                                                "Prevotellaceae","Bifidobacteriaceae","Streptococcaceae")) +
  labs(x= NULL, y = NULL, title= "Kraken")
ggsave("ccvsm_family_kraken.png", k_barplot, width = 8.3, height = 10, dpi = 1200, device = "png")

m_barplot <- metaphlan_0 %>%
  ps_mutate(Outcome = cut(outcome, breaks=2, labels=c("survived", "died"))) %>%
  ps_select(Outcome, Group, sick) %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 8, tax_order = c("Enterobacteriaceae","Enterococcaceae",  
                                                                                "Bacteroidaceae", "Tannerellaceae","Lachnospiraceae",
                                                                                "Prevotellaceae","Bifidobacteriaceae","Streptococcaceae")) +
  labs(x= NULL, y = NULL, title= "Metaphlan sick")
ggsave("ccvsm_family_metaphlan.png", m_barplot, width = 8.3, height = 10, dpi = 1200, device = "png") 

####
#timeline
ccvsm_barplot_12 <- asv_12 %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Genus", bar_width = 0.8, n_taxa = 12, tax_order = c("Escherichia","Enterococcus", "Bacteroides",
                                                                                "Parabacteroides", "Klebsiella", "Shigella",
                                                                                "Prevotella", "Faecalibacterium", "Bifidobacterium", 
                                                                                "Blautia", "Ruminococcus", "Roseburia")) +
  labs(x= NULL, y = NULL, title = "Month 12")
ccvsm_barplot <- asv_0 %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Genus", bar_width = 0.8, n_taxa = 12, tax_order = c("Escherichia","Enterococcus", "Bacteroides",
                                                                                "Parabacteroides", "Klebsiella", "Shigella", 
                                                                                "Prevotella", "Faecalibacterium", "Bifidobacterium", 
                                                                                "Blautia", "Ruminococcus", "Roseburia")) +
  labs(x= NULL, y = NULL, title = "Enrollement")
ccvsm_barplot_1 <- asv_1 %>% ps_filter(month ==1) %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Genus", bar_width = 0.8, n_taxa = 12, tax_order = c("Escherichia","Enterococcus", "Bacteroides",
                                                                                "Parabacteroides", "Klebsiella", "Shigella",
                                                                                "Prevotella", "Faecalibacterium", "Bifidobacterium", 
                                                                                "Blautia", "Ruminococcus", "Roseburia")) +
  labs(x= NULL, y = NULL, title = "Enrollement")
ggsave("ccvsm_genus_month_1.png", ccvsm_barplot_1, width = 3, height =6 , dpi = 1200, device = "png")

#cef use family
asv_sick <- asv_sick %>% ps_mutate(no_ceftri1 = case_when(ceftri1 == 1 ~ "ceftri1", .default = "no ceftri1"))
ceftri_barplot <- asv_sick %>%
  phyloseq::merge_samples(group = "no_ceftri1") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 8, tax_order = c("Enterobacteriaceae","Enterococcaceae",  
                                                                                "Bacteroidaceae", "Tannerellaceae","Lachnospiraceae",
                                                                                "Prevotellaceae","Bifidobacteriaceae","Streptococcaceae"))
ceftri_barplot
ggsave("ceftri_vs_no_ceftri_SM_16S.png", ceftri_barplot, width = 8.3, height = 10, dpi = 1200, device = "png") 

#death for month1
death_barplot_1 <- asv_1 %>% ps_filter(month ==1) %>%
  phyloseq::merge_samples(group = "sm_death12m") %>%
  comp_barplot(tax_level = "Genus", bar_width = 0.8, n_taxa = 12, tax_order = c("Escherichia","Enterococcus", "Bacteroides",
                                                                                "Parabacteroides", "Klebsiella", "Shigella",
                                                                                "Prevotella", "Faecalibacterium", "Bifidobacterium", 
                                                                                "Blautia", "Ruminococcus", "Roseburia")) +
  labs(x= NULL, y = NULL, title = "Month 1 Mortality")
ggsave("death_genus_month_1.pdf", death_barplot_1, width = 3, height =3.5 , dpi = 600, device = "pdf")

library(gridExtra)
both <- grid.arrange(ccvsm_barplot, ccvsm_barplot_12, ncol = 2)
ggsave("ccvsm_genus_enroll_v_12m.png", both, width = 8, height =6 , dpi = 1200, device = "png")



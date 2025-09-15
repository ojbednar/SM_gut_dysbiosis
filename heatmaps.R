
###HEATMAP####
library(microViz)
library(tidyverse)
#antibio2
asv <- readRDS("ndi_asv.rds")
asv_sick <- asv %>% ps_filter(month == 0,sick == TRUE) 
asv_sick <- asv_sick %>% ps_mutate(prior_abx_oral = case_when(amoxyl == 1 ~ 1,
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
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("prior_abx_oral")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
antibio2<- bb_stats_asv_0_group %>% taxatree_stats_get()

#feed
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("feedhour")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
feed<- bb_stats_asv_0_group %>% taxatree_stats_get()

#neu
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("neuabs")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
neu <- bb_stats_asv_0_group %>% taxatree_stats_get()

#wbc
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("wbc")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
wbc <- bb_stats_asv_0_group %>% taxatree_stats_get()

#ho-1
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("ho1_pl1")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
ho <- bb_stats_asv_0_group %>% taxatree_stats_get()

#free hemoglobin
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("hb_se1")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
hb <- bb_stats_asv_0_group %>% taxatree_stats_get()


#hemin
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("hm_se1")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
hm <- bb_stats_asv_0_group %>% taxatree_stats_get()


#hyper_uri
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("uric_se1")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
uri <- bb_stats_asv_0_group %>% taxatree_stats_get()

#heatmap
h <- rbind(antibio2, feed, neu, wbc, hb, hm, ho, uri) # bind all variables vertically
h <- h %>% filter(rank == "Genus")#filter to only Genus
h$taxon <- sub('^G:', '', h$taxon)#remove prefix
h$taxon <- reorder(h$taxon,h$t.statistic)#order by significance


library(viridis)
library(ggplot2)
p<-ggplot(h, aes(term, taxon, fill= t.statistic)) + 
  # scale_x_discrete(limits=(h$term)[order(h$t.statistic)]) + 
  #scale_y_discrete(limits=(h$taxon)[order(h$t.statistic)]) +
  geom_tile() 
p <- p +scale_fill_distiller(palette = 'RdBu') + geom_text(aes(label=if_else(p.value<0.05, "*", " ")), hjust=0.6, vjust=0.75) +  coord_fixed(ratio=0.2) + 
  scale_x_discrete(expand = c(-1,1)) + theme_grey()
p
ggsave("ordered_factors_contributing_w_antibio2.png", plot = p)

#ceftriaxone heatmap
#inpatient ceftriaxone between SM
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("ceftri1")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
cef<- bb_stats_asv_0_group %>% taxatree_stats_get()

h <- cef # bind all variables vertically
h <- h %>% filter(rank == "Genus")#filter to only Genus
h$taxon <- sub('^G:', '', h$taxon)#remove prefix
h$taxon <- reorder(h$taxon,h$t.statistic)#order by significance


library(viridis)
library(ggplot2)
p<-ggplot(h, aes(term, taxon, fill= t.statistic)) + 
  # scale_x_discrete(limits=(h$term)[order(h$t.statistic)]) + 
  #scale_y_discrete(limits=(h$taxon)[order(h$t.statistic)]) +
  geom_tile() 
p <- p +scale_fill_distiller(palette = 'RdBu', limits=c(-12, 12)) + geom_text(aes(label=if_else(p.value<0.05, "*", " ")), hjust=0.6, vjust=0.75) +  coord_fixed(ratio=0.2) + 
  scale_x_discrete(expand = c(-1,1)) + theme_grey()
p
ggsave("inpatient_ceftriaxone.png", plot = p, width = 5, height = 7)

#####
#enr_sma
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("enr_sma")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
sev_anem<- bb_stats_asv_0_group %>% taxatree_stats_get()

#Coma
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("enr_cm")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
coma<- bb_stats_asv_0_group %>% taxatree_stats_get()

#lactic acidosis
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("lacidosis")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
lacid <- bb_stats_asv_0_group %>% taxatree_stats_get()

#haki
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("haki")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
haki <- bb_stats_asv_0_group %>% taxatree_stats_get()

#acidotic
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("acidotic")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
acidot <- bb_stats_asv_0_group %>% taxatree_stats_get()


#intestinal injury
bb_models_asv_0_group <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("tff3_inj")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
tff3 <- bb_stats_asv_0_group %>% taxatree_stats_get()

#heatmap
h <- rbind(sev_anem, coma, lacid, acidot, haki, tff3) # bind all variables vertically
h <- h %>% filter(rank == "Genus")#filter to only Genus
h$taxon <- sub('^G:', '', h$taxon)#remove prefix
h$taxon <- reorder(h$taxon,h$t.statistic)#order by significance


library(viridis)
library(ggplot2)
p<-ggplot(h, aes(term, taxon, fill= t.statistic)) + 
  # scale_x_discrete(limits=(h$term)[order(h$t.statistic)]) + 
  #scale_y_discrete(limits=(h$taxon)[order(h$t.statistic)]) +
  geom_tile() 
p <- p +scale_fill_distiller(palette = 'RdBu') + geom_text(aes(label=if_else(p.value<0.05, "*", " "))) +  coord_fixed(ratio=0.3) + 
  scale_x_discrete(expand = c(-1,1)) + theme_grey()
p
ggsave("ordered_complications_heatmap.pdf", plot = p)


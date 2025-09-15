library(phyloseq)
library(tidyverse)
library(microViz)
library(corncob)

asv <- read.csv("ndi_asv.rds")
asv_0 <- asv %>% ps_filter(month == 0)
asv_cc <- asv %>% ps_filter(month == 0)%>% 
  ps_filter(sick == FALSE)
asv_pfpos <- asv %>% ps_filter(Group == "PfPos")
asv_sm <- asv %>% ps_filter(sick == TRUE)

#SM vs CC
bb_models_asv_0_group <- asv_0 %>% 
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("sick")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
bb_stats_asv_0_group %>% taxatree_stats_get()
tree_group <- bb_stats_asv_0_group %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
tree_group 

#PfPos vs PfNeg
bb_models_asv_0_cc <- asv_cc %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("Group")
  )
bb_stats_asv_0_cc <- taxatree_models2stats(bb_models_asv_0_cc, param = "mu")
bb_stats_asv_0_cc
bb_stats_asv_0_cc %>% taxatree_stats_get()
tree_cc <- bb_stats_asv_0_cc %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 

tree_cc

key <- taxatree_plotkey(
  data = bb_stats_asv_0_cc,
  taxon_renamer = function(x) stringr::str_remove(x, "[PCFGS]: "),
  # 2 lines of conditions below, for filtering taxa to be labelled
  rank == "Family" | rank == "Genus" |rank == "Species", #& prevalence > 0.2,
  p.value < 0.05, 
  !grepl("Kingdom", taxon)
) +
  # add a bit more space for the longer labels by expanding the x axis
  scale_x_continuous(expand = expansion(mult = 0.2))

#pfpos month 12 vs month 0
bb_models_asv_pfpos <- asv_pfpos %>% 
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("month")
  )
bb_stats_asv_pfpos<- taxatree_models2stats(bb_models_asv_pfpos, param = "mu")
bb_stats_asv_pfpos
bb_stats_asv_pfpos %>% taxatree_stats_get()
tree_group <- bb_models_asv_pfpos %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
tree_group 

#SM month 12 vs month 0
bb_models_asv_sm <- asv_sm %>% 
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("month")
  )
bb_stats_asv_sm <- taxatree_models2stats(bb_models_asv_sm, param = "mu")
bb_stats_asv_sm
bb_stats_asv_sm %>% taxatree_stats_get()
tree_group <- bb_models_asv_sm %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
tree_group 

#CC vs SM-CTX
asv_0 <- asv_0 %>% ps_mutate(no_ceftri1 = case_when(ceftri1 == 1 ~ "ceftri1", .default = "no ceftri1"),
                             sm_no_ceftri = case_when(no_ceftri1 == "no ceftri1" & sick == TRUE ~ 1,
                                                      .default = 0)) 
asv_cc_w_SM_no_ceftri <- asv_0 %>% ps_filter(no_ceftri1 == "no ceftri1")

bb_models_asv_0_group <-  asv_cc_w_SM_no_ceftri %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("sick")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
no_ceft_sm_vs_CC <- bb_stats_asv_0_group %>% taxatree_stats_get()

tree_group <- bb_stats_asv_0_group %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
tree_group 

#all mortality
all_mort <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family"),
    variables = c("death")
  )
all_mort <- taxatree_models2stats(all_mort, param = "mu")
all_mort
all_mort_stats <- all_mort %>% taxatree_stats_get()

#tax tree
tree_group <- all_mort %>%
  taxatree_plots(
    node_size_range = c(1, 4), colour_trans = "identity", palette = "Blue-Red 3"
  ) %>%
  # keep only first 4 plots
  .[1-4] %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) #+ labs(title = "death in children with inpat abx")

tree_group 
ggsave("16S_corncob_SM_mort_circle_tree_family.pdf", tree_group, width = 13, height = 5.5, dpi = 600, device = "pdf")

#anemia
anemia <- asv_sick %>% ps_filter(hgbbase <= 7)
ane_mort <- anemia %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    variables = c("death")
  )
ane_mort <- taxatree_models2stats(ane_mort, param = "mu")
ane_mort
ane_mort_stats <- ane_mort %>% taxatree_stats_get()

#tax tree
tree_group <- ane_mort %>%
  taxatree_plots(
    node_size_range = c(1, 4), colour_trans = "identity", palette = "Blue-Red 3"
  ) %>%
  # keep only first 4 plots
  .[1-4] %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) #+ labs(title = "death in children with inpat abx")

tree_group

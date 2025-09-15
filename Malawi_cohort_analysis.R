library(phyloseq)
library(microViz)
library(corncob)
library(ggthemes)
library(tidyverse)

#read in bracken output
#import file
krak_mz <- read.csv("bracken_mpa_merge_MW.csv")
view(krak_mz)

##clean up merged reports
#remove .breport from name
krak_mz <- krak_mz %>% rename_with(~ str_remove(., "\\.breport$"))
view(krak_mz)

#select only species level analysis to replicate OTU table
krak_mz.species <- krak_mz %>% filter(grepl('s__', Classification))
view(krak_mz.species)

#clean up double kingdom lists
clean_taxonomy <- function(x) {
  kingdoms <- str_extract_all(x, "k__[^|]+")[[1]]  # Extract all k__ terms
  
  if (length(kingdoms) > 1) {  # If there are multiple k__ terms
    merged_kingdom <- paste(unique(kingdoms), collapse = "_")  # Merge unique k__ terms with "_"
    x <- str_replace(x, paste(kingdoms, collapse = "\\|"), merged_kingdom)  # Replace in the original string
  }
  
  return(x)
}
# Apply function to vector
krak_mz.species<- krak_mz.species %>%
  mutate(Classification = sapply(Classification, clean_taxonomy))

#make classification rownames
krak_mz.species <- krak_mz.species %>% mutate(taxa= Classification) %>%
  column_to_rownames(var = "taxa")

#create otu table
otu <- subset(krak_mz.species, select = -Classification) %>% as.matrix()
otu_species = otu_table(otu, taxa_are_rows = TRUE)

#create tax table
#tax table
#has to be complicated as not levels were filled in the bracken file
krak_mz.tax <- krak_mz.species %>% select(Classification)

extract_taxonomy <- function(tax_string, prefix) {
  match <- grep(paste0("^", prefix), tax_string)
  if (length(match) > 0) {
    return(gsub(paste0(prefix, "__"), "", tax_string[match]))
  } else {
    return(NA)
  }
}

split_data <- krak_mz.tax %>%
  # Step 1: Split the classification into individual components based on "|"
  separate('Classification', into = paste0("Taxon", 1:7), sep = "\\|", fill = "right", remove = FALSE) %>%
  
  # Step 2: Extract and assign taxonomic ranks using custom function
  rowwise() %>%
  mutate(
    Domain = extract_taxonomy(c_across(starts_with("Taxon")), "k"),
    Phylum = extract_taxonomy(c_across(starts_with("Taxon")), "p"),
    Class = extract_taxonomy(c_across(starts_with("Taxon")), "c"),
    Order = extract_taxonomy(c_across(starts_with("Taxon")), "o"),
    Family = extract_taxonomy(c_across(starts_with("Taxon")), "f"),
    Genus = extract_taxonomy(c_across(starts_with("Taxon")), "g"),
    Species = extract_taxonomy(c_across(starts_with("Taxon")), "s")
  ) %>%
  ungroup() %>%
  select(-starts_with("Taxon")) %>%
  column_to_rownames(var="Classification")%>%
  mutate(across(where(is.character), ~ str_replace_all(.x, "_", " "))) %>% #remove underscores
  as.matrix()
#view(split_data)
tax_krak_mz <- tax_table(split_data)

#now for metadata

mal_meta <- read.csv("malawi_wgs_meta_git.csv")
mal_meta 
rownames(mal_meta) <- mal_meta$seq_ID
mal_meta <- mal_meta %>% select(-seq_ID)
mal_meta <- sample_data(mal_meta)

#now make phyloseq objects
#merge
k_mz_phylo <- merge_phyloseq(otu_species, mal_meta, tax_krak_mz)
k_mz_phylo
#remove any human data that sneaked through previous filtering
k_mz_phylo <-k_mz_phylo %>% tax_select(tax_list = "Homo_sapiens", strict_matches = FALSE, deselect = TRUE)
#now filter to a reasonable amount of reads, fix OTU, and rename rows
k_mz_phylo <- k_mz_phylo %>% tax_filter(min_prevalence = 0.01, min_total_abundance = 1000) %>%
  tax_fix() %>%
  tax_rename(rank = "Species")
k_mz_phylo

#remove other species now if warranted
k_mz_phylo <-k_mz_phylo %>% tax_select(tax_list = "Virus", strict_matches = FALSE, deselect = TRUE)
k_mz_phylo <-k_mz_phylo %>% tax_select(tax_list = "Eukaryota", strict_matches = FALSE, deselect = TRUE)
k_mz_phylo <-k_mz_phylo %>% tax_select(tax_list = "Archae", strict_matches = FALSE, deselect = TRUE)
k_mz_phylo

#diversity
#betadiv
beta_1 <- k_mz_phylo %>% 
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "SM", alpha = 0.6, size = 2) +
  theme_classic(12) +
  coord_fixed(0.5) + scale_color_calc()
beta_1
ggsave("malawi_beta_bacteria.png", beta_1, width = 5, height = 4, dpi = 1200, device = "png")
k_mz_phylo %>% tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("SM", "age_yr_exact", "sex", "prior_antib"), n_perms = 999, seed = 123) %>%
  perm_get()

#alpha div
alpha_stat <-  k_mz_phylo %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_export <- alpha_stat %>% select(SM, shannon_Species, observed_Species)
write_csv(alpha_export, "malawi_kraken_alpha.csv")
#remaining analysis in prism

#barplot
#with specific families
ccvsm_barplot <- k_mz_phylo %>%
  phyloseq::merge_samples(group = "SM") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 8, tax_order = c("Enterobacteriaceae","Enterococcaceae",  
                                                                                "Bacteroidaceae", "Tannerellaceae","Lachnospiraceae",
                                                                                "Prevotellaceae","Bifidobacteriaceae","Streptococcaceae")
  ) +
  labs(x= NULL, y = NULL, title = "sick")
ggsave("ccvsm_family_MW_top.png", ccvsm_barplot, width = 8.3, height = 10, dpi = 1200, device = "png") 

#prep differential abundance analysis
#filter to more abundant bacteria 
k_mz_phylo_filtered <- k_mz_phylo %>% tax_filter(min_prevalence = 0.1, min_total_abundance = 10000)
k_mz_phylo_filtered

ccvsm <- k_mz_phylo_filtered %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    variables = c("SM")
  )
ccvsm <- taxatree_models2stats(ccvsm, param = "mu")
ccvsm
ccvsm_stats <- ccvsm %>% taxatree_stats_get()

#basic table starter
library(gt)
library(gtsummary)
library(gtable)

#family level table
ccvsm_stats_out_family <- ccvsm_stats %>%
  filter(rank == "Family") %>%
  mutate(taxon = str_replace(taxon, "^F:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS Differential Abundance Analysis (Malawi cohort): Family Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

ccvsm_stats_out_family

gtsave(ccvsm_stats_out_family, "DA_malawi_ASvSM_family.docx")
#genus level DA tree plot and table

#SM vs CC tree plot
bb_models_k_mz_phylo_group <- k_mz_phylo_filtered %>% 
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("SM")
  )
bb_stats_k_mz_phylo_group <- taxatree_models2stats(bb_models_k_mz_phylo_group, param = "mu")
bb_stats_k_mz_phylo_group
bb_stats_k_mz_phylo_group %>% taxatree_stats_get()
tree_group <- bb_stats_k_mz_phylo_group %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
tree_group 

key <- taxatree_plotkey(
  data = bb_stats_k_mz_phylo_group,
  taxon_renamer = function(x) stringr::str_remove(x, "[PCFGS]: "),
  # 2 lines of conditions below, for filtering taxa to be labelled
  rank == "Genus", #& prevalence > 0.2,
  p.value < 0.05, 
  !grepl("Kingdom", taxon)
) +
  # add a bit more space for the longer labels by expanding the x axis
  scale_x_continuous(expand = expansion(mult = 0.2))

patchwork::wrap_plots(key, tree_group, nrow = 1, heights = 8:7)
ggsave("malawi_filtered_bacteria_circle_tree_genus.png", tree_group, width = 13, height = 5.5, dpi = 1200, device = "png")

#genus table
ccvsm_stats_out_genus <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS Differential Abundance Analysis (Malawi cohort): Genus Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

ccvsm_stats_out_genus

gtsave(ccvsm_stats_out_genus, "DA_malawi_ASvSM_genus.docx")

#species table
ccvsm_stats_out_species <- ccvsm_stats %>%
  filter(rank == "Species") %>%
  mutate(taxon = str_replace(taxon, "^S:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS Differential Abundance Analysis (Malawi cohort): Species Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

ccvsm_stats_out_species
gtsave(ccvsm_stats_out_species, "DA_malawi_ASvSM_species.docx")

#for species prism plot
malawi_dat <- k_mz_phylo_filtered %>%
  ps_select(SM) %>%
  tax_transform("compositional", rank = "Species") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()
write.csv(malawi_dat, "malawi_species.csv")

#demographics table
library(gt)
library(gtsummary)
library(gtable)

mal_m <- read.csv("malawi_wgs_meta_git.csv")
mal <- mal_m %>% select(SM, ageyrs, Sex, admission_month, hfa_z, wfa_z, prior_antib) %>%
  tbl_summary(by="SM",
              missing_text = "(Missing)",
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})")) %>% add_p()
#statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
#add_p() 
mal %>% as_gt() %>% gtsave("malawi_demographics.docx")

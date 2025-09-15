#Uganda Gene families Maaslin2 analysis

# Maaslin2 with HUMAnN3 uniref90 output - Uganda (NDI) Cohort
# 2025-02-21
# sblehman
###
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install(c("Wrench", "phyloseq", "Maaslin2"))

library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(Wrench)

uniref90_file <- "genefamilies_merge_unstratified_UG.tsv"

abundances <- read.csv(uniref90_file, sep="\t")

abundance <- abundances %>%
  dplyr::rename("Feature" = colnames(.)[1]) %>%
  rename_with(~ str_remove(., "_merge_Abundance.RPKs")) %>%
  filter(Feature != "UNMAPPED") %>%
  filter(Feature != "UNGROUPED")

metadata <- read.csv("./ob_sample_names.csv") %>%
  tibble::column_to_rownames("id")

otu_tab <- abundance %>%
  tibble::column_to_rownames("Feature")

tax_tab <- abundance %>%
  dplyr::select(Feature) %>%
  dplyr::mutate(spec_row = Feature) %>%
  tibble::column_to_rownames(var = "spec_row")

ps <- phyloseq(sample_data(metadata),
               otu_table(otu_tab, taxa_are_rows = TRUE),
               tax_table(as.matrix(tax_tab)))

counts <- as.matrix(data.frame(otu_table(ps))) # table of raw counts 
condition_vec <- sample_data(ps)$sick # sample to metadata (using AS vs SM)

W <- wrench(counts, condition = condition_vec) # fitting Wrench model on "sick"

normalizationFactors <- W$nf # retrieving normalization factors for each sample

# applying normalization factors on data:
renorm <- sweep(counts, 2, normalizationFactors, FUN = '/')
ps_norm <- ps
otu_table(ps_norm) <- otu_table(renorm, taxa_are_rows = TRUE)

# running Maaslin2 LM
fit_data_lm <- Maaslin2(input_data = data.frame(otu_table(ps_norm)),
                        input_metadata = data.frame(sample_data(ps_norm)),
                        normalization = "NONE",
                        output = "./maaslin_output/UG_uniref90_wrench_lm",
                        fixed_effects = 'sick',
                        analysis_method = "LM",
                        cores = 8,
                        plot_scatter = FALSE
)

#######
#Malawi Gene families Maaslin2 analysis
###
# Maaslin2 with HUMAnN3 uniref90 output - Malawi Cohort
# 2025-02-21
# sblehman
###


#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install(c("Wrench", "phyloseq", "Maaslin2"))

library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(Wrench)

uniref90_file <- "genefamilies_merge_unstratified_MW.tsv"

abundances <- read.csv(uniref90_file, sep="\t")

abundance <- abundances %>%
  dplyr::rename("Feature" = colnames(.)[1]) %>%
  rename_with(~ str_remove(., "_merged_Abundance.RPKs")) %>%
  filter(Feature != "UNMAPPED") %>%
  filter(Feature != "UNGROUPED")

metadata <- read.csv("./mz_sample_names.csv") %>%
  tibble::column_to_rownames("id")

otu_tab <- abundance %>%
  tibble::column_to_rownames("Feature")

tax_tab <- abundance %>%
  dplyr::select(Feature) %>%
  dplyr::mutate(spec_row = Feature) %>%
  tibble::column_to_rownames(var = "spec_row")

ps <- phyloseq(sample_data(metadata),
               otu_table(otu_tab, taxa_are_rows = TRUE),
               tax_table(as.matrix(tax_tab)))

counts <- as.matrix(data.frame(otu_table(ps))) # table of raw counts 
condition_vec <- sample_data(ps)$sick # sample to metadata (using AS vs SM)

W <- wrench(counts, condition = condition_vec) # fitting Wrench model on "sick"

normalizationFactors <- W$nf # retrieving normalization factors for each sample

# applying normalization factors on data:
renorm <- sweep(counts, 2, normalizationFactors, FUN = '/')
ps_norm <- ps
otu_table(ps_norm) <- otu_table(renorm, taxa_are_rows = TRUE)

# running Maaslin2 LM
fit_data_lm <- Maaslin2(input_data = data.frame(otu_table(ps_norm)),
                        input_metadata = data.frame(sample_data(ps_norm)),
                        normalization = "NONE",
                        output = "./maaslin_output/MW_uniref90_wrench_lm",
                        fixed_effects = 'sick',
                        analysis_method = "LM",
                        cores = 8,
                        plot_scatter = FALSE
)

###Uganda Pathway Analysis
# Maaslin2 with Wrench normalization and LM on all Uganda NDI HUMAnN output
# pathabundance -> metacyc pathway
# sblehman
# 2025-01-30
###

#if(!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("Maaslin2", "phyloseq", "Wrench"))

library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(Wrench)

pathabundance_file <- "pathabundance_merge_unstratified_UG.tsv"

samples_df <- read.csv("ob_sample_names.csv") %>%
  tibble::column_to_rownames("seq_ID")

metadata <- sample_data(samples_df)

run_Wrench_Maaslin_LM <- function(file, metadata, output_folder_name, mpa_flag=FALSE){
  start <- Sys.time()
  print(paste(start, " Starting ", output_folder_name, "..."))
  if (mpa_flag){
    abundance <- file %>%
      dplyr::rename("Feature" = colnames(.)[1]) %>%
      mutate(across(where(is.numeric), function(x) as.integer(x)))
  } else {
    abundance <- read.csv(file,sep="\t") %>%
      dplyr::rename("Feature" = colnames(.)[1]) %>% # renaming pathway column
      rename_with(~ str_remove(.,".RPKs")) %>%
      rename_with(~ str_remove(.,"_merge_Abundance")) %>% # renaming columns
      mutate(across(where(is.numeric), function(x) as.integer(x))) # makes values simpler to view
  }
  otu_tab <- abundance %>%
    tibble::column_to_rownames("Feature")
  
  tax_tab <- abundance %>%
    dplyr::select(Feature) %>%
    dplyr::mutate(spec_row = Feature) %>%
    tibble::column_to_rownames(var = "spec_row")
  
  ps <- phyloseq(sample_data(metadata),
                 otu_table(otu_tab, taxa_are_rows = TRUE),
                 tax_table(as.matrix(tax_tab))
  )
  
  counts <- as.matrix(data.frame(otu_table(ps))) # table of raw counts 
  condition_vec <- sample_data(ps)$sick # sample to metadata (using CC vs SM)
  
  W <- wrench(counts, condition = condition_vec) # fitting Wrench model on "sick"
  
  normalizationFactors <- W$nf # retrieving normalization factors for each sample
  
  # applying normalization factors on data:
  renorm <- sweep(counts, 2, normalizationFactors, FUN = '/')
  #return(renorm)
  mode(renorm) <- "integer" # makes values simpler to view
  ps_norm <- ps
  otu_table(ps_norm) <- otu_table(renorm, taxa_are_rows = TRUE)
  
  fit_data_lm <- Maaslin2(input_data = data.frame(otu_table(ps_norm)),
                          input_metadata = data.frame(sample_data(ps_norm)),
                          normalization = "NONE",
                          output = paste('./maaslin_output/', output_folder_name,sep=""),
                          fixed_effects = 'sick',
                          analysis_method = "LM"
  )
  stop <- Sys.time()
  print(paste(stop, ": ", output_folder_name, " complete"))
  print(paste("Total runtime: ", stop-start))
}

run_Wrench_Maaslin_LM(pathabundance_file, metadata, "metacyc_pwy_wrench_lm")

###Malawi Pathway Analysis
# Maaslin2 with Wrench normalization and LM on all HUMAnN output
# pathabundance -> metacyc pathway
# 2025-01-30
###

#if(!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("Maaslin2", "phyloseq", "Wrench"))

library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(Wrench)

pathabundance_file <- "pathabundance_merge_unstratified_MW.tsv"

malawi_wgs_meta <- read_csv("malawi_wgs_meta_121724(CM_AS_meta).csv")
samples_df <- read.csv("mz_sample_names.csv")
samples_df <- malawi_wgs_meta[c("seq_ID","age")] %>%
  dplyr::rename("id" = "seq_ID") %>%
  dplyr::inner_join(samples_df, by="id")

metadata <- samples_df %>%
  tibble::column_to_rownames("id") %>%
  sample_data(.)

run_Wrench_Maaslin_LM <- function(file, metadata, output_folder_name, fe="sick", ref=NULL, mpa_flag=FALSE){
  start <- Sys.time()
  print(paste(start, " Starting ", output_folder_name, "..."))
  if (mpa_flag){
    abundance <- file %>%
      dplyr::rename("Feature" = colnames(.)[1]) %>%
      mutate(across(where(is.numeric), function(x) as.integer(x)))
  } else {
    abundance <- read.csv(file,sep="\t") %>%
      dplyr::rename("Feature" = colnames(.)[1]) %>% # renaming pathway column
      rename_with(~ str_remove(.,".RPKs")) %>%
      rename_with(~ str_remove(.,"_merged_Abundance")) %>% # renaming columns
      mutate(across(where(is.numeric), function(x) as.integer(x))) # makes values simpler to view
  }
  otu_tab <- abundance %>%
    tibble::column_to_rownames("Feature")
  
  tax_tab <- abundance %>%
    dplyr::select(Feature) %>%
    dplyr::mutate(spec_row = Feature) %>%
    tibble::column_to_rownames(var = "spec_row")
  
  ps <- phyloseq(sample_data(metadata),
                 otu_table(otu_tab, taxa_are_rows = TRUE),
                 tax_table(as.matrix(tax_tab))
  )
  
  counts <- as.matrix(data.frame(otu_table(ps))) # table of raw counts 
  condition_vec <- sample_data(ps)$sick # sample to metadata (using CC vs SM)
  
  W <- wrench(counts, condition = condition_vec) # fitting Wrench model on "sick"
  
  normalizationFactors <- W$nf # retrieving normalization factors for each sample
  
  # applying normalization factors on data:
  renorm <- sweep(counts, 2, normalizationFactors, FUN = '/')
  #return(renorm)
  mode(renorm) <- "integer" # makes values simpler to view
  ps_norm <- ps
  otu_table(ps_norm) <- otu_table(renorm, taxa_are_rows = TRUE)
  
  fit_data_lm <- Maaslin2(input_data = data.frame(otu_table(ps_norm)),
                          input_metadata = data.frame(sample_data(ps_norm)),
                          normalization = "NONE",
                          output = paste('./maaslin_output/', output_folder_name,sep=""),
                          fixed_effects = fe,
                          analysis_method = "LM",
                          reference = ref
  )
  stop <- Sys.time()
  print(paste(stop, ": ", output_folder_name, " complete"))
  print(paste("Total runtime: ", stop-start))
}

run_Wrench_Maaslin_LM(pathabundance_file, metadata, "metacyc_pwy_wrench_lm")

####
#import post Maaslin tables
malawi_path <- read.csv("malawi_all_path_sawyer.csv")
uganda_path <- read.csv("uganda_all_path_sawyer.csv")

#overlap of pathways for venn diagram
join_path <- full_join(malawi_path, uganda_path, by = "feature")
#positive SM assoc
pos_overlap <- join_path %>% filter(MW_coef > 0 & UG_coef > 0 & MW_qval < 0.05 & UG_qval < 0.05)
pos_MW <- join_path %>% filter(MW_coef > 0 & MW_qval < 0.05) %>%
  filter(!feature %in% pos_overlap$feature)
pos_UG <- join_path %>% filter(UG_coef > 0 & UG_qval < 0.05) %>%
  filter(!feature %in% pos_overlap$feature)
#negative assoc
neg_overlap <- join_path %>% filter(MW_coef < 0 & UG_coef < 0 & MW_qval < 0.05 & UG_qval < 0.05)
neg_MW <- join_path %>% filter(MW_coef < 0 & MW_qval < 0.05) %>%
  filter(!feature %in% neg_overlap$feature)
neg_UG <- join_path %>% filter(UG_coef < 0 & UG_qval < 0.05) %>%
  filter(!feature %in% neg_overlap$feature)

#set theme
theme_set(theme_classic(base_size = 15) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

##Uganda pathways
#lable significant
uganda_path <- uganda_path %>% mutate(diff_reg = case_when(UG_qval >= 0.05 ~ "Not sig",
                                                           UG_coef < 0 ~ "Negative assoc",
                                                           UG_coef > 0 ~ "Positive assoc"),
                                      neg_log10_qval = -log10(UG_qval))
#identify top significant pathways
top_neg <- uganda_path %>%
  filter(diff_reg == "Negative assoc") %>%
  arrange(UG_qval) %>%
  slice_head(n = 3)
top_pos <- uganda_path %>%
  filter(diff_reg == "Positive assoc") %>%
  arrange(UG_qval) %>%
  slice_head(n = 3)
df_labels <- bind_rows(top_pos, top_neg)

#pull out significant and scientifically interesting labels from table
top_with_interest <- uganda_path %>% 
  filter(!is.na(pathway) & pathway != "")
#volcano plot
path_v_ug <- ggplot(uganda_path, aes(x=UG_coef, y=-log10(UG_qval), col = diff_reg)) + 
  geom_vline(xintercept = 0, col = "gray", linetype = 'solid') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "darkgray","darkred"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Negative assoc", "Not sig","Positive assoc")) +
  geom_label_repel(data = subset(top_with_interest, diff_reg == "Negative assoc"),
                   aes(label = pathway, segment.alpha = 0.5, fontface = "bold"), size = 4, point.padding = 0,
                   point.size = 0.1,
                   min.segment.length = 0,
                   box.padding = 0.25,
                   label.padding = 0.2,
                   force_pull = 0,
                   xlim = c(-7, -1),
                   ylim = c(3.5, 8),
                   nudge_x = -2,
                   # direction = "y",
                   max.overlaps = Inf,
                   seed=42) +
  geom_label_repel(data = subset(top_with_interest, diff_reg == "Positive assoc"), 
                   aes(label = pathway, segment.alpha = 0.5, fontface = "bold"), size = 4.5, point.padding = 0,
                   point.size = 0.1,
                   min.segment.length = 0,
                   box.padding = 0.3,
                   label.padding = 0.2,
                   force_pull = 0,
                   xlim = c(4, 11),
                   ylim = c(0, 8),
                   nudge_x= 1,
                   max.overlaps = Inf,
                   seed=42) +
  labs(title="Uganda Pathways")+
  coord_cartesian(xlim = c(-8, 15), ylim = c(0,8), clip = "off")+
  theme(legend.position = "none")
path_v_ug
ggsave("uganda_volcano_pathway.png", plot = path_v_ug, height = 6, width = 9)

#pie charts for categorization 
uganda_sig_SM <- uganda_path %>% filter(UG_qval < 0.05, UG_coef > 0)
table(uganda_sig_SM$category)
uganda_sig_AS <- uganda_path %>% filter(UG_qval < 0.05, UG_coef < 0)
table(uganda_sig_AS$category)


##Malawi pathways
#label significant pathways
malawi_path <- malawi_path %>% mutate(diff_reg = case_when(MW_qval >= 0.05 ~ "Not sig",
                                                           MW_coef < 0 ~ "Negative assoc",
                                                           MW_coef > 0 ~ "Positive assoc"),
                                      neg_log10_qval = -log10(MW_qval))
#assess most significant pathways
library(ggrepel)
top_neg <- malawi_path %>%
  filter(diff_reg == "Negative assoc") %>%
  arrange(MW_qval) %>%
  slice_head(n = 5)
top_pos <- malawi_path %>%
  filter(diff_reg == "Positive assoc") %>%
  arrange(MW_qval) %>%
  slice_head(n = 5)
df_labels <- bind_rows(top_pos, top_neg)

#pull out pathways labeled in table
#second method to pull out labels
top_with_interest <- malawi_path %>% 
  filter(!is.na(pathway) & pathway != "")

#volano plot Malawi paths
path_v_mal <- ggplot(malawi_path, aes(x=MW_coef, y=-log10(MW_qval), col = diff_reg)) + 
  geom_vline(xintercept = 0, col = "gray", linetype = 'solid') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "darkgray","darkred"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Negative assoc", "Not sig","Positive assoc")) +
  geom_label_repel(data = subset(top_with_interest, diff_reg == "Negative assoc"),
                   aes(label = pathway, segment.alpha = 0.5, fontface = "bold"), size = 4, point.padding = 0,
                   point.size = 0.1,
                   min.segment.length = 0,
                   box.padding = 0.4,
                   label.padding = 0.2,
                   force_pull = 0,
                   xlim = c(-16,-1),
                   ylim = c(-1, 18),
                   max.overlaps = Inf) +
  geom_label_repel(data = subset(top_with_interest, diff_reg == "Positive assoc"), 
                   aes(label = pathway, segment.alpha = 0.5, fontface="bold"), size = 4.5, point.padding = 0,
                   point.size = 0.1,
                   min.segment.length = 0,
                   box.padding = 0.4,
                   label.padding = 0.2,
                   force_pull = 0,
                   xlim = c(3.5, 15),
                   ylim = c(-1, 18),
                   max.overlaps = Inf) +
  labs(title="Malawi Pathways") +
  coord_cartesian(xlim = c(-15, 15), ylim = c(0,15), clip = "off")+
  theme(legend.position = "none")
path_v_mal
ggsave("malawi_volcano_pathway.png", plot = path_v_mal, height = 6, width = 9)

#pie chart of categorization
malawi_sig_SM <- malawi_path %>% filter(MW_qval < 0.05, MW_coef > 0)
table(malawi_sig_SM$category)
malawi_sig_AS <- malawi_path %>% filter(MW_qval < 0.05, MW_coef < 0)
table(malawi_sig_AS$category)

####Gene families
#Uganda
uniref90_ug <- read.table("uganda_uniref_sawyer.tsv", sep = "\t",  header = TRUE)
#pull out significant
uganda_uniref <- uniref90_ug %>% mutate(diff_reg = case_when(qval >= 0.01 ~ "Not sig",
                                                             coef < -1 ~ "Negative assoc",
                                                             coef > 1 ~ "Positive assoc",
                                                             .default = "Not sig"),
                                        neg_log10_qval = -log10(qval))

library(ggrepel)
# Top 10 pos/neg assoc for labeling
top_neg <- uganda_uniref %>%
  filter(diff_reg == "Negative assoc") %>%
  arrange(qval) %>%
  slice_head(n = 10) %>%
  mutate(prot = case_when(
    feature == "UniRef90_R5GY18" ~ "Uncharacterized firmicutes protein",
    feature == "UniRef90_R6EWV7" ~ "Uncharacterized firmicutes protein",
    feature == "UniRef90_U2JVJ9" ~ "AAA-ATPase-like domain-containing protein",
    feature == "UniRef90_W1WNM2" ~ "DNA-directed RNA polymerase subunit",
    feature == "UniRef90_A0A0V8QF06" ~ "Glutaconyl-CoA decarboxylase subunit beta",
    #feature == "UniRef90_A8SSZ2" ~ "Coprococcus hypothetical protein",
    #feature == "UniRef90_UPI000E4B4DBC" ~ "Clostridiaceae hypothetical protein",
    #feature == "UniRef90_U2KDJ8" ~ "AAA-ATPase-like domain-containing protein",
    #feature == "UniRef90_UPI000E4707E6" ~ "Butyricicoccus hypothetical protein",
    #feature == "UniRef90_A0A3C0IHR9" ~ "Ruminococcus partial cell surface protein"
  ))

top_pos <- uganda_uniref %>%
  filter(diff_reg == "Positive assoc") %>%
  arrange(qval) %>%
  slice_head(n = 10) %>%
  mutate(prot = case_when(
    feature == "UniRef90_P56580" ~ "PTS system glucitol/sorbitol-specific,EIIB",
    feature == "UniRef90_B1X935" ~ "tRNA methylamino-thiouridine biosynth,MnmC", 
    #full name:tRNA 5-methylaminomethyl-2-thiouridine biosynthesis,MnmC
    feature == "UniRef90_E1IZ42" ~ "Escherichia hypothetical protein",
    feature == "UniRef90_E2QQW0" ~ "Racemase",
    feature == "UniRef90_P0AC74" ~ "Evolved beta-galactosidase,subunit beta",
    #feature == "UniRef90_A0A376W0Y4" ~ "Phosphoglycerol transferase I",
    #feature == "UniRef90_Q46948" ~ "Protein/nucleic acid deglycase 3",
    #feature == "UniRef90_A0A377D7P4" ~ "Ornithine decarboxylase",
    #feature == "UniRef90_A0A271R0A6" ~ "Escherichia coli IS6 family transposase",
    feature == "UniRef90_A0A376YCG2" ~ "Formate acetyltransferase 1"
  ))

#pull out interest
interest <- uganda_uniref %>%
  mutate(prot = case_when(
    feature == "UniRef90_A0A395WRJ3" ~ "formate dehydrogenase-N subunit alpha",
    feature == "UniRef90_A0A377K0Z6" ~ "Aminomethyltransferase(GCS protein)",
    feature == "UniRef90_UPI000C19FF2A" ~ "ferric receptor FhuE",
    feature == "UniRef90_P43019" ~ "Superoxide dismutase [Mn]",
    feature == "UniRef90_A0A377W6T3" ~ "CFA/I fimbrial,subunit-C",
    feature == "UniRef90_P45475" ~ "Ubiquinone biosynthesis protein,UbiV",
    feature == "UniRef90_P0AET0" ~ "Acid stress chaperone,HdeA",
    feature == "UniRef90_A8ADV1" ~ "NADH-quinone oxidoreductase,subunit-A",
  )) %>%
  filter(!is.na(prot) & prot != "")

# Combine for labeling
df_labels <- bind_rows(top_pos, top_neg, interest)
#volcano plot
uniref_v_ug<- ggplot(uganda_uniref, aes(x=coef, y=-log10(qval), col = diff_reg)) + 
  geom_vline(xintercept = 0, col = "gray", linetype = 'solid') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1), col = "gray", linetype = 'dashed') +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("#00AFBB", "darkgray","darkred"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Negative assoc", "Not sig","Positive assoc")) +
  labs(title="Uganda Gene Families") +
  geom_label_repel(data = subset(df_labels, diff_reg == "Negative assoc"),
                   aes(label = prot, segment.alpha = 0.5, alpha = 0.001, fontface = "bold"), size = 4, point.padding = 0,
                   point.size = 0.1,
                   min.segment.length = 0,
                   box.padding = 0.2,
                   label.padding = 0.2,
                   force_pull = 0,
                   hjust = 0.5,
                   ylim = c(10, 18),
                   direction = "y",
                   #xlim = c(-Inf, 0),
                   max.overlaps = Inf) +
  geom_label_repel(data = subset(df_labels, diff_reg == "Positive assoc"), 
                   aes(label = prot, segment.alpha = 0.5, alpha = 0.001, fontface = "bold"), size = 4, point.padding = 0,
                   point.size = 0.1,
                   min.segment.length = 0,
                   box.padding = 0.4,
                   label.padding = 0.2,
                   force = 15,
                   direction = "both",
                   hjust = 0.5,
                   xlim = c(5.5, 13),
                   ylim = c(-1, 18),
                   #nudge_x = 2,
                   max.overlaps = Inf,
                   seed = 42) +
  coord_cartesian(xlim = c(-10, 15), ylim=c(0,18), clip = "off")+
  theme(legend.position = "none")
uniref_v_ug
ggsave("uganda_genefamilies_volcano_pathway.png", plot = uniref_v_ug, height = 6, width = 10)

##Malawi Gene Families
uniref90_mal <- read.table("malawi_uniref90.tsv", sep = "\t",  header = TRUE)
#identify significant
malawi_uniref <- uniref90_mal %>% mutate(diff_reg = case_when(qval >= 0.01 ~ "Not sig",
                                                              coef < -1 ~ "Negative assoc",
                                                              coef > 1 ~ "Positive assoc",
                                                              .default = "Not sig"),
                                         neg_log10_qval = -log10(qval))
# Top 10 pos/neg assoc for labeling
top_neg <- malawi_uniref %>%
  filter(diff_reg == "Negative assoc") %>%
  arrange(qval) %>%
  slice_head(n = 10) %>%
  mutate(prot = case_when(
    feature == "UniRef90_R7K1W2" ~ "Clostridium uPF0210 protein",
    feature == "UniRef90_A0A173SMI5" ~ "DNA-directed RNA polymerase subunit P",
    feature == "UniRef90_A0A329TV07" ~ "Faecalibacterium prausnitzii acid-shock protein",
    feature == "UniRef90_A0A3E2W4I6" ~ "Gamma-glutamyltranspeptidase",
    feature == "UniRef90_R6K1C0" ~ "Histidinol dehydrogenase",
    #feature == "UniRef90_A0A174H2X3" ~ "Lachnospira pectinoschiza Na+/alanine symporter",
    #feature == "UniRef90_A0A174WK26" ~ "Prevotella copri Glutamate synthase [NADPH] large chain precursor",
    #feature == "UniRef90_A0A174XU46" ~ "Prevotella copri Ferrous iron transport protein B",
    #feature == "UniRef90_A0A174Y3U7" ~ "Prevotella copri Uncharacterised protein",
    #feature == "UniRef90_A0A176U726" ~ "RNA polymerase sigma factor, sigma-70 family"
  ))

top_pos <- malawi_uniref %>%
  filter(diff_reg == "Positive assoc") %>%
  arrange(qval) %>%
  slice_head(n = 10) %>%
  mutate(prot = case_when(
    feature == "UniRef90_D4HBJ9" ~ "Adenosine deaminase",
    feature == "UniRef90_U7M7P1" ~ "UDP-N-acetylmuramoylalanyl ligase",
    feature == "UniRef90_A0A376VJP0" ~ "KHG/KDPG aldolase",
    feature == "UniRef90_Q6A5S0" ~ "Probable potassium transport system,Kup",
    feature == "UniRef90_Q6A7L9" ~ "Proline--tRNA ligase",
    #feature == "UniRef90_A0A1X3INE4" ~ "Escherichia coli hypothetical protein",
    #feature == "UniRef90_A0A075M9R5" ~ "Escherichia coli uncharacterized protein",
    #feature == "UniRef90_D4HBD3" ~ "Electron transfer flavoprotein, carnitine metabolism (FixB protein)",
    #feature == "UniRef90_Q6A6K3" ~ "Large ribosomal subunit protein uL1",
    #feature == "UniRef90_A0A2H5C945" ~ "Eschecherichia coli uncharacterized protein"
  ))

#pull out interest
interest <- malawi_uniref %>%
  mutate(prot = case_when(
    feature == "UniRef90_A0A395WRJ3" ~ "formate dehydrogenase-N,subunit-alpha",
    feature == "UniRef90_A0A377K0Z6" ~ "Aminomethyltransferase(GCS protein)",
    feature == "UniRef90_UPI000C19FF2A" ~ "ferric receptor FhuE",
    feature == "UniRef90_P43019" ~ "Superoxide dismutase [Mn]",
    feature == "UniRef90_A0A377W6T3" ~ "CFA/I fimbrial,subunit-C",
    feature == "UniRef90_P45475" ~ "Ubiquinone biosynthesis,UbiV",
    feature == "UniRef90_P0AET0" ~ "Acid stress chaperone,HdeA",
    feature == "UniRef90_A8ADV1" ~ "NADH-quinone oxidoreductase,subunit-A",
  )) %>%
  filter(!is.na(prot) & prot != "")

# Combine for labeling
df_labels <- bind_rows(top_pos, top_neg, interest)


uniref_v_mal <- ggplot(malawi_uniref, aes(x=coef, y=-log10(qval), col = diff_reg)) + 
  geom_vline(xintercept = 0, col = "gray", linetype = 'solid') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1), col = "gray", linetype = 'dashed') +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("#00AFBB", "darkgray","darkred"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Negative assoc", "Not sig","Positive assoc")) +
  labs(title="Malawi Gene Families") +
  geom_label_repel(data = subset(df_labels, diff_reg == "Negative assoc"),
                   aes(label = prot, segment.alpha = 0.5, alpha = 0.001, fontface = "bold"), size = 4, point.padding = 0,
                   point.size = 0.1,
                   min.segment.length = 0,
                   box.padding = 0.2,
                   label.padding = 0.2,
                   force_pull = 0,
                   hjust = 0.5,
                   ylim = c(12.5, 18),
                   #xlim = c(-Inf, 0),
                   direction = "y",
                   max.overlaps = Inf) +
  geom_label_repel(data = subset(df_labels, diff_reg == "Positive assoc"), 
                   aes(label = prot, segment.alpha = 0.5, alpha = 0.001, fontface = "bold"), size = 4, point.padding = 0,
                   point.size = 0.1,
                   min.segment.length = 0,
                   box.padding = 0.4,
                   label.padding = 0.2,
                   force_pull = 0,
                   hjust = 0.5,
                   xlim = c(4, 12),
                   ylim = c(-1, 20),
                   direction = "both",
                   max.overlaps = Inf) +
  coord_cartesian(xlim = c(-10, 15), ylim=c(0,18), clip = "off")+
  theme(legend.position = "none")
uniref_v_mal
ggsave("malawi_genefamilies_volcano_pathway.png", plot = uniref_v_mal, height = 6, width = 10)




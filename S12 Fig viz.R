library(here)
library(grid)

#load functions
source(here::here("000_qxone_cluster_polysolver_v3.R"))

#set seed
set.seed(2025)

#Load samples
data_dir <- file.path("/mnt/depts/dept04/compbio/projects/SingleCell/data/ddpcr")
fig4_rep1 <- run_silently(analysis_wrapper(analysis_workbook = file.path(data_dir,"S12 Fig_data.xlsx"),
                                           raw_sheet = "raw",
                                           cluster_sheet = "cluster",
                                           probe_order = "VIC_Cy5_FAM_Cy5.5") )

#Format the data to long for plot
linkage_metadata_df <- primer_order_picker("VIC_Cy5_FAM_Cy5.5")
linked <- linkage_metadata_df %>% filter(experimental_linkage == "linked" | droplet_type == "singlet") %>%
  pull(fragment_combination) %>% gsub("channel","prob",.)

#Retrieve sample description 3 (this is the independent replicate)

meta_data <- read_xlsx(file.path(data_dir,"S12 Fig_data.xlsx"),
                      sheet="raw") %>% 
  janitor::clean_names() %>% 
  select(well, 
         sample = sample_description_1, # note this needs formatting
         sample_dilution = sample_description_2, # note this needs formatting
         sample_description_3) %>% distinct()

#Format for rep1
fig4_rep1_l <- fig4_rep1 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         well_row = gsub("[[:digit:]]","",well),
         well_row_num = match(well_row, LETTERS) ) %>% 
  left_join(meta_data) %>%
  # subset to these columns when building tmp metadata
  # select(starts_with("well_"),
  #        sample_combo = sample,
  #        sample_dilution_description = sample_dilution,
  #        sample_description_3) %>% 
  rename(sample_combo = sample,
         sample_dilution_description = sample_dilution) %>% 
  mutate(rep = case_when(sample_description_3 == "indep prep rep 1" ~ "rep1",
                         sample_description_3 == "indep prep rep2" ~ "rep2"),
         sample_description = case_when(sample_combo == "S2S3" & sample_dilution_description != "100% pCMVEGFP" ~ "pdeltaWPRE_MfeI",
                                        sample_combo == "S1S3" & sample_dilution_description != "100% pCMVEGFP" ~ "pdeltaCMV_MfeI",
                                        sample_combo == "S2S3" & sample_dilution_description == "100% pCMVEGFP" ~ "pCMVEGFP_SphI",
                                        sample_combo == "S1S3" & sample_dilution_description == "100% pCMVEGFP" ~ "pCMVEGFP_SphI",
                            ),
         sample_dilution = case_when( ((well_col %in% c("01", "05")) | 
                                      (well_row == "E" & well_col %in% c("04","08") ) | 
                                      (well_row == "H" & well_col %in% c("04","08") ) ) ~ "dil1", 
                                      
                                      ((well_col %in% c("02","06")) |
                                      (well_row == "D" & well_col %in% c("04","08") ) |
                                      (well_row == "G" & well_col %in% c("04","08") ) ) ~ "dil2",
                                      
                                      ((well_col %in% c("03","07")) |
                                      (well_row == "C" & well_col %in% c("04","08") ) |
                                      (well_row == "F" & well_col %in% c("04","08") ) ) ~ "dil3",
                                      
                                      ((well_row == "A" & well_col  %in% c("04","08") ) |
                                      (well_row == "B" & well_col  %in% c("04","08") ) ) ~ "NTC",
                                      
                                      ),
         
         sample = case_when(well_row == "A" & well_col %in% c("01","02","03","05","06","07") ~ "S8",
                            well_row == "B" & well_col %in% c("01","02","03","05","06","07") ~ "S7",
                            well_row == "C" & well_col %in% c("01","02","03","05","06","07") ~ "S6",
                            well_row == "D" & well_col %in% c("01","02","03","05","06","07") ~ "S5",
                            well_row == "E" & well_col %in% c("01","02","03","05","06","07") ~ "S4",
                            well_row == "F" & well_col %in% c("01","02","03","05","06","07") ~ "S3",
                            well_row == "G" & well_col %in% c("01","02","03","05","06","07") ~ "S2",
                            well_row == "H" & well_col %in% c("01","02","03","05","06","07") ~ "S1",
                            well_row %in% c("C","D","E") & well_col %in% c("04","08") ~ "S10",
                            well_row %in% c("F","G","H") & well_col %in% c("04","08") ~ "S9",
                            TRUE ~ "BLANK"
                            )
         ) %>% 
  #select(starts_with("well_"), sample, sample_description, sample_dilution, sample_dilution_description, rep, everything() )
  #select(-sample_combo, -sample_description_3) %>%
  select(well, sample, sample_description, sample_dilution, sample_dilution_description, rep, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 

#Combine
fig4_comb <- fig4_rep1_l

#Void wells where accepted_droplets is <10k
fig4_comb <- fig4_comb %>% 
  mutate(across(.cols = starts_with("prob_"),
                .fns = ~ if_else(accepted_droplets < 10000,NA,.x) ) )

#Get mean value and nonmissing n
fig4_comb_s <- fig4_comb %>% 
  pivot_longer(cols=starts_with("prob_"), names_to="frag_type", values_to="frag_prob") %>% 
  group_by(sample, sample_dilution, frag_type) %>% 
  summarize(mean_frag_prob = if(all(is.na(frag_prob))) NA_real_ else mean(frag_prob, na.rm=T),
            non_miss_frag_prob = sum(!is.na(frag_prob))
  ) %>% 
  ungroup() 

#Get max FPD and join
fig4_fpd <- fig4_comb %>% 
  select(sample, sample_dilution, max_frag_val) %>% 
  group_by(sample, sample_dilution) %>% 
  summarize(max_frag_val = round(mean(max_frag_val, na.rm=T),0) ) %>% 
  ungroup()

fig4_comb_s <- fig4_comb_s %>% 
  left_join(fig4_fpd, by=c("sample", "sample_dilution")) %>% 
  filter(sample != "Buffer") %>% 
  mutate(sample_dlution = factor(sample_dilution,
                                 levels = paste0("dil",1:3) ),
         #Sample dilutions from original data prep
         sample_dilution_conc = case_when(sample_dilution == "dil1" ~333.33,
                                          sample_dilution == "dil2" ~111.11,
                                          sample_dilution == "dil3" ~37.04,
                                          TRUE ~ NA),
         sample = factor(sample, levels=paste0("S",1:10))
  )

custom_labels1 <- c("S1" = "A. ", 
                    "S2" = "B. ", 
                    "S3" = "C. ", 
                    "S4" = "D. ",
                    "S5" = "E. ", 
                    "S6" = "F. ")

pout1 <- prob_barplot(dfl = fig4_comb_s %>% filter(sample %in% paste0("S",1:5)), 
                      xval = "sample_dilution_conc", 
                      xlabtext = "Plasmid Concentration (copies/ul)", 
                      yval = "mean_frag_prob", 
                      ylabtext = "Probability estimate", 
                      fillval = "frag_type", 
                      filllabtext = "Fragment type", 
                      showlabs=FALSE, #hiding max FPD when combining replicates
                      labval = "max_frag_val",
                      facetval = "sample", 
                      x_breaks = sapply(c(333.33, 111.11, 37.04),round, 0),
                      facet_col = 5
) + 
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0)) + 
  facet_wrap(~sample, labeller = labeller(sample = custom_labels1), nrow = 1)+
  theme(axis.text.x = element_text(size=rel(0.9),
                                   angle=90,
                                   vjust=0.5,
                                   hjust=1) )


custom_labels2 <- c("S6" = "G. ", 
                    "S7" = "H. ", 
                    "S8" = "I. ", 
                    "S9" = "J. ",
                    "S10" = "K. ")

pout2 <- prob_barplot(dfl = fig4_comb_s %>% filter(sample %in% paste0("S",6:10)), 
                      xval = "sample_dilution_conc", 
                      xlabtext = "Plasmid Concentration (copies/ul)", 
                      yval = "mean_frag_prob", 
                      ylabtext = "Probability estimate", 
                      fillval = "frag_type", 
                      filllabtext = "Fragment type", 
                      showlabs=FALSE, #hiding max FPD when combining replicates
                      labval = "max_frag_val",
                      facetval = "sample", 
                      x_breaks = sapply(c(333.33, 111.11, 37.04),round, 0),
                      facet_col = 5
) +  facet_wrap(~sample, labeller = labeller(sample = custom_labels2), nrow = 1)+
  theme(axis.text.x = element_text(size=rel(0.9),
                                   angle=90,
                                   vjust=0.5,
                                   hjust=1) )

#Modify layout
mod <- ggpubr::ggarrange(pout1 + rremove("ylab") + rremove("xlab"), 
                         pout2 + rremove("ylab"), 
                         nrow=2, 
                         common.legend = TRUE, 
                         legend = "right" ) 

mod_annot <- annotate_figure(mod, left=textGrob("Probability estimate", rot=90, vjust=1, gp=gpar(cex=1.1) ) ) 

#Save vector and hires raster
ggsave(plot = mod_annot, filename="fig4_viz_delta_review.svg", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8.5, height=4)

ggsave(plot = mod_annot, filename="fig4_viz_delta_review.png", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8.5, height=4, 
       dpi=600,
       bg="white")

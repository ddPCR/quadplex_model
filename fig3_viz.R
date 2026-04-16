library(here)
library(grid)

#load functions
source(here::here("000_qxone_cluster_polysolver_v3.R"))

#set seed
set.seed(2025)

#Load samples
fig3_rep1 <- run_silently(analysis_wrapper(analysis_workbook = here::here("plasmid_data.xlsx"),
                                           raw_sheet = "fig_2_3_rep1_raw",
                                           cluster_sheet = "fig_2_3_rep1_cluster",
                                           probe_order = "VIC_Cy5_FAM_Cy5.5") )

fig3_rep2 <- run_silently(analysis_wrapper(analysis_workbook = here::here("plasmid_data.xlsx"),
                                           raw_sheet = "fig_2_3_rep2_raw",
                                           cluster_sheet = "fig_2_3_rep2_cluster",
                                           probe_order = "VIC_Cy5_FAM_Cy5.5") )

fig3_rep3 <- run_silently(analysis_wrapper(analysis_workbook = here::here("plasmid_data.xlsx"),
                                           raw_sheet = "fig_2_3_rep3_raw",
                                           cluster_sheet = "fig_2_3_rep3_cluster",
                                           probe_order = "VIC_Cy5_FAM_Cy5.5") )

fig3_rep4 <- run_silently(analysis_wrapper(analysis_workbook = here::here("plasmid_data.xlsx"),
                                           raw_sheet = "fig_2_3_rep4_raw",
                                           cluster_sheet = "fig_2_3_rep4_cluster",
                                           probe_order = "VIC_Cy5_FAM_Cy5.5") )


fwritep(x = fig3_rep1, file = here::here("tmp_files", "fig3_rep1.csv"))
fwritep(x = fig3_rep2, file = here::here("tmp_files", "fig3_rep2.csv"))
fwritep(x = fig3_rep3, file = here::here("tmp_files", "fig3_rep3.csv"))
fwritep(x = fig3_rep4, file = here::here("tmp_files", "fig3_rep4.csv"))

#Format the data to long for plot
linkage_metadata_df <- primer_order_picker("VIC_Cy5_FAM_Cy5.5")
linked <- linkage_metadata_df %>% filter(experimental_linkage == "linked" | droplet_type == "singlet") %>%
  pull(fragment_combination) %>% gsub("channel","prob",.)

#Format for rep1
fig3_rep1_l <- fig3_rep1 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col)) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 

#Format for rep2
fig3_rep2_l <- fig3_rep2 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         sample_dilution = if_else(well_col_num != 12 , paste0("dil",well_col_num), NA) ) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 

#Format for rep3
#Column 12 has samples but dropping due to only 11 dilutions in the original series
fig3_rep3_l <- fig3_rep3 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         sample_dilution = if_else(well_col_num != 12 , paste0("dil",well_col_num), NA),
         sample = if_else(well_col_num == 12, "Buffer", sample)) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num)

#Format for rep4
fig3_rep4_l <- fig3_rep4 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         sample_dilution = if_else(well_col_num != 12 , paste0("dil",well_col_num), NA),
         sample = if_else(well_col_num == 12, "Buffer", sample)) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num)

#Combine
fig3_comb <- rbindlist(list(rep1=fig3_rep1_l,
                            rep2=fig3_rep2_l,
                            rep3=fig3_rep3_l,
                            rep4=fig3_rep4_l),idcol="rep")

#Void wells where accepted_droplets is <10k
fig3_comb <- fig3_comb %>% 
  mutate(across(.cols = starts_with("prob_"),
                .fns = ~ if_else(accepted_droplets < 10000,NA,.x) ) )

fwritep(x = fig3_comb, file = here::here("tmp_files", "fig3_rep_comb.csv"))

#Get mean value and nonmissing n
fig3_comb_s <- fig3_comb %>% 
  pivot_longer(cols=starts_with("prob_"), names_to="frag_type", values_to="frag_prob") %>% 
  group_by(sample, sample_dilution, frag_type) %>% 
  summarize(mean_frag_prob = if(all(is.na(frag_prob))) NA_real_ else mean(frag_prob, na.rm=T),
            non_miss_frag_prob = sum(!is.na(frag_prob))) %>% 
  ungroup() 

#Get max FPD and join
fig3_fpd <- fig3_comb %>% 
  select(sample, sample_dilution, max_frag_val) %>% 
  group_by(sample, sample_dilution) %>% 
  summarize(max_frag_val = round(mean(max_frag_val, na.rm=T),0) ) %>% 
  ungroup()

fig3_comb_s <- fig3_comb_s %>% 
  left_join(fig3_fpd, by=c("sample", "sample_dilution")) %>% 
  filter(sample != "Buffer") %>% 
  mutate(sample_dlution = factor(sample_dilution,
                                 levels = paste0("dil",1:11) ),
         #Sample dilutions from original data prep
         sample_dilution_conc = case_when(sample_dilution == "dil1" ~5000.000,
                                          sample_dilution == "dil2" ~2500.000,
                                          sample_dilution == "dil3" ~1250.000,
                                          sample_dilution == "dil4" ~625.000,
                                          sample_dilution == "dil5" ~312.500,
                                          sample_dilution == "dil6" ~156.250,
                                          sample_dilution == "dil7" ~78.125,
                                          sample_dilution == "dil8" ~39.062,
                                          sample_dilution == "dil9" ~19.531,
                                          sample_dilution == "dil10" ~9.766,
                                          sample_dilution == "dil11" ~4.883,
                                          TRUE ~ NA)
  )

#fig3_plot 
custom_labels <- c("S1" = "A.", 
                   "S2" = "B.", 
                   "S3" = "C.", 
                   "S4" = "D.")

pout1 <- prob_barplot(
  dfl = fig3_comb_s %>% filter(sample %in% paste0("S",1:4)), 
  xval = "sample_dilution_conc", 
  xlabtext = "Plasmid Concentration (copies/ul)", 
  yval = "mean_frag_prob", 
  ylabtext = "Probability estimate", 
  fillval = "frag_type", 
  filllabtext = "Fragment type", 
  showlabs = FALSE, 
  labval = "max_frag_val",
  facetval = "sample", 
  x_breaks = sapply(c(5000.000,2500.000,1250.000,625.000,312.500,156.250,78.125,39.062,19.531,9.766,4.883), round, 0)
) + 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.0)) + 
  theme(axis.text.x = element_text(size = rel(0.9),
                                   angle = 90,
                                   vjust = 0.5,
                                   hjust = 1)) +
  facet_wrap(~sample, labeller = labeller(sample = custom_labels), nrow = 1)


custom_labels2 <- c("S5" = "E.", 
                    "S6"= "F.", 
                    "S7"= "G. ", 
                    "S8"= "H.")

pout2 <- prob_barplot(dfl = fig3_comb_s %>% filter(sample %in% paste0("S",5:8)), 
                      xval = "sample_dilution_conc", 
                      xlabtext = "Plasmid Concentration (copies/ul)", 
                      yval = "mean_frag_prob", 
                      ylabtext = "Probability estimate", 
                      fillval = "frag_type", 
                      filllabtext = "Fragment type", 
                      showlabs=FALSE, #hiding max FPD when combining replicates
                      labval = "max_frag_val",
                      facetval = "sample", 
                      x_breaks = sapply(c(5000.000,2500.000,1250.000,625.000,312.500,156.250,78.125,39.062,19.531,9.766,4.883),round, 0)
) + theme(axis.text.x = element_text(size=rel(0.9),
                                     angle=90,
                                     vjust=0.5,
                                     hjust=1) )+
  facet_wrap(~sample, labeller = labeller(sample = custom_labels2), nrow = 1)


#Modify layout
mod <- ggpubr::ggarrange(pout1 + rremove("ylab") + rremove("xlab"), 
                         pout2 + rremove("ylab"), 
                         nrow=2, 
                         common.legend = TRUE, 
                         legend = "right" ) 

mod_annot <- annotate_figure(mod, left=textGrob("Probability estimate", rot=90, vjust=1, gp=gpar(cex=1.1) ) ) 


#Save vector and hires raster
ggsave(plot = mod_annot, filename="fig3.svg", path = here::here("manuscript_figures"), create.dir = TRUE, width = 10, height=6)

ggsave(plot = mod_annot, filename="fig3.png", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8, height=6, 
       dpi=600,
       bg="white")

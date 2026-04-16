library(here)
library(grid)

#load functions
source(here::here("000_qxone_cluster_polysolver_v3.R"))

#set seed
set.seed(2025)

#Load samples
fig4_rep1 <- run_silently(analysis_wrapper(analysis_workbook = here::here("plasmid_data.xlsx"),
                                           raw_sheet = "fig_4_rep1_raw",
                                           cluster_sheet = "fig_4_rep1_cluster",
                                           probe_order = "VIC_Cy5_FAM_Cy5.5") )

fig4_rep2 <- run_silently(analysis_wrapper(analysis_workbook = here::here("plasmid_data.xlsx"),
                                           raw_sheet = "fig_4_rep2_raw",
                                           cluster_sheet = "fig_4_rep2_cluster",
                                           probe_order = "VIC_Cy5_FAM_Cy5.5") )

fig4_rep3 <- run_silently(analysis_wrapper(analysis_workbook = here::here("plasmid_data.xlsx"),
                                           raw_sheet = "fig_4_rep3_raw",
                                           cluster_sheet = "fig_4_rep3_cluster",
                                           probe_order = "VIC_Cy5_FAM_Cy5.5") )

fwritep(x = fig4_rep1, file = here::here("tmp_files", "fig4_rep1.csv"))
fwritep(x = fig4_rep2, file = here::here("tmp_files", "fig4_rep2.csv"))
fwritep(x = fig4_rep3, file = here::here("tmp_files", "fig4_rep3.csv"))

#Format the data to long for plot
linkage_metadata_df <- primer_order_picker("VIC_Cy5_FAM_Cy5.5")
linked <- linkage_metadata_df %>% filter(experimental_linkage == "linked" | droplet_type == "singlet") %>%
  pull(fragment_combination) %>% gsub("channel","prob",.)

#Format for rep1
fig4_rep1_l <- fig4_rep1 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         well_row = gsub("[[:digit:]]","",well),
         well_row_num = match(well_row, LETTERS),
         sample_dilution = ifelse(well_row_num != 8, paste0("dil",well_row_num), NA) ) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 

#Format for rep2 
fig4_rep2_l <- fig4_rep2 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         well_row = gsub("[[:digit:]]","",well),
         well_row_num = match(well_row, LETTERS),
         sample_dilution = ifelse(well_row_num != 8, paste0("dil",well_row_num), NA) ) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 

#Format for rep3
fig4_rep3_l <- fig4_rep3 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         well_row = gsub("[[:digit:]]","",well),
         well_row_num = match(well_row, LETTERS),
         sample_dilution = ifelse(well_row_num != 8, paste0("dil",well_row_num), NA) ) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 

#Combine
fig4_comb <- rbindlist(list(rep1=fig4_rep1_l,
                            rep2=fig4_rep2_l,
                            rep3=fig4_rep3_l),idcol="rep")

#Void wells where accepted_droplets is <10k
fig4_comb <- fig4_comb %>% 
  mutate(across(.cols = starts_with("prob_"),
                .fns = ~ if_else(accepted_droplets < 10000,NA,.x) ) )

fwritep(x = fig4_comb, file = here::here("tmp_files", "fig4_rep_comb.csv"))

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
                                 levels = paste0("dil",1:7) ),
         #Sample dilutions from original data prep
         sample_dilution_conc = case_when(sample_dilution == "dil1" ~5000.000,
                                          sample_dilution == "dil2" ~1785.714,
                                          sample_dilution == "dil3" ~637.755,
                                          sample_dilution == "dil4" ~227.770,
                                          sample_dilution == "dil5" ~81.346,
                                          sample_dilution == "dil6" ~29.052,
                                          sample_dilution == "dil7" ~10.37,
                                          TRUE ~ NA),
         sample = factor(sample, levels=paste0("S",1:12))
  )

custom_labels1 <- c("S1" = "A. ", 
                   "S2" = "B. ", 
                   "S3" = "C. ", 
                   "S4" = "D. ",
                   "S5" = "E. ", 
                   "S6" = "F. ")

#fig3_plot 
pout1 <- prob_barplot(dfl = fig4_comb_s %>% filter(sample %in% paste0("S",1:6)), 
                      xval = "sample_dilution_conc", 
                      xlabtext = "Plasmid Concentration (copies/ul)", 
                      yval = "mean_frag_prob", 
                      ylabtext = "Probability estimate", 
                      fillval = "frag_type", 
                      filllabtext = "Fragment type", 
                      showlabs=FALSE, #hiding max FPD when combining replicates
                      labval = "max_frag_val",
                      facetval = "sample", 
                      x_breaks = sapply(c(5000.000, 1785.714, 637.755, 227.770, 81.346, 29.052, 10.37),round, 0),
                      facet_col = 6
) + 
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0)) + 
  facet_wrap(~sample, labeller = labeller(sample = custom_labels1), nrow = 1)+
  theme(axis.text.x = element_text(size=rel(0.9),
                                   angle=90,
                                   vjust=0.5,
                                   hjust=1) )

custom_labels2 <- c("S7" = "G. ", 
                    "S8" = "H. ", 
                    "S9" = "I. ", 
                    "S10" = "J. ",
                    "S11" = "K. ", 
                    "S12" = "L. ")
pout2 <- prob_barplot(dfl = fig4_comb_s %>% filter(sample %in% paste0("S",7:12)), 
                      xval = "sample_dilution_conc", 
                      xlabtext = "Plasmid Concentration (copies/ul)", 
                      yval = "mean_frag_prob", 
                      ylabtext = "Probability estimate", 
                      fillval = "frag_type", 
                      filllabtext = "Fragment type", 
                      showlabs=FALSE, #hiding max FPD when combining replicates
                      labval = "max_frag_val",
                      facetval = "sample", 
                      x_breaks = sapply(c(5000.000, 1785.714, 637.755, 227.770, 81.346, 29.052, 10.37),round, 0),
                      facet_col = 6
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
ggsave(plot = mod_annot, filename="fig4_viz.svg", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8.5, height=4)

ggsave(plot = mod_annot, filename="fig4_viz.png", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8.5, height=4, 
       dpi=600,
       bg="white")

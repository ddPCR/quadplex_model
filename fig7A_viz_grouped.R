library(here)
library(grid)

#load functions
source(here::here("000_qxone_cluster_polysolver_v3.R"))

#set seed
set.seed(2025)

#Load samples
fig7_rep1 <- run_silently(analysis_wrapper(analysis_workbook = here::here("aav_data_5Kth.xlsx"),
                                           raw_sheet = "fig_7_14JUL25_raw",
                                           cluster_sheet = "fig_7_14JUL25_cluster",
                                           probe_order = "VIC_Cy5_Cy5.5_FAM") )

fig7_rep2 <- run_silently(analysis_wrapper(analysis_workbook = here::here("aav_data_5Kth.xlsx"),
                                           raw_sheet = "fig_7_04AUG25_raw",
                                           cluster_sheet = "fig_7_04AUG25_cluster",
                                           probe_order = "VIC_Cy5_Cy5.5_FAM") )

fig7_rep3 <- run_silently(analysis_wrapper(analysis_workbook = here::here("aav_data_5Kth.xlsx"),
                                           raw_sheet = "fig_7_06AUG25_raw",
                                           cluster_sheet = "fig_7_06AUG25_cluster",
                                           probe_order = "VIC_Cy5_Cy5.5_FAM") )

fwritep(x = fig7_rep1, file = here::here("tmp_files", "fig7_rep1.csv"))
fwritep(x = fig7_rep2, file = here::here("tmp_files", "fig7_rep2.csv"))
fwritep(x = fig7_rep3, file = here::here("tmp_files", "fig7_rep3.csv"))


#Format the data to long for plot
linkage_metadata_df <- primer_order_picker("VIC_Cy5_Cy5.5_FAM")
linked <- linkage_metadata_df %>% filter(experimental_linkage == "linked" | droplet_type == "singlet") %>%
  pull(fragment_combination) %>% gsub("channel","prob",.)

#Format for rep1
fig7_rep1_l <- fig7_rep1 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         well_row = gsub("[[:digit:]]","",well),
         well_row_num = match(well_row, LETTERS),
         sample = paste0(sample,";", sample_dilution), #keep treatment information
         sample_dilution = ifelse(!well_row_num %in% c(8), paste0("dil",well_row_num), NA) ) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 

#Format for rep2
fig7_rep2_l <- fig7_rep2 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         well_row = gsub("[[:digit:]]","",well),
         well_row_num = match(well_row, LETTERS),
         sample = paste0(sample,";", sample_dilution), #keep treatment information
         sample_dilution = ifelse(!well_row_num %in% c(8), paste0("dil",well_row_num), NA) ) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 

#Format for rep3 
fig7_rep3_l <- fig7_rep3 %>% 
  mutate(well_col = gsub("[[:alpha:]]","",well),
         well_col_num = as.numeric(well_col),
         well_row = gsub("[[:digit:]]","",well),
         well_row_num = match(well_row, LETTERS),
         sample = paste0(sample,";", sample_dilution), #keep treatment information
         sample_dilution = ifelse(!well_row_num %in% c(8), paste0("dil",well_row_num), NA) ) %>%
  select(well, sample, sample_dilution, all_of(linked),
         prob_1111 = est_across_interest, 
         intact_conc_copies_ul_estimate = theoretical_conc_copies_ul_cluster_channel_1111,
         max_frag_val,
         well_col,
         well_col_num,
         accepted_droplets) %>%
  arrange(well_col_num) 


#Combine
fig7_comb <- rbindlist(list(rep1=fig7_rep1_l,
                            rep2=fig7_rep2_l,
                            rep3=fig7_rep3_l),idcol="rep")


#Void wells where accepted_droplets is <10k
fig7_comb <- fig7_comb %>% 
  mutate(across(.cols = starts_with("prob_"),
                .fns = ~ if_else(accepted_droplets < 10000,NA,.x) ) )

#Get mean value and nonmissing n
fig7_comb_s <- fig7_comb %>% 
  pivot_longer(cols=starts_with("prob_"), names_to="frag_type", values_to="frag_prob") %>% 
  group_by(sample, frag_type) %>% 
  summarize(mean_frag_prob = if(all(is.na(frag_prob))) NA_real_ else mean(frag_prob, na.rm=T),
            non_miss_frag_prob = sum(!is.na(frag_prob))
  ) %>% 
  ungroup() 

#Get max FPD and join
fig7_fpd <- fig7_comb %>% 
  select(sample, max_frag_val) %>% 
  group_by(sample) %>% 
  summarize(max_frag_val = round(mean(max_frag_val, na.rm=T),0) ) %>% 
  ungroup()


fig7_comb_s_format <- fig7_comb_s %>% 
  left_join(fig7_fpd, by=c("sample")) %>% 
  filter(!sample %in% c("Buffer;NA","blank;NA") ) %>% 
  mutate(sample_type_extract = str_extract(sample, "\\;.*$"),
         sample_type = gsub(";","",sample_type_extract),
         
         sample_type = ifelse(sample_type=="NA", "Plasmid", sample_type),
         
         sample_extract = str_extract(sample, ".*\\;"),
         sample = gsub("\\;","",sample_extract),
         sample = ifelse(sample == "pDDK003_MfeI_04MAR25", "pDDK003", sample))


prob_barplot1 <- function(dfl,  xval,   yval,  ylabtext, fillval, filllabtext, showlabs = TRUE, labval = NULL
){  
  my_pal <- c(
    "prob_1000" = "#1B9E77",
    "prob_0100" = "#93752C",
    "prob_0010" = "#BD6332",
    "prob_0001" = "#7E6EA2",
    "prob_1001" = "#B3499C",
    "prob_0110" = "#CF3F76",
    "prob_0011" = "#7D8F31",
    "prob_1010" = "#A0A811",
    "prob_1110" = "#E0A604",
    "prob_1011" = "#B78415",
    "prob_0111" = "#8E7037",
    "prob_1111" = "#666666"
  )
  
  s1_test <- ggplot(dfl) +
    geom_bar(aes(x = .data[[xval]],
                 y = .data[[yval]],
                 fill = .data[[fillval]]),
             position = "stack",
             stat = "identity", na.rm=T) +
    scale_fill_manual(values = my_pal) +
    labs(y=ylabtext,
         x="",
         fill=filllabtext) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(
            size=12, 
            face="bold",
            hjust=0,
            margin=margin(t=0, r=0, b=5, l=0, unit="pt")
          ),
          panel.spacing.y = unit(1,"lines")
    )
  
  if(showlabs){
    s1_test <- s1_test + geom_label(aes(x=.data[[xval]], y=1.25, label=.data[[labval]] ), 
                                    hjust=0.5, 
                                    show.legend = F, 
                                    fill="white") 
    
  }
  
  return(s1_test)
  
}




pout1_threshold <- prob_barplot1(
  dfl = fig7_comb_s_format %>%
    mutate(sample_type = as.character(sample_type),  # convert factor to character
           sample_type = ifelse(sample_type == "0 min UV", "control", sample_type),
           sample_type = factor(sample_type,
                                levels = c("Plasmid", "control","15 min UV", "30 min UV", "60 min UV"))),
  xval = "sample_type",
  yval = "mean_frag_prob", 
  ylabtext = "Probability estimate", 
  fillval = "frag_type", 
  filllabtext = "Fragment type", 
  showlabs = FALSE, 
  labval = "max_frag_val"
)

print(pout1_threshold)

ggsave(plot = pout1_threshold, filename="fig7_threshold_corrected.svg", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8.5, height=4)

ggsave(plot = pout1_threshold, filename="fig7_threshold_corrected.png", path = here::here("manuscript_figures"), create.dir = TRUE, width = 8.5, height=4, 
       dpi=600,
       bg="white")

library(here)
library(grid)
library(scales)
library(broom)

#load functions
source(here::here("000_qxone_cluster_polysolver_v3.R"))

#set seed
set.seed(2025)

#Load tmp files
#fig3
fig3_rep1 <- fread(file = here::here("tmp_files", "fig3_rep1.csv")) %>% mutate(rep="rep1")
fig3_rep2 <- fread(file = here::here("tmp_files", "fig3_rep2.csv")) %>% mutate(rep="rep2")
fig3_rep3 <- fread(file = here::here("tmp_files", "fig3_rep3.csv")) %>% mutate(rep="rep3")
fig3_rep4 <- fread(file = here::here("tmp_files", "fig3_rep4.csv")) %>% mutate(rep="rep4")
fig3_df <- rbindlist(list(fig3_rep1,fig3_rep2,fig3_rep3,fig3_rep4), use.names=TRUE)

#fig4
fig4_rep1 <- fread(file = here::here("tmp_files", "fig4_rep1.csv")) %>% mutate(rep="rep1")
fig4_rep2 <- fread(file = here::here("tmp_files", "fig4_rep2.csv")) %>% mutate(rep="rep2")
fig4_rep3 <- fread(file = here::here("tmp_files", "fig4_rep3.csv")) %>% mutate(rep="rep3")
fig4_df <- rbindlist(list(fig4_rep1,fig4_rep2,fig4_rep3), use.names=TRUE)

#formatted fig3
fig3_comb <- fread(file = here::here("tmp_files", "fig3_rep_comb.csv"))
#formatted fig4
fig4_comb <- fread(file = here::here("tmp_files", "fig4_rep_comb.csv"))

#Prep data joining empty droplets
#drop problem wells based on result data
fig3_comb_empty <- fig3_comb %>% 
  left_join(fig3_df %>% select(rep,well,count_channel_0000), by=c("rep","well") ) %>%
  filter(accepted_droplets > 10000)
fig4_comb_empty <- fig4_comb %>%
  left_join(fig4_df %>% select(rep,well,count_channel_0000), by=c("rep","well") ) %>% 
  filter(accepted_droplets > 10000)

#Isolate expected intact and degraded samples from fig3 experiment
fig3_comb_type <- fig3_comb_empty %>% 
  mutate(scenario = case_when(sample == "S1" ~ "intact",
                              sample == "S8" ~ "deg",
                              TRUE ~ NA_character_) ) %>% 
  filter(scenario %in% c("intact","deg")) %>% 
  mutate(expected_final_conc_copies_u_l = case_when(sample_dilution == "dil1" ~5000.000,
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
                                                    TRUE ~ NA))


fig4_comb_type <- fig4_comb_empty %>%
  mutate(scenario = case_when(sample %in% c("S6","S12") ~ "intact",
                              sample == "S10" ~ "d50",
                              sample == "S8" ~ "d14",
                              sample == "S7" ~ "deg",
                              TRUE ~ NA_character_) ) %>% 
  filter(scenario %in% c("intact","deg","d50","d14")) %>%
  mutate(expected_final_conc_copies_u_l = case_when(sample_dilution == "dil1" ~5000.000,
                                                    sample_dilution == "dil2" ~1785.714,
                                                    sample_dilution == "dil3" ~637.755,
                                                    sample_dilution == "dil4" ~227.770,
                                                    sample_dilution == "dil5" ~81.346,
                                                    sample_dilution == "dil6" ~29.052,
                                                    sample_dilution == "dil7" ~10.37,
                                                    TRUE ~ NA))

#Summarize expected intact across experiments
#assign a placeholder -1 when intact concentration (copies/ul) could not be estimated
summ_df <- rbindlist(list(fig3=fig3_comb_type,
                          fig4=fig4_comb_type), idcol="experiment") %>%
  mutate(intact_conc_copies_ul_estimate = if_else(is.na(intact_conc_copies_ul_estimate), -1, intact_conc_copies_ul_estimate)) %>% 
  drop_na(prob_1111) %>%
  mutate(experiment_label = case_when(experiment == "fig3" ~ "Plasmid RE Digest",
                                      experiment == "fig4" ~ "Plasmid Spike-in"))

######## LOESS on degraded samples ######

#Format deg only
summ_df_bound_deg <- summ_df %>% 
  filter(scenario == "deg") %>% 
  filter(count_channel_0000 > 0) %>%
  mutate(count_channel_0000 = as.numeric(count_channel_0000)) %>% 
  as.data.frame() %>%
  mutate(prob_1111 = -1 * prob_1111, #invert axis to make elbow shape when estimates fall below zero
         prob_1111 = rescale(prob_1111, to=c(0,1))) #scale to 1

model = loess(prob_1111 ~ count_channel_0000, data = summ_df_bound_deg, span=0.75)

summ_df_bound_deg <- augment(model, summ_df_bound_deg)

summ_df_bound_deg <- summ_df_bound_deg %>% arrange(count_channel_0000, prob_1111)

summ_df_bound_deg$slope <- c(NA, diff(summ_df_bound_deg$.fitted)/diff(summ_df_bound_deg$count_channel_0000))
summ_df_bound_deg$slope_dx <- c(NA, round(diff(summ_df_bound_deg$slope), 5)) #Use lower precision to get range of points
summ_df_bound_deg$change <- ifelse(summ_df_bound_deg$slope_dx == max(summ_df_bound_deg$slope_dx, na.rm = T), "change", "")

#Get change points
elbow_points <- summ_df_bound_deg %>% filter(change == "change") 

empty_droplets_plot_a <- summ_df_bound_deg %>%
  ggplot(aes(x=count_channel_0000,
             y=prob_1111,
             fill= log2(expected_final_conc_copies_u_l) ) ) +
  geom_point(shape=21) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       midpoint = median(log(summ_df$expected_final_conc_copies_u_l) ) ) +
  scale_y_continuous(labels = scales::percent) +
  geom_line(aes(count_channel_0000,.fitted), color="#fc8d59",  linetype=2) +
  geom_vline(xintercept = elbow_points$count_channel_0000, linetype=2, color="#a1d76a") +
  annotate("text", x = elbow_points$count_channel_0000[1], y = 0.75, label ="\nElbow point", size = 4, angle=90) +
  theme_bw() +
  labs(fill = "Expected Conc. \n (copies/ul) \n (log2-scaled)",
       y = "% Intact (Scaled Estimate)",
       x = "Empty Droplet Count")

edpa_leg <- get_legend(empty_droplets_plot_a + theme(legend.key.size = unit(0.4, "cm"),
                                                     legend.text = element_text(size=rel(0.75)),
                                                     legend.title = element_text(size=rel(0.75)) )
)

edpa_noleg <- empty_droplets_plot_a + rremove("xlab") + #ggtitle("A") + 
  theme(legend.position="none")
#edpb_noleg <- empty_droplets_plot_b + rremove("ylab") + rremove("xlab") + ggtitle("B") + theme(legend.position="none")
#edpc_noleg <- empty_droplets_plot_c + rremove("ylab") + rremove("xlab") + ggtitle("C") + theme(legend.position="none")

mod <- ggpubr::ggarrange(edpa_noleg,
                         #edpb_noleg,
                         #edpc_noleg,
                         nrow=1) 

mod_annot <- annotate_figure(mod, bottom=textGrob("Empty Droplet Count", gp=gpar(cex=1.0) ) ) 

mod_leg <- ggpubr::ggarrange(edpa_leg,
                             #edpc_leg,
                             ncol=1)

mod_comb <- ggpubr::ggarrange(mod_annot,
                              mod_leg,
                              ncol=2,
                              widths=c(4,1),
                              align="hv")

ggsave(plot = mod_comb, filename="fig5_viz_a_only.svg", path = here::here("manuscript_figures"), create.dir = TRUE, width = 5, height=4)

ggsave(plot = mod_comb, filename="fig5_viz_a_only.png", path = here::here("manuscript_figures"), create.dir = TRUE, width = 5, height=4,
       dpi=600,
       bg="white")

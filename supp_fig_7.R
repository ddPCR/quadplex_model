library(tidyverse)

expected_probs <- list(
  S1 = c(prob_1111 = 0, prob_0110 = 0.5, prob_1001=0.5),
  S2 = c(prob_1111 = 0.25, prob_0110 = 0.375, prob_1001=0.375),
  S3 = c(prob_1111 = 0.428571429, prob_0110 = 0.285714286, prob_1001=0.285714286),
  S4 = c(prob_1111 = 0.666666667, prob_0110 = 0.166666667, prob_1001=0.166666667),
  S5 = c(prob_1111 = 0.818181818,prob_0110 = 0.090909091, prob_1001=0.090909091),
  S6 = c(prob_1111 = 1,prob_0110 = 0, prob_1001=0),
  S7 = c(prob_1111 = 0, prob_0100 = 0.25, prob_0010 = 0.25,prob_1000=0.25,prob_0001=0.25),
  S8 = c(prob_1111 = 0.142857143, prob_0100 = 0.214285714, prob_0010 = 0.214285714,prob_1000=0.214285714,prob_0001=0.214285714),
  S9 = c(prob_1111 = 0.272727273, prob_0100 = 0.181818182, prob_0010 = 0.181818182,prob_1000=0.181818182,prob_0001=0.181818182),
  S10 = c(prob_1111 = 0.5, prob_0100 = 0.125, prob_0010 = 0.125,prob_1000=0.125,prob_0001=0.125),
  S11 = c(prob_1111 = 0.692307692, prob_0100 = 0.076923077, prob_0010 = 0.076923077,prob_1000=0.076923077,prob_0001=0.076923077),
  S12 = c(prob_1111 = 1, prob_0100 = 0, prob_0010 = 0,prob_1000=0,prob_0001=0)
)



fig4_comb_filtered <- fig4_comb %>%
  filter(accepted_droplets >= 10000)


nrow(fig4_comb) - nrow(fig4_comb_filtered)


fragment_recovery_long <- pmap_dfr(fig4_comb_filtered, function(...) {
  row <- list(...)
  sample_id <- row$sample
  well_id <- row$well
  rep_id <- row$rep
  sample_dilution_id <- row$sample_dilution
  max_frag_val_id <- row$max_frag_val
  intact_conc_copies_ul_estimate_id <- row$intact_conc_copies_ul_estimate 
  expected <- expected_probs[[sample_id]]
  
  # skip rows without expected fragments
  if (is.null(expected)) return(NULL)
  
  # observed values for expected fragments
  observed <- sapply(names(expected), function(col) as.numeric(row[[col]]))
  
  # calculate percent recovery per fragment
  percent_recovery_per_fragment <- (observed / expected) * 100
  
  tibble(
    rep = rep_id,
    well = well_id,
    sample = sample_id,
    sample_dilution = sample_dilution_id,
    max_frag_val = max_frag_val_id,
    intact_conc_copies_ul_estimate = intact_conc_copies_ul_estimate_id,
    fragment = names(expected),
    observed = observed,
    expected = unname(expected),
    percent_recovery = percent_recovery_per_fragment
  )
})


fragment_recovery_long <- fragment_recovery_long %>%
  mutate(
    sample_dilution = factor(sample_dilution,
                             levels = unique(sample_dilution)[
                               order(as.numeric(str_extract(unique(sample_dilution), "\\d+")))
                             ])
  )


fragment_recovery_long_filtered <- fragment_recovery_long %>%
  filter(!expected == 0)

summary_stats <- fragment_recovery_long_filtered %>%
  group_by(sample
           #,fragment
  ) %>%
  summarise(
    mean_percent_recovery = mean(percent_recovery, na.rm = TRUE),
    sd_percent_recovery = sd(percent_recovery, na.rm = TRUE),
    rsd_percent_recovery = ifelse(mean_percent_recovery != 0,
                                  (sd_percent_recovery / mean_percent_recovery) * 100,
                                  NA_real_)
  ) %>%
  ungroup()


dilution_vals <- list(
  dil1 = "5000",
  dil2 = "1786",
  dil3 = "638",
  dil4 = "228",
  dil5 = "81",
  dil6 = "29",
  dil7 = "10"
)
# convert list → named character vector
dilution_labels <- unlist(dilution_vals)

fragment_recovery_long_filtered$conc <-
  dilution_vals[fragment_recovery_long_filtered$sample_dilution] |> unlist()

fragment_recovery_long_filtered <- fragment_recovery_long_filtered %>%
  mutate(
    conc = factor(
      conc,
      levels = as.numeric(unlist(dilution_vals))  # 5000 → 5
    )
  )


summary_stats <- summary_stats %>%
  mutate(sample = factor(sample, levels = paste0("S", 1:12)))

ggplot() +
  geom_col(data = summary_stats, aes(x = sample, y = mean_percent_recovery), fill = "skyblue", alpha = 0.7) +
  geom_errorbar(data = summary_stats, aes(x = sample, 
                                          ymin = mean_percent_recovery - sd_percent_recovery,
                                          ymax = mean_percent_recovery + sd_percent_recovery),
                width = 0.2, color = "blue") +
  geom_jitter(data = fragment_recovery_long_filtered, 
              aes(x = sample, y = percent_recovery, color = conc),  
              width = 0.2, size = 2, alpha = 0.8) +
 # scale_y_continuous(limits = c(-0.50, 130))+
  geom_hline(
    yintercept = 100,
    linetype = "dotted",
    color = "grey70",
    linewidth = 0.8
  ) +
  labs(
    title = "Percent Recovery per Sample",
    x = "Sample",
    y = "Percent Recovery (%)",
    color = "Sample concentration copies/ul"
  ) +
  theme_minimal()


ggsave(
  filename = "fig4_supp1.png",
  plot = last_plot(),   
  width = 180,          
  height = 100,        
  units = "mm",
  dpi = 600,
  bg = "white"
) 


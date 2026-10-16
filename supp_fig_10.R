library(tidyverse)

expected_probs <- list(
  S1 = c(prob_1111 = 1),
  S2 = c(prob_0100 = 0.5, prob_1011 = 0.5),
  S3 = c(prob_0110 = 0.5, prob_1001 = 0.5),
  S4 = c(prob_1110 = 0.5, prob_0001 = 0.5),
  S5 = c(prob_0100 = (1/3), prob_0010 = (1/3), prob_1001 = (1/3)),
  S6 = c(prob_0100 = (1/3), prob_1010 = (1/3), prob_0001 = (1/3)),
  S7 = c(prob_0110 = (1/3), prob_1000 = (1/3), prob_0001 = (1/3)),
  S8 = c(prob_0100 = 0.25, prob_0010 = 0.25, prob_0001 = 0.25, prob_1000 = 0.25)
)


fig3_comb_filtered <- fig3_comb %>%
  filter(accepted_droplets >= 10000)


nrow(fig3_comb) - nrow(fig3_comb_filtered)


fragment_recovery_long <- pmap_dfr(fig3_comb_filtered, function(...) {
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
  
  # build tibble in long format
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

summary_stats <- fragment_recovery_long %>%
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
  dil2 = "2500",
  dil3 = "1250",
  dil4 = "625",
  dil5 = "312",
  dil6 = "156",
  dil7 = "78",
  dil8 = "39",
  dil9 = "20",
  dil10 = "10",
  dil11 = "5"
)
# convert list → named character vector
dilution_labels <- unlist(dilution_vals)

fragment_recovery_long$conc <-
  dilution_vals[fragment_recovery_long$sample_dilution] |> unlist()

fragment_recovery_long <- fragment_recovery_long %>%
  mutate(
    conc = factor(
      conc,
      levels = as.numeric(unlist(dilution_vals)) 
    )
  )


ggplot() +
  geom_col(data = summary_stats, aes(x = sample, y = mean_percent_recovery), fill = "skyblue", alpha = 0.7) +
  geom_errorbar(data = summary_stats, aes(x = sample, 
                                          ymin = mean_percent_recovery - sd_percent_recovery,
                                          ymax = mean_percent_recovery + sd_percent_recovery),
                width = 0.2, color = "blue") +
  geom_jitter(data = fragment_recovery_long, 
              aes(x = sample, y = percent_recovery, color = conc), 
              width = 0.2, size = 2, alpha = 0.8) +
  scale_y_continuous(limits = c(-0.50, 130))+
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
  filename = "fig3_supp1.png",
  plot = last_plot(),   
  width = 180,         
  height = 100,        
  units = "mm",
  dpi = 600,
  bg = "white"
) 

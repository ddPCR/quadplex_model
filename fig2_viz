library(tidyverse)
library(readxl)
library(here)

file_path <- here("plasmid_data.xlsx")

rep1 <- read_excel(file_path, sheet = 1)
rep2 <- read_excel(file_path, sheet = 3)
pattern <- tibble(
  `Sample description 2` = rep(paste0("dil", 1:12), each = 4)
)

rep2 <- rep2 %>%
  mutate(
    `Sample description 2` = rep(pattern$`Sample description 2`, times = 8)
  )

reps <- list(rep1 = rep1, rep2 = rep2)
reps_subset <- lapply(names(reps), function(nm) {
  reps[[nm]] %>%
    filter(`Sample description 1` == "S1") %>%
    filter(`Conc(copies/µL)` != "No Call") %>%
    mutate(
      `Conc(copies/µL)` = as.numeric(`Conc(copies/µL)`),
      replicate = nm
    ) %>%
    select(
      Well,
      replicate,
      `Sample description 1`,
      `Sample description 2`,
      Target,
      `Conc(copies/µL)`,
      `DyeName(s)`
    )
})

combined_df <- bind_rows(reps_subset)

dilution_lookup <- data.frame(
  `Sample description 2` = paste0("dil", 1:11),
  Dilution = 5000 / 2^(0:10),
  check.names = FALSE
)

combined_df <- combined_df %>%
  left_join(dilution_lookup, by = "Sample description 2")


combined_df <- combined_df %>%
  mutate(Recovery = (`Conc(copies/µL)` / Dilution) * 100)

combined_df <- combined_df %>%
  mutate(
    Recovery = ifelse(
      !is.na(Dilution) & Dilution > 0,
      (`Conc(copies/µL)` / Dilution) * 100,
      NA
    ))


library(ggpubr)
ggplot(combined_df, aes(x = Dilution, y = Recovery, color = `DyeName(s)`)) +
  geom_point(size = 2, alpha = 0.6) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1.2,
    aes(group = `DyeName(s)`)
  ) +
  geom_hline(
    yintercept = 100,
    linetype = "dotted",
    color = "grey70",
    linewidth = 0.8
  ) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 200))+
  scale_color_brewer(palette = "PuOr") +
  labs(
    x = "Expected copies/ul (log scale)",
    y = "% Recovery"
  ) +
  theme_pubr(border = FALSE) +
  theme(legend.title = element_blank())

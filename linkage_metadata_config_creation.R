library(tidyverse)
library(combinat)
library(data.table)
library(here)

#Make dataframe guide
try1 <- 1:4
ptry1 <- permn(try1)
ptry1_df <- as.data.frame( do.call(rbind, ptry1) )

#Theoretical 4-channel
ptry1_df_format <- ptry1_df %>%
  mutate(lookup = paste0(V1,V2,V3,V4), 
         primer_set_1 = case_when(V1 == 1 ~ "FAM",
                                  V1 == 2 ~ "VIC",
                                  V1 == 3 ~ "Cy5",
                                  V1 == 4 ~ "Cy5.5"),
         primer_set_2 = case_when(V2 == 1 ~ "FAM",
                                  V2 == 2 ~ "VIC",
                                  V2 == 3 ~ "Cy5",
                                  V2 == 4 ~ "Cy5.5"),
         primer_set_3 = case_when(V3 == 1 ~ "FAM",
                                  V3 == 2 ~ "VIC",
                                  V3 == 3 ~ "Cy5",
                                  V3 == 4 ~ "Cy5.5"),
         primer_set_4 = case_when(V4 == 1 ~ "FAM",
                                  V4 == 2 ~ "VIC",
                                  V4 == 3 ~ "Cy5",
                                  V4 == 4 ~ "Cy5.5"))

#Primer order: 1234
linkage_metadata_df1 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                               "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                               "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                               "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                               "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                               "channel_1100","doublet","linked",NA,NA, #fragment 1, 2, 1+2
                               "channel_1010","doublet","unlinked",NA,NA, #fragment 1, 3
                               "channel_1001","doublet","unlinked",NA,NA, #fragment 1, 4
                               "channel_0110","doublet","linked",NA,NA, #fragment 2, 3, 2+3
                               "channel_0101","doublet","unlinked",NA,NA, #fragment 2, 4, 
                               "channel_0011","doublet","linked",NA,NA, #fragment 3, 4, 3+4
                               "channel_1110","triplet","linked","prob_1100","prob_0110", #fragment 1, 2, 1+2, 2+3, 3, 1+2+3
                               "channel_0111","triplet","linked","prob_0110","prob_0011", #fragment 2, 3, 2+3, 3+4, 4, 2+3+4
                               "channel_1011","triplet","semi_linked","prob_0011","prob_1010", #fragment 1, 3, 3+4, 4
                               "channel_1101","triplet","semi_linked","prob_1100","prob_0101", #fragment 1, 2, 1+2, 4
) %>%
  mutate(lookup=1234) %>%
  as.data.frame()


#Primer order 1243
linkage_metadata_df2 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                "channel_1110","triplet","semi_linked","prob_1100","prob_0110", # 1, 2, 1+2, 3
                                "channel_0111","triplet","linked","prob_0110","prob_0011", # 2, 3, 2+3, 4 3+4, 2+3+4
                                "channel_1011","triplet","semi_linked","prob_0011","prob_1010", # 1, 3, 4, 3+4
                                "channel_1101","triplet","linked","prob_1100","prob_0101", # 1, 2, 1+2, 4, 2+4, 1+2+4
) %>%
  mutate(lookup=1243) %>% 
  as.data.frame()


#Primer order 1423
linkage_metadata_df3 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                "channel_1110","triplet","semi_linked","prob_0110","prob_1100", # 1, 2, 3, 2+3
                                "channel_0111","triplet","linked","prob_0110","prob_0011", # 2, 3, 2+3, 4, 2+4, 2+3+4 
                                "channel_1011","triplet","semi_linked","prob_1001","prob_1010", # 1, 3, 4, 1+4
                                "channel_1101","triplet","linked","prob_1001","prob_0101", # 1, 2, 4, 1+4, 2+4, 1+2+4
) %>%
  mutate(lookup=1423) %>%
  as.data.frame()


#Primer order 4123
linkage_metadata_df4 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                "channel_1110","triplet","linked","prob_0110","prob_1100", # 1, 2, 3, 1+2, 2+3, 1+2+3
                                "channel_0111","triplet","semi_linked","prob_0110","prob_0011", # 2, 3, 2+3, 4
                                "channel_1011","triplet","semi_linked","prob_1001","prob_1010", # 1, 3, 4, 1+4 
                                "channel_1101","triplet","linked","prob_1001","prob_1100", # 1, 2, 1+2, 4, 1+4, 1+2+4
) %>%
  mutate(lookup=4123) %>%
  as.data.frame()


#Primer order 4132
linkage_metadata_df5 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                "channel_1110","triplet","linked","prob_0110","prob_1010", # 1, 2, 3, 1+3, 2+3, 1+2+3
                                "channel_0111","triplet","semi_linked","prob_0110","prob_0011", # 2, 3, 2+3, 4
                                "channel_1011","triplet","linked","prob_1001","prob_1010", # 1, 3, 1+3, 4, 1+4, 1+3+4
                                "channel_1101","triplet","semi_linked","prob_1001","prob_1100", # 1, 2, 4, 1+4
) %>%
  mutate(lookup=4132) %>%
  as.data.frame()


#Primer order 1432
linkage_metadata_df6 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                "channel_1110","triplet","semi_linked","prob_0110","prob_1010", # 1, 2, 3, 2+3
                                "channel_0111","triplet","linked","prob_0110","prob_0011", # 2, 3, 2+3, 4, 3+4, 2+3+4
                                "channel_1011","triplet","linked","prob_1001","prob_0011", # 1, 3, 4, 1+4, 3+4, 1+3+4
                                "channel_1101","triplet","semi_linked","prob_1001","prob_1100", # 1, 2, 4, 1+4
) %>%
  mutate(lookup=1432) %>%
  as.data.frame()


#Primer order 1342 
linkage_metadata_df7 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                "channel_1110","triplet","semi_linked","prob_1010","prob_1100", # 1, 2, 3, 1+3
                                "channel_0111","triplet","linked","prob_0101","prob_0011", # 2, 3, 4, 2+4, 3+4, 2+3+4
                                "channel_1011","triplet","linked","prob_1010","prob_0011", # 1, 3, 1+3, 4, 3+4, 1+3+4
                                "channel_1101","triplet","semi_linked","prob_0101","prob_1100", # 1, 2, 4, 2+4
) %>%
  mutate(lookup=1342) %>%
  as.data.frame()


#Primer order 1324 
linkage_metadata_df8 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                "channel_1110","triplet","linked","prob_1010","prob_0110", # 1, 2, 3, 1+2, 2+3, 1+2+3
                                "channel_0111","triplet","linked","prob_0101","prob_0110", # 2, 3, 2+3, 4, 2+4, 2+3+4
                                "channel_1011","triplet","semi_linked","prob_1010","prob_1001", # 1, 3, 1+3, 4
                                "channel_1101","triplet","semi_linked","prob_0101","prob_1100", # 1, 2, 4, 2+4 
) %>%
  mutate(lookup=1324) %>%
  as.data.frame()


#Primer order 3124
linkage_metadata_df9 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                "channel_1110","triplet","linked","prob_1010","prob_1100", # 1, 2, 1+2, 3, 1+3, 1+2+3
                                "channel_0111","triplet","semi_linked","prob_0101","prob_0011", # 2, 3, 4, 2+4
                                "channel_1011","triplet","semi_linked","prob_1010","prob_1001", # 1, 3, 4, 1+3
                                "channel_1101","triplet","linked","prob_1100","prob_0101", # 1, 2, 1+2, 4, 2+4, 1+2+4
) %>%
  mutate(lookup=3124) %>%
  as.data.frame()


#Primer order 3142 
linkage_metadata_df10 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                "channel_1110","triplet","semi_linked","prob_1010","prob_1100", # 1, 2, 3, 1+3
                                "channel_0111","triplet","semi_linked","prob_0101","prob_0011", # 2, 3, 4, 2+4
                                "channel_1011","triplet","linked","prob_1010","prob_1001", # 1, 3, 1+3, 4, 1+4, 1+3+4
                                "channel_1101","triplet","linked","prob_1001","prob_0101", # 1, 2, 4, 1+4, 2+4, 1+2+4
) %>%
  mutate(lookup=3142) %>%
  as.data.frame()


#Primer order 3412
linkage_metadata_df11 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                 "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                 "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                 "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                 "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                 "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                 "channel_1110","triplet","semi_linked","prob_1100","prob_1010", # 1, 2, 1+2, 3
                                 "channel_0111","triplet","semi_linked","prob_0011","prob_0101", # 2, 3, 4, 3+4
                                 "channel_1011","triplet","linked","prob_0011","prob_1001", # 1, 3, 4, 1+4, 3+4, 1+3+4
                                 "channel_1101","triplet","linked","prob_1001","prob_1100", # 1, 2, 1+2, 4, 1+4, 1+2+4
) %>%
  mutate(lookup=3412) %>%
  as.data.frame()


#Primer order 4312
linkage_metadata_df12 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                 "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                 "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                 "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                 "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                 "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                 "channel_1110","triplet","linked","prob_1100","prob_1010", # 1, 2, 1+2, 3, 1+3, 1+2+3
                                 "channel_0111","triplet","semi_linked","prob_0011","prob_0101", # 2, 3, 4, 3+4
                                 "channel_1011","triplet","linked","prob_0011","prob_1010", # 1, 3, 1+3, 4, 3+4, 1+3+4
                                 "channel_1101","triplet","semi_linked","prob_1100","prob_0101", # 1, 2, 1+2, 4
) %>%
  mutate(lookup=4312) %>%
  as.data.frame()


#Primer order 4321
linkage_metadata_df13 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                 "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                 "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                 "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                 "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                 "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                 "channel_1110","triplet","linked","prob_1100","prob_0110", # 1, 2, 1+2, 3, 2+3, 1+2+3
                                 "channel_0111","triplet","linked","prob_0011","prob_0110", # 2, 3, 2+3, 4, 3+4, 2+3+4
                                 "channel_1011","triplet","semi_linked","prob_0011","prob_1010", # 1, 3, 4, 3+4
                                 "channel_1101","triplet","semi_linked","prob_1100","prob_0101", # 1, 2, 1+2, 4
) %>%
  mutate(lookup=4321) %>%
  as.data.frame()


#Primer order 3421
linkage_metadata_df14 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                 "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                 "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                 "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                 "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                 "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                 "channel_1110","triplet","semi_linked","prob_1100","prob_0110", # 1, 2, 1+2, 3
                                 "channel_0111","triplet","linked","prob_0011","prob_0101", # 2, 3, 4, 2+4, 3+4, 2+3+4
                                 "channel_1011","triplet","semi_linked","prob_0011","prob_1010", # 1, 3, 4, 3+4
                                 "channel_1101","triplet","linked","prob_1100","prob_0101", # 1, 2, 1+2, 4, 2+4, 1+2+4
) %>%
  mutate(lookup=3421) %>%
  as.data.frame()


#Primer order 3241
linkage_metadata_df15 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                 "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                 "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                 "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                 "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                 "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                 "channel_1110","triplet","semi_linked","prob_0110","prob_1010", # 1, 2, 3, 2+3
                                 "channel_0111","triplet","linked","prob_0110","prob_0101", # 2, 3, 2+3, 4, 2+4, 2+3+4
                                 "channel_1011","triplet","semi_linked","prob_1001","prob_1010", # 1, 3, 4, 1+4
                                 "channel_1101","triplet","linked","prob_1001","prob_0101", # 1, 2, 4, 1+4, 2+4, 1+2+4
) %>%
  mutate(lookup=3241) %>%
  as.data.frame()


#Primer order 3214
linkage_metadata_df16 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                 "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                 "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                 "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                 "channel_0101","doublet","unlinked",NA,NA,  # 2, 4
                                 "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                 "channel_1110","triplet","linked","prob_1100","prob_0110", # 1, 2, 3, 1+2, 2+3, 1+2+3
                                 "channel_0111","triplet","semi_linked","prob_0110","prob_0101", # 2, 3, 4, 2+3
                                 "channel_1011","triplet","semi_linked","prob_1001","prob_0011", # 1, 3, 4, 1+4
                                 "channel_1101","triplet","linked","prob_1001","prob_1100", # 1, 2, 4, 1+2, 1+4, 1+2+4
) %>%
  mutate(lookup=3214) %>%
  as.data.frame()


#Primer order 2314
linkage_metadata_df17 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                               "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                               "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                               "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                               "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                               "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                               "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                               "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                               "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                               "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                               "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                               "channel_1110","triplet","linked","prob_1010","prob_0110", # 1, 2, 3, 1+3, 2+3, 1+2+3
                               "channel_0111","triplet","semi_linked","prob_0110","prob_0101", # 2, 3, 2+3, 4
                               "channel_1011","triplet","linked","prob_1010","prob_1001", # 1, 3, 1+3, 4, 1+4, 1+3+4
                               "channel_1101","triplet","semi_linked","prob_1001","prob_1100", # 1, 2, 4, 1+4
) %>%
  mutate(lookup=2314) %>%
  as.data.frame()


#Primer order 2341
linkage_metadata_df18 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                 "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                 "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                 "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                 "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                 "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                 "channel_1110","triplet","semi_linked","prob_0110","prob_1010", # 1, 2, 3, 2+3
                                 "channel_0111","triplet","linked","prob_0110","prob_0011", # 2, 3, 2+3, 4, 3+4, 2+3+4
                                 "channel_1011","triplet","linked","prob_1001","prob_0011", # 1, 3, 4, 1+4, 3+4, 1+3+4
                                 "channel_1101","triplet","semi_linked","prob_1001","prob_1100", # 1, 2, 4, 1+4
) %>%
  mutate(lookup=2341) %>%
  as.data.frame()


#Primer order 2431
linkage_metadata_df19 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                 "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                 "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                 "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                 "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                 "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                 "channel_1110","triplet","semi_linked","prob_1010","prob_1100", # 1, 2, 3, 1+3 
                                 "channel_0111","triplet","linked","prob_0101","prob_0011", # 2, 3, 4, 2+4, 3+4, 2+3+4
                                 "channel_1011","triplet","linked","prob_1010","prob_0011", # 1, 3, 1+3, 4, 3+4, 1+3+4
                                 "channel_1101","triplet","semi_linked","prob_0101","prob_1100", # 1, 2, 4, 2+4
) %>%
  mutate(lookup=2431) %>%
  as.data.frame()


#Primer order 4231
linkage_metadata_df20 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                 "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                 "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                 "channel_0110","doublet","linked",NA,NA, # 2, 3, 2+3
                                 "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                 "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                 "channel_1110","triplet","linked","prob_1010","prob_0110", # 1, 2, 3, 1+3, 2+3, 1+2+3
                                 "channel_0111","triplet","linked","prob_0110","prob_0101", # 2, 3, 2+3, 4, 2+4, 2+3+4
                                 "channel_1011","triplet","semi_linked","prob_1010","prob_1001", # 1, 3, 1+3, 4
                                 "channel_1101","triplet","semi_linked","prob_0101","prob_1100", # 1, 2, 4, 2+4
) %>%
  mutate(lookup=4231) %>%
  as.data.frame()


#Primer order 4213
linkage_metadata_df21 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                 "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                 "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                 "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                 "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                 "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                 "channel_1110","triplet","linked","prob_1010","prob_1100", # 1, 2, 1+2, 3, 1+3, 1+2+3
                                 "channel_0111","triplet","semi_linked","prob_0101","prob_0011", # 2, 3, 4, 2+4
                                 "channel_1011","triplet","semi_linked","prob_1010","prob_1001", # 1, 3, 1+3, 4
                                 "channel_1101","triplet","linked","prob_0101","prob_1100", # 1, 2, 1+2, 4, 2+4, 1+2+4
) %>%
  mutate(lookup=4213) %>%
  as.data.frame()


#Primer order 2413
linkage_metadata_df22 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","unlinked",NA,NA, # 1, 2
                                 "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                 "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                 "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                 "channel_0101","doublet","linked",NA,NA, # 2, 4, 2+4
                                 "channel_0011","doublet","unlinked",NA,NA, # 3, 4
                                 "channel_1110","triplet","semi_linked","prob_1010","prob_1100", # 1, 2, 3, 1+3
                                 "channel_0111","triplet","semi_linked","prob_0101","prob_0011", # 2, 3, 4, 2+4
                                 "channel_1011","triplet","linked","prob_1010","prob_1001", # 1, 3, 1+3, 4, 1+4, 1+3+4
                                 "channel_1101","triplet","linked","prob_0101","prob_1001", # 1, 2, 4, 1+4, 2+4, 1+2+4
) %>%
  mutate(lookup=2413) %>%
  as.data.frame()


#Primer order 2143
linkage_metadata_df23 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                 "channel_1010","doublet","unlinked",NA,NA, # 1, 3
                                 "channel_1001","doublet","linked",NA,NA, # 1, 4, 1+4
                                 "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                 "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                 "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                 "channel_1110","triplet","semi_linked","prob_1100","prob_1010", # 1, 2, 1+2, 3
                                 "channel_0111","triplet","semi_linked","prob_0011","prob_0101", # 2, 3, 4, 3+4
                                 "channel_1011","triplet","linked","prob_0011","prob_1001", # 1, 3, 4, 1+3, 1+4, 1+3+4
                                 "channel_1101","triplet","linked","prob_1100","prob_1001", # 1, 2, 1+2, 4, 1+4, 1+2+4
) %>%
  mutate(lookup=2143) %>%
  as.data.frame()


#Primer order 2134
linkage_metadata_df24 <- tribble(~fragment_combination, ~droplet_type, ~experimental_linkage, ~semi_linked_prob, ~unlinked_prob,
                                 "channel_1000","singlet","unlinked",NA,NA, #fragment 1
                                 "channel_0100","singlet","unlinked",NA,NA, #fragment 2
                                 "channel_0010","singlet","unlinked",NA,NA, #fragment 3
                                 "channel_0001","singlet","unlinked",NA,NA, #fragment 4
                                 "channel_1100","doublet","linked",NA,NA, # 1, 2, 1+2
                                 "channel_1010","doublet","linked",NA,NA, # 1, 3, 1+3
                                 "channel_1001","doublet","unlinked",NA,NA, # 1, 4
                                 "channel_0110","doublet","unlinked",NA,NA, # 2, 3
                                 "channel_0101","doublet","unlinked",NA,NA, # 2, 4
                                 "channel_0011","doublet","linked",NA,NA, # 3, 4, 3+4
                                 "channel_1110","triplet","linked","prob_1100","prob_1010", # 1, 2, 1+2, 3, 1+3, 1+2+3
                                 "channel_0111","triplet","semi_linked","prob_0011","prob_0101", # 2, 3, 4, 3+4
                                 "channel_1011","triplet","linked","prob_0011","prob_1010", # 1, 3, 1+3, 4, 3+4, 1+3+4
                                 "channel_1101","triplet","semi_linked","prob_1100","prob_1001", # 1, 2, 1+2, 4
) %>%
  mutate(lookup=2134) %>%
  as.data.frame()

#collect linkages
linkage_metadata_df <- rbindlist( (mget(ls(pattern="linkage_metadata_df"))) )

linkage_metadata_df_format <- linkage_metadata_df %>% 
  mutate(lookup = as.character(lookup)) %>%
  left_join(ptry1_df_format, by="lookup") %>%
  select(-c("V1","V2","V3","V4")) %>% 
  mutate(primer_order = paste(primer_set_1, primer_set_2, primer_set_3, primer_set_4, sep = "_")) %>%
  select(-lookup, -starts_with("primer_set"))

fwrite(linkage_metadata_df_format, file.path(here::here("000_linkage_metadata_config.csv")))

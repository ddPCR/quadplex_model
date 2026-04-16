library(tidyverse)
library(data.table)
library(readxl)
library(glue)
library(janitor)
library(polynom)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)
library(svglite)

#' primer order picker function
#'
#' @param primer_order character string primer order of format "FAM_VIC_Cy5_Cy5.5"
#'
#' @return primer linkages for cluster data combinations
#' @export FALSE
#'
#' @examples
primer_order_picker <- function(primer_order){
  
  #Load lookup
  if(!file.exists( here::here("000_linkage_metadata_config.csv") )){
    source( here::here("linkage_metadata_config_creation.R") )
    primer_lookup_df <- fread( here::here("000_linkage_metadata_config.csv") )
  } else {
    primer_lookup_df <- fread( here::here("000_linkage_metadata_config.csv") )
  }
  
  #Filter lookup
  primer_lookup_df_filt <- primer_lookup_df %>% 
    dplyr::filter(primer_order == !!primer_order) %>%
    dplyr::select(-primer_order) %>%
    mutate_if(is.character, list(~na_if(.,""))) %>%
    as.data.frame()
  
  return(primer_lookup_df_filt)
  
}

#' Linkage picking lookup function
#'
#' @param channel the combination of cluster from the linkage_df, eg channel_1000
#' @param linkage_df the lookup table of relationships between fragments and clusters 
#' @note the user would create this data based on the experiment based on how many fragments are linked together in the plate
#'
#' @return a string that is used to determine which function to run for each scenario based on linkage
#' @export FALSE
#'
#' @examples
linkage_picker <- function(channel, linkage_df){
  
  linkage_select <- linkage_df %>% 
    dplyr::filter(fragment_combination == !!channel)
  
  linkage_value <- case_when(linkage_select$droplet_type == "singlet" ~ "polysolver",
                             linkage_select$droplet_type == "doublet" & 
                               linkage_select$experimental_linkage == "linked" ~ "polysolver_double",
                             linkage_select$droplet_type == "doublet" & 
                               linkage_select$experimental_linkage == "unlinked" ~ "polysolver_double_unlinked",
                             linkage_select$droplet_type == "triplet" & 
                               linkage_select$experimental_linkage == "linked" ~ "polysolver_triple_linked",
                             linkage_select$droplet_type == "triplet" & 
                               linkage_select$experimental_linkage == "semi_linked" ~ "polysolver_triple_semi_linked",
                             linkage_select$droplet_type == "triplet" & 
                               linkage_select$experimental_linkage == "unlinked" ~ "polysolver_triple_unlinked"
  )
  
  return(linkage_value)
  
}

#' semi-linked probability picker
#'
#' @param channel which cluster combination is being looked up
#' @param linkage_df experimental metadata
#'
#' @return value of semi-linked channel probability prob_xxxx
#' @export
#'
#' @examples
prob_picker <- function(channel, linkage_df){
  
  linkage_select <- linkage_df %>% 
    dplyr::filter(fragment_combination == !!channel) %>%
    dplyr::pull("semi_linked_prob")
  
  unlinked_select <- linkage_df %>% 
    dplyr::filter(fragment_combination == !!channel) %>%
    dplyr::pull("unlinked_prob")
  
  return(list(linkage_select,unlinked_select))
  
}


#Calculate max fragment
#' max fragment function
#'
#' @param lambda calculated mean value lambda
#' @param threshold poison cdf
#' @param n_test event to test, default 50
#' @note n_test arbitrarily will test out to 50 as max FPD
#'
#' @return selected n max fragments per droplet
#' @export
#'
#' @examples
max_frag <- function(lambda, threshold, n_test=50){
  
  pdist <- ppois(seq(n_test), lambda)
  
  #pdist <- pbinom(seq(n_test), m, 1/AD)
  
  names(pdist) <- seq(n_test)
  
  pdist_cut <- pdist[pdist >= threshold]
  
  n_select <- as.numeric( names(pdist_cut)[which.min(pdist_cut)] )
  
  #Establish a minimum for multiplex of 2 FPD
  #n_select <- ifelse(n_select <= 1, 2, n_select)
  
  return(n_select)
  
}

#Using tlin polysolver global fx
#' polynomial equation solver using Poisson distribution
#'
#' @param m trials
#' @param AD accepted droplets
#' @param target detection channel 
#' @param n maximal number of fragments per droplet
#' @note modified to handle zero detected events to write probability 0
#'
#' @return root of equation for target
#' @export FALSE
#'
#' @examples
polysolver <- function(m, AD, target, n) {
  
  if(target == 0){
    root1 <- 0
  } else{
    
    error_flag <- FALSE
    
    tryCatch({
      
      p <- as.polynomial(c(-1, dbinom(seq(n), size = m, prob = 1 / AD) * AD / target)) #note tlin only solves using the binomial, not poisson
      p
      pz <- solve(p)
      good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
      root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
      stopifnot(length(root1) == 1)
      
    }, error = function(e){
      error_flag <<- TRUE
      message("* Caught an error during polysolve ", target)
    })
    
    if (error_flag) {
      root1 <- NA_real_
    }
    
  }
  return(root1)
}

#' Double amplicon polysolver Jake version
#'
#' @param m trials
#' @param AD accepted droplets
#' @param n maximal number of fragments per droplet
#' @param ch1_prob ch1 singlet probability
#' @param ch2_prob ch2 singlet probability
#' @param ch1ch2_target number of events in ch1ch2 channel
#'
#' @return root of equation for target
#' @export FALSE
#'
#' @examples
polysolver_double <- function(m, AD, n, 
                              ch1_prob, 
                              ch2_prob, 
                              ch1ch2_target){
  
  #If the cluster target is a zero, write out a zero
  if(ch1ch2_target == 0){
    root1 <- 0
  } else if( n == 1 ){
    #If there is measured double fragment target and FPD is 1 then there is no grid expansion
    root1 <- polysolver(m=m,AD=AD,target=ch1ch2_target,n=n)
    
    } else {
    
    error_flag <- FALSE
    
    tryCatch({
    
      #Set up events
      num_events <- n
      num_events_vec <- 0:num_events
      num_events_vec_list <- replicate(3, num_events_vec, simplify = FALSE) #ch1 ch2 ch1ch2 style channels
      
      #Set up event combinations
      #Filter and format event combinations to match physical limits
      #add row_id for join (not used in this version)
      double_grid <- expand.grid(num_events_vec_list) %>% 
        mutate(sum=Var1+Var2+Var3) %>%
        dplyr::filter(sum <= num_events, sum > 0) %>% 
        dplyr::rename(ch1=Var1,ch2=Var2,ch1ch2=Var3) %>% 
        mutate(zero_flag = ifelse( (ch1==0&ch1ch2==0) | (ch2==0&ch1ch2==0), TRUE, FALSE) ) %>%
        dplyr::filter(!zero_flag) %>%
        arrange(ch1ch2) %>%
        mutate(row_id = 1:nrow(.) )
      
      #Derive polynomial coefficients 
      double_grid$poly_coef <- ch1_prob^double_grid$ch1 * ch2_prob^double_grid$ch2 * dbinom(double_grid$sum, size = m, prob = 1 / AD)
      
      sum_coef <- double_grid %>%
        dplyr::group_by(ch1ch2) %>%
        dplyr::summarize(total_coef = sum(poly_coef)) %>%
        dplyr::pull(total_coef)
      
      #Correct for ch1ch2 target observed
      adjusted_sum_coef = sum_coef * AD / ch1ch2_target
      
      adjusted_sum_coef[1] <- adjusted_sum_coef[1] - 1 
      
      p <- as.polynomial(adjusted_sum_coef)
      pz <- solve(p)
      good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
      root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
      stopifnot(length(root1) == 1)
    
    }, error = function(e){
      error_flag <<- TRUE
      message("* Caught an error during polysolve ", ch1ch2_target)
    })
    
    if (error_flag) {
      root1 <- NA_real_
    }
    
  }
    
  return(root1)
}

#' Double amplicon unlinked channels polysolver 
#'
#' @param m trials
#' @param AD accepted droplets
#' @param n maximal number of fragments per droplet
#' @param ch1_prob ch1 singlet probability
#' @param ch3_prob ch3 singlet probability
#' @param ch1ch3_target number of events in ch1ch3 channel
#' @note assumes always separate fragments (ch1 FAM/ ch3 Cy5, etc)
#'
#' @return root of equation for target
#' @export FALSE
#'
#' @examples
polysolver_double_unlinked <- function(m, AD, n, 
                                       ch1_prob, 
                                       ch3_prob, 
                                       ch1ch3_target){
  
  if(ch1ch3_target %in% c(0,1) ){
    root1 <- 0 #Cluster Target detected for an unlinked double cannot be 1 or 0, would be impossible
  } else if(n <= 1){
    root1 <- 0 #FPD for an unlinked double must be greater than 1, otherwise would be impossible
  } else{
    
    error_flag <- FALSE
    
    tryCatch({
    
    #Set up events
    num_events <- n
    num_events_vec <- 0:num_events
    num_events_vec_list <- replicate(3, num_events_vec, simplify = FALSE) 
    
    #Set up event combinations
    #Filter and format event combinations to match physical limits
    #add row_id for join (not used in this version)
    double_grid <- expand.grid(num_events_vec_list) %>% 
      mutate(sum=Var1+Var2+Var3) %>%
      dplyr::filter(sum <= num_events, sum > 0) %>% 
      dplyr::rename(ch1=Var1,ch3=Var2,ch1ch3=Var3) %>% 
      mutate(zero_flag = ifelse( ch1ch3 != 0 | (ch1 == 0 | ch3 == 0 ), TRUE, FALSE) ) %>% #ch1ch3 sum is unlinked 
      dplyr::filter(!zero_flag) %>%
      arrange(ch1ch3) %>%
      mutate(row_id = 1:nrow(.) )
    
    #Derive polynomial coefficients 
    double_grid$poly_coef <- ch1_prob^double_grid$ch1 * ch3_prob^double_grid$ch3 * dbinom(double_grid$sum, size = m, prob = 1 / AD)
    
    root1 <- sum(double_grid$poly_coef) * AD
    
    }, error = function(e){
      error_flag <<- TRUE
      message("* Caught an error during polysolve ", ch1ch3_target)
      
    })
    
    if (error_flag) {
      root1 <- NA_real_
    }
    
    # sum_coef <- double_grid %>%
    #   #dplyr::group_by(ch1ch3) %>% #cannot group by combination 
    #   dplyr::group_by(sum) %>%
    #   dplyr::summarize(total_coef = sum(poly_coef)) %>%
    #   dplyr::pull(total_coef)
    # 
    # #assume no 0 order and single 1st order due to unlinked double
    # sum_coef <- c(0, 1, sum_coef)
    # 
    # #Correct for ch1ch3 target observed
    # adjusted_sum_coef = sum_coef * AD / ch1ch3_target
    # 
    # adjusted_sum_coef[1] <- adjusted_sum_coef[1] - 1 
    # 
    # p <- as.polynomial(adjusted_sum_coef)
    # pz <- solve(p)
    # good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
    # root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
    # stopifnot(length(root1) == 1)
    
  }
  
  return(root1)
  
}

#' Triple amplicon polysolver Jake version
#'
#' @param m trials
#' @param AD accepted droplets
#' @param n maximal number of fragments per droplet
#' @param ch1_prob ch1 singlet probability
#' @param ch2_prob ch2 singlet probability
#' @param ch3_prob ch3 singlet probability
#' @param ch1ch2_prob ch1ch2 doublet probability
#' @param ch2ch3_prob ch2ch3 doublet probability
#' @param ch1ch2ch3_target number of events in ch1ch2ch3 channel
#' @param linkage_type linked semi_linked unlinked type
#'
#' @return root of equation for target
#' @export FALSE
#'
#' @examples
polysolver_triple <- function(m, AD, n, 
                              ch1_prob, 
                              ch2_prob, 
                              ch3_prob, 
                              ch1ch2_prob, 
                              ch2ch3_prob, 
                              ch1ch2ch3_target,
                              linkage_type){
  
  if(ch1ch2ch3_target == 0){
    root1 = 0
  } else if(n==1 && linkage_type == "linked"){
    
    #If there is measured linked triple fragment target and FPD is 1 then there is no grid expansion
    root1 <- polysolver(m=m,AD=AD,target=ch1ch2ch3_target,n=n)
    
  } else{
  
  #Set up events
  num_events <- n
  num_events_vec <- 0:num_events
  num_events_vec_list <- replicate(6, num_events_vec, simplify = FALSE)
  
  error_flag <- FALSE
  
  tryCatch({
  
  #Set up event combinations
  #Filter and format event combinations to match physical limits
  #add row_id for join (not used in this version)
  if(linkage_type=="linked"){
    
    triple_grid <- expand.grid(num_events_vec_list) %>%
      mutate(sum=Var1+Var2+Var3+Var4+Var5+Var6) %>%
      dplyr::filter(sum <= num_events, sum > 0) %>%
      dplyr::rename(ch1=Var1,ch2=Var2,ch3=Var3,ch1ch2=Var4,ch2ch3=Var5,ch1ch2ch3=Var6) %>%
      mutate(good_flag = ifelse((ch1>0|ch1ch2>0|ch1ch2ch3>0)&(ch2>0|ch1ch2>0|ch2ch3>0|ch1ch2ch3>0)&(ch3>0|ch2ch3>0|ch1ch2ch3>0), TRUE, FALSE) ) %>%
      dplyr::filter(good_flag) %>%
      dplyr::arrange(ch1ch2ch3) %>%
      mutate(row_id = 1:nrow(.) )
    
    triple_grid$poly_coef = ch1_prob^triple_grid$ch1 * 
      ch2_prob^triple_grid$ch2 * 
      ch3_prob^triple_grid$ch3 * 
      ch1ch2_prob^triple_grid$ch1ch2 * 
      ch2ch3_prob^triple_grid$ch2ch3 * 
      dbinom(triple_grid$sum, size = m, prob = 1 / AD)
    
    sum_coef <- triple_grid %>%
      dplyr::group_by(ch1ch2ch3) %>%
      dplyr::summarize(total_coef = sum(poly_coef)) %>%
      dplyr::pull(total_coef)
    
    #Adjust coefficients
    adjusted_sum_coef = sum_coef * AD / ch1ch2ch3_target
    
    adjusted_sum_coef[1] = adjusted_sum_coef[1] - 1 
    
    p <- as.polynomial(adjusted_sum_coef)
    pz <- solve(p)
    good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
    root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
    stopifnot(length(root1) == 1)
  
  }
  
  if(linkage_type=="semi_linked"){
    
    triple_grid <- expand.grid(num_events_vec_list) %>%
      mutate(sum=Var1+Var2+Var3+Var4+Var5+Var6) %>%
      dplyr::filter(sum <= num_events, sum > 0) %>%
      dplyr::rename(ch1=Var1,ch2=Var2,ch3=Var3,ch1ch2=Var4,ch2ch3=Var5,ch1ch2ch3=Var6) %>%
      mutate(keep_flag = ifelse( 
        
        (ch1ch2==0 & ch2ch3 == 0 & ch1ch2ch3 == 0) & (ch1 > 0 & ch2 > 0 & ch3 > 0) |
          (ch2ch3 == 0 & ch1ch2ch3 == 0) & (ch1 >= 0 & ch2 >= 0 & ch3 > 0) & (ch1ch2 > 0) #updated logic 19SEPT2024
        
        ,TRUE,FALSE) ) %>%
      dplyr::filter(keep_flag) %>%
      dplyr::arrange(ch1ch2) %>%
      mutate(row_id = 1:nrow(.) )
    
    triple_grid$poly_coef <- ch1_prob^triple_grid$ch1 * 
      ch2_prob^triple_grid$ch2 * 
      ch3_prob^triple_grid$ch3 * 
      ch1ch2_prob^triple_grid$ch1ch2 *
      dbinom(triple_grid$sum, size = m, prob = 1 / AD)
    
    root1 <- sum(triple_grid$poly_coef) * AD
    
    return(root1)
    
    # sum_coef <- triple_grid %>%
    #   #dplyr::group_by(ch1ch3ch4) %>% #cannot group by combination 
    #   dplyr::group_by(sum) %>%
    #   dplyr::summarize(total_coef = sum(poly_coef)) %>%
    #   dplyr::pull(total_coef)
    # 
    # #assume no 0 order and single 1st order due to unlinked triple
    # sum_coef <- c(0, 1, sum_coef)
    # 
    # #Correct for ch1ch3 target observed
    # adjusted_sum_coef = sum_coef * AD / ch1ch2ch3_target
    # 
    # adjusted_sum_coef[1] <- adjusted_sum_coef[1] - 1 
    # 
    # p <- as.polynomial(adjusted_sum_coef)
    # pz <- solve(p)
    # good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
    # root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
    # stopifnot(length(root1) == 1)
    
  }
  
  if(linkage_type == "unlinked"){
    
    triple_grid <- expand.grid(num_events_vec_list) %>%
      mutate(sum=Var1+Var2+Var3+Var4) %>%
      dplyr::filter(sum <= num_events, sum > 0) %>%
      dplyr::rename(ch1=Var1,ch2=Var2,ch3=Var3,ch1ch2=Var4,ch2ch3=Var5,ch1ch2ch3=Var6) %>%
      mutate(keep_flag = ifelse( 
        
        (ch1ch2==0 & ch2ch3 == 0 & ch1ch2ch3 == 0) & (ch1 > 0 & ch2 > 0 & ch3 > 0) 
        
        ,TRUE,FALSE) ) %>%
      dplyr::filter(keep_flag) %>%
      dplyr::arrange(ch1) %>%
      mutate(row_id = 1:nrow(.) )
    
    triple_grid$poly_coef <- ch1_prob^triple_grid$ch1 * 
      ch2_prob^triple_grid$ch2 * 
      ch3_prob^triple_grid$ch3 * 
      ch1ch2_prob^triple_grid$ch1ch2 *
      dbinom(triple_grid$sum, size = m, prob = 1 / AD)
    
    root1 <- sum(triple_grid$poly_coef) * AD
    
    return(root1)
    
    # sum_coef <- triple_grid %>%
    #   #dplyr::group_by(ch1ch3ch4) %>% #cannot group by combination 
    #   dplyr::group_by(sum) %>%
    #   dplyr::summarize(total_coef = sum(poly_coef)) %>%
    #   dplyr::pull(total_coef)
    # 
    # #assume no 0 order and single 1st order due to unlinked triple
    # sum_coef <- c(0, 1, sum_coef)
    # 
    # #Correct for ch1ch3 target observed
    # adjusted_sum_coef = sum_coef * AD / ch1ch2ch3_target
    # 
    # adjusted_sum_coef[1] <- adjusted_sum_coef[1] - 1 
    # 
    # p <- as.polynomial(adjusted_sum_coef)
    # pz <- solve(p)
    # good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
    # root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
    # stopifnot(length(root1) == 1)
    
  }
    
  }, error = function(e){
    error_flag <<- TRUE
    message("* Caught an error during polysolve ", ch1ch2ch3_target)
    
  })
  
  if (error_flag) {
    root1 <- NA_real_
  }
  
  }
  
  return(root1)
  
}

#' Triple amplicon semi linked channels polysolver 
#'
#' @param m trials
#' @param AD accepted droplets
#' @param n maximal number of fragments per droplet
#' @param ch1_prob ch1 singlet probability
#' @param ch2_prob ch2 singlet probability
#' @param ch3_prob ch3 singlet probability
#' @param ch1ch2_prob ch1ch2 probability
#' @param ch2ch3_prob ch2ch3 probability NA
#' @param ch1ch2ch3_target number of events in ch1ch2ch3 channel
#' @param ch1ch2_prob NA
#'
#' @note assumes ch1, ch2, ch3, ch1ch2 fragments assort independently
#'
#' @return root of equation for target
#' @export FALSE
#'
#' @examples
polysolver_triple_semi_linked <- function(m, AD, n, 
                                       ch1_prob, 
                                       ch2_prob, 
                                       ch3_prob,
                                       ch1ch2_prob, 
                                       ch2ch3_prob = NA,
                                       ch1ch2ch3_target){
  
  if(ch1ch2ch3_target == 0){
    root1 = 0
  } else{
    
    #Set up events
    num_events <- n
    num_events_vec <- 0:num_events
    num_events_vec_list <- replicate(6, num_events_vec, simplify = FALSE)
    
    #Set up event combinations
    #Filter and format event combinations to match physical limits
    #add row_id for join (not used in this version)
    triple_grid <- expand.grid(num_events_vec_list) %>%
      mutate(sum=Var1+Var2+Var3+Var4+Var5+Var6) %>%
      dplyr::filter(sum <= num_events, sum > 0) %>%
      dplyr::rename(ch1=Var1,ch2=Var2,ch3=Var3,ch1ch2=Var4,ch2ch3=Var5,ch1ch2ch3=Var6) %>%
      mutate(keep_flag = ifelse( 
        
        (ch1ch2==0 & ch2ch3 == 0 & ch1ch2ch3 == 0) & (ch1 > 0 & ch2 > 0 & ch3 > 0) |
        (ch2ch3 == 0 & ch1ch2ch3 == 0) & (ch1 == 0 & ch2 == 0 & ch3 > 0) & (ch1ch2 > 0) 
        
        ,TRUE,FALSE) ) %>%
      dplyr::filter(keep_flag) %>%
      dplyr::arrange(ch1ch2) %>%
      mutate(row_id = 1:nrow(.) )
    
    triple_grid$poly_coef <- ch1_prob^triple_grid$ch1 * 
      ch2_prob^triple_grid$ch2 * 
      ch3_prob^triple_grid$ch3 * 
      ch1ch2_prob^triple_grid$ch1ch2 *
      dbinom(triple_grid$sum, size = m, prob = 1 / AD)
    
    sum_coef <- triple_grid %>%
      #dplyr::group_by(ch1ch3ch4) %>% #cannot group by combination 
      dplyr::group_by(sum) %>%
      dplyr::summarize(total_coef = sum(poly_coef)) %>%
      dplyr::pull(total_coef)
    
    #assume no 0 order and single 1st order due to unlinked triple
    sum_coef <- c(0, 1, sum_coef)
    
    #Correct for ch1ch3 target observed
    adjusted_sum_coef = sum_coef * AD / ch1ch2ch3_target
    
    adjusted_sum_coef[1] <- adjusted_sum_coef[1] - 1 
    
    p <- as.polynomial(adjusted_sum_coef)
    pz <- solve(p)
    good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
    root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
    stopifnot(length(root1) == 1)
    
  }
  
  return(root1)
  
}

#' Triple amplicon unlinked channels polysolver 
#'
#' @param m trials
#' @param AD accepted droplets
#' @param n maximal number of fragments per droplet
#' @param ch1_prob ch1 singlet probability
#' @param ch2_prob ch2 singlet probability
#' @param ch3_prob ch3 singlet probability
#' @param ch1ch2_prob ch1ch2 probability NA
#' @param ch2ch3_prob ch2ch3 probability NA
#' @param ch1ch2ch3_target number of events in ch1ch3ch4 channel
#'
#' @note assumes ch1, ch2, ch3 fragments assort independently
#'
#' @return root of equation for target
#' @export FALSE
#'
#' @examples
polysolver_triple_unlinked <- function(m, AD, n, 
                                       ch1_prob, 
                                       ch2_prob, 
                                       ch3_prob, 
                                       ch1ch2_prob = NA,
                                       ch2ch3_prob = NA,
                                       ch1ch2ch3_target){
  
  if(ch1ch2ch3_target == 0){
    root1 = 0
  } else{
    
    #Set up events
    num_events <- n
    num_events_vec <- 0:num_events
    num_events_vec_list <- replicate(6, num_events_vec, simplify = FALSE)
    
    #Set up event combinations
    #Filter and format event combinations to match physical limits
    #add row_id for join (not used in this version)
    triple_grid <- expand.grid(num_events_vec_list) %>%
      mutate(sum=Var1+Var2+Var3+Var4) %>%
      dplyr::filter(sum <= num_events, sum > 0) %>%
      dplyr::rename(ch1=Var1,ch2=Var2,ch3=Var3,ch1ch2=Var4,ch2ch3=Var5,ch1ch2ch3=Var6) %>%
      mutate(keep_flag = ifelse( 
        
        (ch1ch2==0 & ch2ch3 == 0 & ch1ch2ch3 == 0) & (ch1 > 0 & ch2 > 0 & ch3 > 0) 
        
        ,TRUE,FALSE) ) %>%
      dplyr::filter(keep_flag) %>%
      dplyr::arrange(ch1) %>%
      mutate(row_id = 1:nrow(.) )
    
    triple_grid$poly_coef <- ch1_prob^triple_grid$ch1 * 
      ch2_prob^triple_grid$ch2 * 
      ch3_prob^triple_grid$ch3 * 
      ch1ch2_prob^triple_grid$ch1ch2 *
      dbinom(triple_grid$sum, size = m, prob = 1 / AD)
    
    sum_coef <- triple_grid %>%
      #dplyr::group_by(ch1ch3ch4) %>% #cannot group by combination 
      dplyr::group_by(sum) %>%
      dplyr::summarize(total_coef = sum(poly_coef)) %>%
      dplyr::pull(total_coef)
    
    #assume no 0 order and single 1st order due to unlinked triple
    sum_coef <- c(0, 1, sum_coef)
    
    #Correct for ch1ch3 target observed
    adjusted_sum_coef = sum_coef * AD / ch1ch2ch3_target
    
    adjusted_sum_coef[1] <- adjusted_sum_coef[1] - 1 
    
    p <- as.polynomial(adjusted_sum_coef)
    pz <- solve(p)
    good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
    root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
    stopifnot(length(root1) == 1)
    
  }
  
  return(root1)
  
}

#' analysis wrapper
#'
#' @param analysis_workbook excel file with raw and cluster data
#' @param raw_sheet name of raw data sheet
#' @param cluster_sheet name of cluster data sheet
#' @param probe_order order of probes 
#' @note if running different probe orders on the same plate then run the subset of wells one per wrapper call
#'
#' @return a dataframe of predicted intact 
#' @export FALSE
#'
#' @examples
analysis_wrapper <- function(analysis_workbook,
                             raw_sheet, 
                             cluster_sheet,
                             probe_order){
  
  #Single Linkage pattern
  linkage_metadata_df <- primer_order_picker(probe_order)
  
  #File for experimental data
  excel_file <- file.path(analysis_workbook)
  
  #Read in xlsx sheet 1 (QXONE_rawdata)
  raw_data <- read_xlsx(excel_file, sheet=raw_sheet)
  
  #Lauren derivations
  raw_data_simple <- raw_data %>% 
    janitor::clean_names() %>%
    select(well, conc_copies_ul = conc_copies_m_l, target, dye_name_s, accepted_droplets, positives, negatives) %>% 
    mutate(
      conc_copies_ul = gsub("\\,","",conc_copies_ul),
      conc_copies_ul = as.numeric(conc_copies_ul),
      lambda_conc_ul =  -log(negatives/accepted_droplets), #this lambda is for any target+ cluster
      droplet_vol_ul = lambda_conc_ul / conc_copies_ul
    ) %>%
    pivot_wider(id_cols=c("well","accepted_droplets"), names_from = "target", values_from=c("conc_copies_ul", "lambda_conc_ul", "droplet_vol_ul", "positives", "negatives")) 
  
  raw_well_sample <- raw_data %>%
    janitor::clean_names() %>%
    select(well, sample = sample_description_1, sample_dilution = sample_description_2) %>% distinct() 
  
  #Read in xlsx sheet 2 (QXONE_clusterdata)
  cluster_data <- read_xlsx(excel_file, sheet=cluster_sheet)
  
  #droplet volume columns
  drop_vol_cols <- colnames(raw_data_simple)[ which(grepl("droplet_vol_ul", colnames(raw_data_simple) ) ) ]
  
  #droplet volume columns are the same but have different patterns of missingness depending on sample
  #choose volume column with minimal missingness
  drop_vol_df <- raw_data_simple %>% select(all_of(drop_vol_cols))
  na_drop_vol_count <- colSums(is.na(drop_vol_df))
  
  #drop_vol_col <- drop_vol_cols[1] #droplet volumes per target are similar per well
  drop_vol_col <- names(which.min(na_drop_vol_count))
  
  #Lauren derivations
  #Calculate theoretical lambda and m values
  cluster_data_simple <- cluster_data %>% 
    janitor::clean_names() %>%
    dplyr::select(well, starts_with("target"), count) %>%
    left_join(raw_data_simple %>% 
                select(well,
                       accepted_droplets,
                       all_of ( drop_vol_col )
                ), by="well") %>%
    mutate(
      lambda_conc_ul_cluster = -log( (accepted_droplets - count) / accepted_droplets),
      conc_copies_ul_cluster = lambda_conc_ul_cluster / !!sym(drop_vol_col),
      channel_barcode = paste0("channel_",target_1,target_2,target_3,target_4) 
    ) %>%
    dplyr::select(-starts_with("target")) %>%
    pivot_wider(id_cols=c("well","accepted_droplets",all_of( drop_vol_col ) ), 
                names_from = "channel_barcode", 
                values_from=c("count","lambda_conc_ul_cluster","conc_copies_ul_cluster"), 
                values_fill = 0 #Assume a missing cluster value is a true 0
    ) %>%
    mutate(
      lambda_theoretical = -log(count_channel_0000 / accepted_droplets),
      m_theoretical = floor(accepted_droplets * lambda_theoretical)
    )
  
  #Join dataframes
  if(any(cluster_data_simple$count_channel_0000 == 0)){
    
    message("No empty droplets in at least one well, substituting max lambda for that sample.")
    
    combined_df <- cluster_data_simple %>%
      dplyr::select(-accepted_droplets, -all_of(drop_vol_col)) %>%
      left_join(raw_data_simple, by="well") %>%
      as.data.frame() %>%
      left_join(raw_well_sample, by="well") %>%
      group_by(sample) %>%
      mutate(lambda_theoretical = ifelse(is.infinite(lambda_theoretical),
                                         NA_real_,
                                         lambda_theoretical),
             lambda_theoretical = ifelse(is.na(lambda_theoretical),
                                         max(lambda_theoretical,na.rm = TRUE),
                                         lambda_theoretical),
             m_theoretical = floor(accepted_droplets * lambda_theoretical)
      ) %>%
      ungroup()
    
  } else{
    
    combined_df <- cluster_data_simple %>%
      dplyr::select(-accepted_droplets, -all_of(drop_vol_col)) %>%
      left_join(raw_data_simple, by="well") %>%
      as.data.frame() %>%
      left_join(raw_well_sample, by="well") 
    
  }
  
  #Calculate maximum FPD per well
  #Run the single fragment probabilities
  combined_df_solve_single <- combined_df %>%
    
    mutate(
      
      max_frag_val = list(lambda_theoretical, 0.99) %>% purrr::pmap_dbl(.f = max_frag),
      
      prob_1000 = list(m_theoretical, accepted_droplets, count_channel_1000, max_frag_val) %>% purrr::pmap_dbl(.f = polysolver),
      prob_0100 = list(m_theoretical, accepted_droplets, count_channel_0100, max_frag_val) %>% purrr::pmap_dbl(.f = polysolver),
      prob_0010 = list(m_theoretical, accepted_droplets, count_channel_0010, max_frag_val) %>% purrr::pmap_dbl(.f = polysolver),
      prob_0001 = list(m_theoretical, accepted_droplets, count_channel_0001, max_frag_val) %>% purrr::pmap_dbl(.f = polysolver)
      
    )
  
  #Run the double fragment probabilities
  combined_df_solve_double <- combined_df_solve_single %>% 
    mutate(
      
      #Linkage pick using linkage_metadata
      prob_1100 = list(m_theoretical, accepted_droplets, max_frag_val, prob_1000, prob_0100, count_channel_1100) %>% purrr::pmap_dbl(.f = get(linkage_picker("channel_1100", linkage_df = linkage_metadata_df))),
      prob_0011 = list(m_theoretical, accepted_droplets, max_frag_val, prob_0010, prob_0001, count_channel_0011) %>% purrr::pmap_dbl(.f = get(linkage_picker("channel_0011", linkage_df = linkage_metadata_df))),
      prob_0110 = list(m_theoretical, accepted_droplets, max_frag_val, prob_0100, prob_0010, count_channel_0110) %>% purrr::pmap_dbl(.f = get(linkage_picker("channel_0110", linkage_df = linkage_metadata_df))),
      prob_1010 = list(m_theoretical, accepted_droplets, max_frag_val, prob_1000, prob_0010, count_channel_1010) %>% purrr::pmap_dbl(.f = get(linkage_picker("channel_1010", linkage_df = linkage_metadata_df))),
      prob_1001 = list(m_theoretical, accepted_droplets, max_frag_val, prob_1000, prob_0001, count_channel_1001) %>% purrr::pmap_dbl(.f = get(linkage_picker("channel_1001", linkage_df = linkage_metadata_df))),
      prob_0101 = list(m_theoretical, accepted_droplets, max_frag_val, prob_0100, prob_0001, count_channel_0101) %>% purrr::pmap_dbl(.f = get(linkage_picker("channel_0101", linkage_df = linkage_metadata_df))),
      
    )
  
  #Run the triple fragment probabilities
  combined_df_solve_triple <- combined_df_solve_double %>%
    mutate(
      
      #Linkage pick using linkage_metadata
      prob_1110 = list(m_theoretical, 
                       accepted_droplets, 
                       max_frag_val, 
                       ch1_prob = prob_1000, 
                       ch2_prob = prob_0100, 
                       ch3_prob = prob_0010, 
                       ch1ch2_prob = get( prob_picker("channel_1110", linkage_metadata_df)[[1]] ), 
                       ch2ch3_prob = get( prob_picker("channel_1110", linkage_metadata_df)[[2]] ),
                       ch1ch2ch3_target = count_channel_1110,
                       linkage_type = gsub("polysolver_triple_", "", linkage_picker("channel_1110", linkage_df = linkage_metadata_df) ) ) %>% 
        purrr::pmap_dbl(.f = polysolver_triple ),
      prob_0111 = list(m_theoretical, 
                       accepted_droplets, 
                       max_frag_val, 
                       ch1_prob = prob_0100, 
                       ch2_prob = prob_0010, 
                       ch3_prob = prob_0001, 
                       ch1ch2_prob = get( prob_picker("channel_0111", linkage_metadata_df)[[1]] ),
                       ch2ch3_prob = get( prob_picker("channel_0111", linkage_metadata_df)[[2]] ), 
                       ch1ch2ch3_target = count_channel_0111,
                       linkage_type = gsub("polysolver_triple_", "", linkage_picker("channel_0111", linkage_df = linkage_metadata_df) ) ) %>% 
        purrr::pmap_dbl(.f = polysolver_triple ),
      prob_1011 = list(m_theoretical, 
                       accepted_droplets, 
                       max_frag_val, 
                       ch1_prob = prob_1000, 
                       ch2_prob = prob_0010, 
                       ch3_prob = prob_0001,
                       ch1ch2_prob = get( prob_picker("channel_1011", linkage_metadata_df)[[1]] ), 
                       ch2ch3_prob = get( prob_picker("channel_1011", linkage_metadata_df)[[2]] ),
                       ch1ch2ch3_target = count_channel_1011,
                       linkage_type = gsub("polysolver_triple_", "", linkage_picker("channel_1011", linkage_df = linkage_metadata_df) ) ) %>% 
        purrr::pmap_dbl(.f = polysolver_triple ),
      prob_1101 = list(m_theoretical, 
                       accepted_droplets, 
                       max_frag_val, 
                       ch1_prob = prob_1000, 
                       ch2_prob = prob_0100, 
                       ch3_prob = prob_0001, 
                       ch1ch2_prob = get( prob_picker("channel_1101", linkage_metadata_df)[[1]] ),
                       ch2ch3_prob = get( prob_picker("channel_1101", linkage_metadata_df)[[2]] ), 
                       ch1ch2ch3_target = count_channel_1101,
                       linkage_type = gsub("polysolver_triple_", "", linkage_picker("channel_1101", linkage_df = linkage_metadata_df) ) ) %>% 
        purrr::pmap_dbl(.f = polysolver_triple )
      
    )
  
  
  #Run the estimation summary
  combined_df_solve_estimate <- combined_df_solve_triple %>%
    rowwise() %>%
    
    mutate(
      
      #Calculation of probability sum using only linked probabilities of interest from model 
      prob_across_interest = sum(
        c_across(
          all_of( 
            
            linkage_metadata_df %>% filter( droplet_type == "singlet" | ( (droplet_type != "singlet") & (experimental_linkage == "linked") ) ) %>% pull("fragment_combination") %>%
              gsub("channel_","prob_",.)
            
          ) 
        ),
        na.rm = FALSE), #carry errors through to the sum
      
      # unlinked, semi-linked calculations are droplet estimate values 
      est_across_interest = 1 - prob_across_interest, # use allowed linked probabilities in linkage_metadata_df
      
      theoretical_count_channel_1111 = accepted_droplets * est_across_interest * lambda_theoretical, # use estimate to get number of droplets containing intact plasmid
      
      theoretical_lambda_conc_ul_cluster_channel_1111 = -log( (accepted_droplets - theoretical_count_channel_1111) / accepted_droplets ), # calculate intact plasmid per droplet
      
      theoretical_conc_copies_ul_cluster_channel_1111 = theoretical_lambda_conc_ul_cluster_channel_1111 / !!sym(drop_vol_col) # calculate intact plasmid concentration
    ) %>%
    ungroup()
  
  return(combined_df_solve_estimate)
  
}

#' fwrite with path
#'
#' @param x dataframe
#' @param file path to write to
#' @param ... pass to fwrite
#'
#' @return writes out a file
#' @export FALSE
#'
#' @examples
fwritep <- function(x, file, ...){
  
  dir_path <- dirname(file)
  
  if(dir_path != "." && !dir.exists(dir_path)){
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
  }
  
  fwrite(x, file=file, ...)
  
}


#' silence functions
#'
#' @param expr function call
#' @param return_error_obj return error
#' @note did this because unlinked, semilinked do not evaluate with this model, so error calls happen
#'
#' @return returned expr output
#' @export FALSE
#'
#' @examples
run_silently <- function(expr, return_error_obj = FALSE) {
  # Get the unevaluated expression
  expr_to_eval <- substitute(expr)
  # set up env
  eval_env <- parent.frame()
  
  actual_result <- NULL
  error_occurred <- NULL
  
  # Get messages
  utils::capture.output({ 
    utils::capture.output({ 
      
      # message warning error trycatch and handle
      tryCatch({
      
      actual_result <- withCallingHandlers( 
        
        #Eval
        eval(expr_to_eval, envir = eval_env),
        
        #Message call
        message = function(m) {
          #Silence write to log
          if(isTRUE(getOption("show.error.messages"))){}
          # use handler
          invokeRestart("muffleMessage")
          
        },
        
        #Warning call
        warning = function(w) {
          #Silence write to log
          if(isTRUE(getOption("show.error.messages"))){}
          # use handler
          invokeRestart("muffleWarning")
        }
      
      )}, error = function(e) {
        error_occurred <<- e # Assign to error_occurred in run_silently env
      })
      
    }, type = "message", file = nullfile()) # stder
  }, type = "output", file = nullfile())    # stdout
  
  if (!is.null(error_occurred)) {
    if (return_error_obj) {
      return(invisible(error_occurred))
    } else {
      return(invisible(NULL))
    }
  } else {
    return(actual_result)
  }
}

#' Make barplot
#'
#' @param dfl 
#' @param xval 
#' @param yval 
#'
#' @return
#' @export
#'
#' @examples
prob_barplot <- function(dfl, xval, xlabtext, yval, ylabtext, fillval, filllabtext, showlabs=TRUE, labval,
                         facetval, 
                         x_breaks, facet_col = 4){
  
  #palette
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
  
  facet_fml <- as.formula( paste0("~",facetval) )
  
  s1_test <- ggplot(dfl) +
    geom_bar(aes(x = .data[[xval]],
                 y = .data[[yval]],
                 fill = .data[[fillval]]),
             position = "stack",
             stat = "identity", na.rm=T) +
    scale_fill_manual(values = my_pal) +
    scale_x_log10(breaks = x_breaks) +
    labs(y=ylabtext,
         x=xlabtext,
         fill=filllabtext) +
    facet_wrap(facet_fml, ncol=facet_col) +
    #scale_y_continuous(breaks = c(-1.5, -1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.5)) +
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

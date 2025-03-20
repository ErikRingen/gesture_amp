clean_gesture_data <- function(data_file) {
  
  d <- readxl::read_xlsx(data_file)
  
  d <- d %>% 
    mutate(
      amp = ifelse(AMPLITUDE %in% c("NA", "Unknown"), NA, AMPLITUDE),
      distance = ifelse(`Distance (Physical distance),0=close contact; 1=arm-length distance; 2= distance>1m)` == "NA", NA, `Distance (Physical distance),0=close contact; 1=arm-length distance; 2= distance>1m)`),
      sender_id = `ID (signaller of gesture)`,
      
      reciever_id = ifelse(sapply(strsplit(`Dyad name`,"@"), `[`, 1) == `ID (signaller of gesture)`, sapply(strsplit(`Dyad name`,"@"), `[`, 2), sapply(strsplit(`Dyad name`,"@"), `[`, 1)),
      
      dyad_id = `Dyad name`,
      interaction_id = as.character(`Interaction ID`),
      gesture_id = as.character(`GESTURE TYPE`),
      rank_diff = `DSIst (Rank difference)`,
      bond_strength = `Elost (bond strength)`,
      modality = MODALITY,
      species = SPECIES
    ) %>% 
    mutate(amp = factor(amp, levels = c("Minimum", "Medium", "Maximum"), labels = c("Low", "Medium", "High"), ordered = T),
           distance = factor(distance, levels = c("0", "1", "2"), labels = c("close contact", "arm-length", "> 1m"), ordered = T),
           modality = fct_recode(MODALITY,
                                 Visible = "Visual"
           ))
  
  ## Get individual-level average for within/between comparisons ####
  d_SID <- d %>% 
    group_by(sender_id) %>% 
    summarise(mean_rank_diff = mean(rank_diff, na.rm = T),
              mean_bond_strength = mean(bond_strength, na.rm = T))
  
  d <- d %>% 
    left_join(d_SID) %>% 
    mutate(rank_diff_z = rank_diff - mean_rank_diff,
           bond_strength_z = bond_strength - mean_bond_strength) %>% 
    filter(!(is.na(rank_diff)) & !(is.na(amp)) & !(is.na(modality)))
  
  return(d[,14:ncol(d)])
}

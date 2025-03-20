require(tidyverse)
require(brms)
require(patchwork)

plot_bws_linpred <- function(m_rank_bws, m_bond_bws, d_gesture){
  
### First, rank difference #######
### Within subjects rank diff ##
  d_pid <- d_gesture %>%
    group_by(sender_id) %>% 
    summarise(
      mean_rank_diff = unique(mean_rank_diff),
      low_bond_z = 0 - sd(rank_diff_z, na.rm = T)*2,
      high_bond_z = low_bond_z * -1,
      species = unique(species)
    )
  
  newdata <- expand.grid(
    rank_diff_z = c(0, seq(from = mean(d_gesture$rank_diff_z, na.rm = T) - sd(d_gesture$rank_diff_z, na.rm = T)*2, to = mean(d_gesture$rank_diff_z, na.rm = T) + sd(d_gesture$rank_diff_z, na.rm = T)*2, length.out = 50 )),
    sender_id = d_pid$sender_id,
    modality = unique(d_gesture$modality)
  ) %>%
    left_join(d_pid) %>% 
    # only make predictions within the range observed for the individual
    filter(rank_diff_z <= abs(low_bond_z)) %>% 
    mutate(rank_diff_abs = rank_diff_z + mean_rank_diff,
           pred_id = as.character(1:n()))
  
  preds <- posterior_linpred(m_rank_bws,
                             newdata = newdata,
                             re.form = ~ (1 + mean_rank_diff + rank_diff_z|modality*species) + (1 + rank_diff_z|sender_id)
  )
  
  preds_long <- as.data.frame(preds) %>% 
    set_names(1:ncol(preds)) %>% 
    mutate(samp = 1:n()) %>% 
    pivot_longer(-samp, names_to = "pred_id", values_to = "est")
  
  pred_summary <- preds_long %>% 
    group_by(pred_id) %>% 
    summarise(mean_pred = mean(est)) %>% 
    left_join(newdata) %>% 
    mutate(class = "within-subjects")    
  
### Between subjects rank diff effects ##
d_species <- d_gesture %>% 
  group_by(species) %>% 
  summarise(min_rank_diff = min(mean_rank_diff),
            max_rank_diff = max(mean_rank_diff))

newdata_between <- expand.grid(
  mean_rank_diff = seq(from = min(d_gesture$mean_rank_diff), to = max(d_gesture$mean_rank_diff), length.out = 50),
  rank_diff_z = 0,
  modality = unique(d_gesture$modality),
  species = unique(d_gesture$species)
) %>% 
  left_join(d_species) %>% 
  filter(mean_rank_diff >= min_rank_diff & mean_rank_diff <= max_rank_diff) %>% 
  mutate(pred_id = as.character(1:n()))


preds_between <- posterior_linpred(m_rank_bws,
                                   newdata = newdata_between,
                                   re.form = ~ (1 + mean_rank_diff + rank_diff_z|modality*species)
)


preds_between_long <- as.data.frame(preds_between) %>% 
  set_names(1:ncol(preds_between)) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "pred_id", values_to = "est") 

pred_between_summary <- preds_between_long %>% 
  group_by(pred_id) %>% 
  summarise(mean_pred = mean(est)) %>% 
  left_join(newdata_between) %>% 
  mutate(class = "between-subjects")

pred_summary_mean <- pred_summary %>% 
  filter(rank_diff_z == 0)

## Plot between and within together ####
p_rank <- ggplot(pred_summary, aes(x = rank_diff_abs, y = mean_pred, group = sender_id, color = class)) +
  facet_grid(species ~ modality) + 
  geom_line(alpha = 0.5) +
  geom_point(data = pred_summary_mean, aes(x = rank_diff_abs, y = mean_pred, group = sender_id), fill = "#DC7075", color = "black", pch = 21, size = 1, alpha = 0.5) +
  geom_line(data = pred_between_summary, aes(x = mean_rank_diff, y = mean_pred, group = NA, color = class), lwd = 1, alpha = 0.6) +
  scale_color_manual(values = c("black", "#DC7075")) +
  xlab("Rank Difference (sender - reciever)") +
  ylab("Gesture Amplitude (latent scale)") +
  theme_bw(base_size = 15) +
  theme(
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    legend.title = element_blank(),
    legend.position = "top"
  )
######################################
### Same for bond strength ##########
d_pid <- d_gesture %>%
  group_by(sender_id) %>% 
  summarise(
    mean_bond_strength = unique(mean_bond_strength),
    low_bond_z = 0 - sd(bond_strength_z, na.rm = T)*2,
    high_bond_z = low_bond_z * -1,
    species = unique(species)
  )

newdata <- expand.grid(
  bond_strength_z = c(0, seq(from = mean(d_gesture$bond_strength_z, na.rm = T) - sd(d_gesture$bond_strength_z, na.rm = T)*2, to = mean(d_gesture$bond_strength_z, na.rm = T) + sd(d_gesture$bond_strength_z, na.rm = T)*2, length.out = 50 )),
  sender_id = d_pid$sender_id,
  modality = unique(d_gesture$modality)
) %>%
  left_join(d_pid) %>% 
  # only make predictions within the range observed for the individual
  filter(bond_strength_z <= abs(low_bond_z)) %>% 
  mutate(bond_strength_abs = bond_strength_z + mean_bond_strength,
         pred_id = as.character(1:n()))

preds <- posterior_linpred(m_bond_bws,
                           newdata = newdata,
                           re.form = ~ (1 + mean_bond_strength + bond_strength_z|modality*species) + (1 + bond_strength_z|sender_id)
)

preds_long <- as.data.frame(preds) %>% 
  set_names(1:ncol(preds)) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "pred_id", values_to = "est")

pred_summary <- preds_long %>% 
  group_by(pred_id) %>% 
  summarise(mean_pred = mean(est)) %>% 
  left_join(newdata) %>% 
  mutate(class = "within-subjects")

### Mean rank diff effects #####
d_species <- d_gesture %>% 
  group_by(species) %>% 
  summarise(min_bond_strength = min(mean_bond_strength),
            max_bond_strength = max(mean_bond_strength))

newdata_between <- expand.grid(
  mean_bond_strength = seq(from = min(d_gesture$mean_bond_strength), to = max(d_gesture$mean_bond_strength), length.out = 50),
  bond_strength_z = 0,
  modality = unique(d_gesture$modality),
  species = unique(d_gesture$species)
) %>% 
  left_join(d_species) %>% 
  filter(mean_bond_strength >= min_bond_strength & mean_bond_strength <= max_bond_strength) %>% 
  mutate(pred_id = as.character(1:n()))


preds_between <- posterior_linpred(m_bond_bws,
                                   newdata = newdata_between,
                                   re.form = ~ (1 + mean_bond_strength + bond_strength_z|modality*species)
)


preds_between_long <- as.data.frame(preds_between) %>% 
  set_names(1:ncol(preds_between)) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "pred_id", values_to = "est") 

pred_between_summary <- preds_between_long %>% 
  group_by(pred_id) %>% 
  summarise(mean_pred = mean(est)) %>% 
  left_join(newdata_between) %>% 
  mutate(class = "between-subjects")

pred_summary_mean <- pred_summary %>% 
  filter(bond_strength_z == 0)

p_bond <- ggplot(pred_summary, aes(x = bond_strength_abs, y = mean_pred, group = sender_id, color = class)) +
  facet_grid(species ~ modality) + 
  geom_line(alpha = 0.5) +
  geom_point(data = pred_summary_mean, aes(x = bond_strength_abs, y = mean_pred, group = sender_id), fill = "#549C9E", color = "black", pch = 21, size = 1, alpha = 0.5) +
  geom_line(data = pred_between_summary, aes(x = mean_bond_strength, y = mean_pred, group = NA, color = class), lwd = 1, alpha = 0.6) +
  scale_color_manual(values = c("black", "#549C9E")) +
  xlab("Bond Strength") +
  ylab("Gesture Amplitude (latent scale)") +
  theme_bw(base_size = 15) +
  theme(
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    legend.title = element_blank(),
    legend.position = "top"
  )

p_both <- p_rank / p_bond

ggsave("figures/plot_bws_linpred.pdf", plot = p_both, width = 8.5, height = 11)
return(p_rank / p_bond)
}

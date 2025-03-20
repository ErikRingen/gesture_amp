require(tidyverse)
require(brms)
require(patchwork)
require(rethinking)

plot_bs_smd <- function(m_rank_bws, m_bond_bws, d_gesture){

### Between subjects #####################################  
# Define sequence of predictors
mean_bond <-  mean(d_gesture$mean_bond_strength, na.rm = T)
sd_bond <- sd(d_gesture$mean_bond_strength, na.rm = T)

mean_rank <-  mean(d_gesture$mean_rank_diff, na.rm = T)
sd_rank <- sd(d_gesture$mean_rank_diff, na.rm = T)
  
bond_seq <- c(mean_bond, mean_bond + sd_bond)
rank_seq <- c(mean_rank, mean_rank + sd_rank)
  
#### Predict for rank difference #####
newdata <- expand.grid(
  mean_rank_diff = rank_seq,
  rank_diff_z = 0,
  modality = unique(d_gesture$modality[!(is.na(d_gesture$modality))]),
  species = unique(d_gesture$species)
)

newdata$pred_row <- 1:nrow(newdata)

pred <- posterior_linpred(
  m_rank_bws,
  newdata = newdata,
  re_formula = ~ (1 + mean_rank_diff|modality*species)
)

n_samps <- brms::ndraws(m_rank_bws)

pred <- array(pred, dim = c(n_samps, nrow(newdata)), dimnames = list( samps = 1:n_samps, pred_row = 1:nrow(newdata) ))

# Pivot array long
pred_rank_B_long <- pred %>% 
  cubelyr::as.tbl_cube(met_name = "est") %>% 
  as_tibble() %>% 
  left_join(newdata) %>% 
  select(-rank_diff_z)

pred_rank_B_diff <- pred_rank_B_long %>% 
  group_by(samps, modality, species) %>% 
  mutate(diff = est[mean_rank_diff == (mean_rank + sd_rank)] - est[mean_rank_diff == mean_rank],
         smd = diff * (sqrt(3)/pi) ) 

post_probs <- pred_rank_B_diff %>% 
  group_by(modality, species) %>% 
  summarise(PP = PP(diff), x_pos = ifelse(mean(diff) > 0, 0.9, -0.9)) %>% 
  mutate(PP_label = paste0("PP = ", PP))

p_rank_bs <- ggplot(pred_rank_B_diff, aes(x = diff, y = modality)) +
  facet_grid(species ~ .) +
  geom_hline(yintercept = c(1:4), alpha = 0.6) +
  geom_density_ridges2(scale = 0.8, rel_min_height = 0.01, fill = "#DC7075", color = NA) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_text(data = post_probs, aes(x = x_pos, y = modality, label = PP_label), vjust = -3) +
  annotate("blank", x = 0, y = 4) +
  xlab("Standardized Mean Diff. (+ 1SD Rank Diff)") +
  ylab("") +
  labs(fill = "Amplitude") +
  labs(color = "Amplitude") +
  ggtitle("Between-subjects") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.border = element_rect(colour = "black"))

#### Predict for bond strength #####
newdata <- expand.grid(
  mean_bond_strength = bond_seq,
  bond_strength_z = 0,
  modality = unique(d_gesture$modality[!(is.na(d_gesture$modality))]),
  species = unique(d_gesture$species)
)

newdata$pred_row <- 1:nrow(newdata)

pred <- posterior_linpred(
  m_bond_bws,
  newdata = newdata,
  re_formula = ~ (1 + mean_bond_strength|modality*species)
)

n_samps <- brms::ndraws(m_bond_bws)

pred <- array(pred, dim = c(n_samps, nrow(newdata)), dimnames = list( samps = 1:n_samps, pred_row = 1:nrow(newdata) ))

# Pivot array long
pred_bond_B_long <- pred %>% 
  cubelyr::as.tbl_cube(met_name = "est") %>% 
  as_tibble() %>% 
  left_join(newdata)

pred_bond_B_diff <- pred_bond_B_long %>% 
  group_by(samps, modality, species) %>% 
  mutate(diff = est[mean_bond_strength == (mean_bond + sd_bond)] - est[mean_bond_strength == mean_bond],
         smd = diff * (sqrt(3)/pi) ) 

post_probs <- pred_bond_B_diff %>% 
  group_by(modality, species) %>% 
  summarise(PP = PP(diff), x_pos = ifelse(mean(diff) > 0, 0.9, -0.9)) %>% 
  mutate(PP_label = paste0("PP = ", PP))

#post_probs$PP_label <- ifelse(round(post_probs$PP, 3) == 1, "PP ≈ 1", post_probs$PP_label)

p_bond_bs <- ggplot(pred_bond_B_diff, aes(x = diff, y = modality)) +
  facet_grid(species ~ .) +
  geom_hline(yintercept = c(1:4), alpha = 0.6) +
  geom_density_ridges2(scale = 0.8, rel_min_height = 0.01, fill = "#549C9E", color = NA) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_text(data = post_probs, aes(x = x_pos, y = modality, label = PP_label), vjust = -3) +
  annotate("blank", x = 0, y = 4) +
  xlab("Standardized Mean Diff. (+ 1SD Bond Strength)") +
  ylab("") +
  labs(fill = "Amplitude") +
  labs(color = "Amplitude") +
  ggtitle("Between-subjects") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.border = element_rect(colour = "black"))

#########################################
##### Within subjects ###################
# Define sequence of predictors
mean_bond <-  0
sd_bond <- sd(d_gesture$bond_strength_z, na.rm = T)

mean_rank <-  0
sd_rank <- sd(d_gesture$rank_diff_z, na.rm = T)

bond_seq <- c(mean_bond, mean_bond + sd_bond)
rank_seq <- c(mean_rank, mean_rank + sd_rank)

#### Predict for rank difference #####
newdata <- expand.grid(
  rank_diff_z = rank_seq,
  mean_rank_diff = 0,
  modality = unique(d_gesture$modality[!(is.na(d_gesture$modality))]),
  species = unique(d_gesture$species)
)

newdata$pred_row <- 1:nrow(newdata)

pred <- posterior_linpred(
  m_rank_bws,
  newdata = newdata,
  re_formula = ~ (1 + rank_diff_z|modality*species)
)

n_samps <- brms::ndraws(m_rank_bws)

pred <- array(pred, dim = c(n_samps, nrow(newdata)), dimnames = list( samps = 1:n_samps, pred_row = 1:nrow(newdata) ))

# Pivot array long
pred_rank_B_long <- pred %>% 
  cubelyr::as.tbl_cube(met_name = "est") %>% 
  as_tibble() %>% 
  left_join(newdata)

pred_rank_B_diff <- pred_rank_B_long %>% 
  group_by(samps, modality, species) %>% 
  mutate(diff = est[rank_diff_z == (mean_rank + sd_rank)] - est[rank_diff_z == mean_rank],
         smd = diff * (sqrt(3)/pi) ) 

post_probs <- pred_rank_B_diff %>% 
  group_by(modality, species) %>% 
  summarise(PP = PP(diff), x_pos = ifelse(mean(diff) > 0, 0.55, -0.9)) %>% 
  mutate(PP_label = paste0("PP = ", PP))

p_rank_ws <- ggplot(pred_rank_B_diff, aes(x = diff, y = modality)) +
  facet_grid(species ~ .) +
  geom_hline(yintercept = c(1:4), alpha = 0.6) +
  geom_density_ridges2(scale = 0.8, rel_min_height = 0.01, fill = "#DC7075", color = NA) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_text(data = post_probs, aes(x = x_pos, y = modality, label = PP_label), vjust = -3) +
  annotate("blank", x = 0, y = 4) +
  xlab("Standardized Mean Diff. (+ 1SD Rank Diff)") +
  ylab("") +
  labs(fill = "Amplitude") +
  labs(color = "Amplitude") +
  ggtitle("Within-subjects") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.border = element_rect(colour = "black"))

#### Predict for bond strength #####
newdata <- expand.grid(
  bond_strength_z = bond_seq,
  mean_bond_strength = 0,
  modality = unique(d_gesture$modality[!(is.na(d_gesture$modality))]),
  species = unique(d_gesture$species)
)

newdata$pred_row <- 1:nrow(newdata)

pred <- posterior_linpred(
  m_bond_bws,
  newdata = newdata,
  re_formula = ~ (1 + bond_strength_z|modality*species)
)

n_samps <- brms::ndraws(m_bond_bws)

pred <- array(pred, dim = c(n_samps, nrow(newdata)), dimnames = list( samps = 1:n_samps, pred_row = 1:nrow(newdata) ))

# Pivot array long
pred_bond_B_long <- pred %>% 
  cubelyr::as.tbl_cube(met_name = "est") %>% 
  as_tibble() %>% 
  left_join(newdata)

pred_bond_B_diff <- pred_bond_B_long %>% 
  group_by(samps, modality, species) %>% 
  mutate(diff = est[bond_strength_z == (mean_bond + sd_bond)] - est[bond_strength_z == mean_bond],
         smd = diff * (sqrt(3)/pi) ) 

post_probs <- pred_bond_B_diff %>% 
  group_by(modality, species) %>% 
  summarise(PP = PP(diff), x_pos = ifelse(mean(diff) > 0, 0.55, -0.9)) %>% 
  mutate(PP_label = paste0("PP = ", PP))

post_probs$PP_label <- ifelse(round(post_probs$PP, 3) == 1, "PP ≈ 1", post_probs$PP_label)

p_bond_ws <- ggplot(pred_bond_B_diff, aes(x = diff, y = modality)) +
  facet_grid(species ~ .) +
  geom_hline(yintercept = c(1:4), alpha = 0.6) +
  geom_density_ridges2(scale = 0.8, rel_min_height = 0.01, fill = "#549C9E", color = NA) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_discrete(expand = c(0.00, 0)) + 
  geom_text(data = post_probs, aes(x = x_pos, y = modality, label = PP_label), vjust = -3) +
  annotate("blank", x = 0, y = 4) +
  xlab("Standardized Mean Diff. (+ 1SD Bond Strength)") +
  ylab("") +
  labs(fill = "Amplitude") +
  labs(color = "Amplitude") +
  ggtitle("Within-subjects") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.border = element_rect(colour = "black"))


#########################
p_all <- (p_rank_ws + p_bond_ws) / (p_rank_bs + p_bond_bs) 

ggsave("figures/plot_bws_smd.pdf", plot = p_all, width = 11, height = 8.5)
return(p_all)
}

require(tidyverse)
require(brms)
require(patchwork)
require(rethinking)

plot_ws_ordpred <- function(m_rank_bws, m_bond_bws, d_gesture){

# Define sequence of predictors
mean_bond_z <-  mean(d_gesture$bond_strength_z, na.rm = T)
sd_bond_z <- sd(d_gesture$bond_strength_z, na.rm = T)

mean_rank_z <-  mean(d_gesture$rank_diff_z, na.rm = T)
sd_rank_z <- sd(d_gesture$rank_diff_z, na.rm = T)
  
bond_seq <- seq(from = mean_bond_z - 2*sd_bond_z, to = mean_bond_z + sd_bond_z*2, length.out = 30) - mean_bond_z
rank_seq <- seq(from = mean_rank_z - 2*sd_rank_z, to = mean_rank_z + sd_rank_z*2, length.out = 30)
  
#### Predict for rank difference #####
newdata <- expand.grid(
  rank_diff_z = rank_seq,
  mean_rank_diff = 0,
  modality = unique(d_gesture$modality[!(is.na(d_gesture$modality))]),
  species = unique(d_gesture$species)
)

newdata$pred_row <- 1:nrow(newdata)

pred <- posterior_epred(
  m_rank_bws,
  newdata = newdata,
  re_formula = ~ (1 + rank_diff_z|modality*species)
)

n_samps <- brms::ndraws(m_rank_bws)

pred <- array(pred, dim = c(n_samps, nrow(newdata), 3), dimnames = list( samps = 1:n_samps, pred_row = 1:nrow(newdata), amp = c("Low", "Medium", "High") ))

# Pivot array long
pred_rank_B_long <- pred %>% 
  cubelyr::as.tbl_cube(met_name = "est") %>% 
  as_tibble() %>% 
  left_join(newdata) 

pred_rank_B_summary <- pred_rank_B_long %>% 
  group_by(amp, rank_diff_z, modality, species) %>% 
  summarise(med = median(est),
            lower = PI(est, prob = 0.95)[1],
            upper = PI(est, prob = 0.95)[2])

pred_rank_B_summary$amp <- factor(pred_rank_B_summary$amp, levels = c("Low", "Medium", "High"))

pred_rank_B_summary$rank_z_score <- (pred_rank_B_summary$rank_diff_z - mean_rank_z) / sd_rank_z

p_rank <- ggplot(pred_rank_B_summary, aes(x = rank_z_score, y = med, color = amp, fill = amp)) +
  facet_grid(species ~ modality) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, color = NA) + 
  geom_line() + 
  scale_color_manual(values = (c("#F4B68E", "#DC7075", "#6F274A"))) +
  scale_fill_manual(values = (c("#F4B68E", "#DC7075", "#6F274A"))) +
  xlab("Rank Diff. (within-subjects z-score)") +
  ylab("Pr(Amplitude)") +
  labs(fill = "Amplitude") +
  labs(color = "Amplitude") +
  ggtitle("") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.border = element_rect(colour = "black"))

### Predict for bond strength ######
newdata <- expand.grid(
  bond_strength_z = bond_seq,
  mean_bond_strength = 0,
  modality = unique(d_gesture$modality[!(is.na(d_gesture$modality))]),
  species = unique(d_gesture$species)
)

newdata$pred_row <- 1:nrow(newdata)

pred <- posterior_epred(
  m_bond_bws,
  newdata = newdata,
  re_formula = ~ (1 + bond_strength_z|modality*species)
)

n_samps <- brms::ndraws(m_bond_bws)

pred <- array(pred, dim = c(n_samps, nrow(newdata), 3), dimnames = list( samps = 1:n_samps, pred_row = 1:nrow(newdata), amp = c("Low", "Medium", "High") ))

# Pivot array long
pred_bond_B_long <- pred %>% 
  cubelyr::as.tbl_cube(met_name = "est") %>% 
  as_tibble() %>% 
  left_join(newdata) 

pred_bond_B_summary <- pred_bond_B_long %>% 
  group_by(amp, bond_strength_z, modality, species) %>% 
  summarise(med = median(est),
            lower = PI(est, prob = 0.95)[1],
            upper = PI(est, prob = 0.95)[2])

pred_bond_B_summary$amp <- factor(pred_bond_B_summary$amp, levels = c("Low", "Medium", "High"))

pred_bond_B_summary$bond_z_score <- (pred_bond_B_summary$bond_strength_z - mean_bond_z) / sd_bond_z

p_bond <- ggplot(pred_bond_B_summary, aes(x = bond_z_score, y = med, color = amp, fill = amp)) +
  facet_grid(species ~ modality) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, color = NA) + 
  geom_line() + 
  scale_color_manual(values = (c("#A5DAC1", "#549C9E", "#123F59"))) +
  scale_fill_manual(values = (c("#A5DAC1", "#549C9E", "#123F59"))) +
  xlab("Bond Strength (within-subjects z-score)") +
  ylab("Pr(Amplitude)") +
  labs(fill = "Amplitude") +
  labs(color = "Amplitude") +
  ggtitle("") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"),
        panel.border = element_rect(colour = "black"))


p_both <- p_rank / p_bond

ggsave("figures/plot_ws_ordpred.pdf", plot = p_both, width = 8.5, height = 11)
return(p_rank / p_bond)
}

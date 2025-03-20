library(targets)
source("R/functions.R")

# Dependencies for this pipeline
tar_option_set(packages = c(
  "tidyverse",
  "brms",
  "patchwork",
  "rethinking"
))

list(
  # Extract raw data
  tar_target(gesture_data_file,
             "data_raw/GENTY-APES_GESTURE_AMPLITUDE_DATA.xlsx", format = "file"),
  
  # Clean data
  tar_target(d_gesture,
             clean_gesture_data(gesture_data_file)),
  
  # Bond strength model, between and within individual
  tar_target(m_bond_bws,
             fit_bond_bws(d_gesture)),
  
  # Rank diff model, between and within individual
  tar_target(m_rank_bws,
             fit_rank_bws(d_gesture)),
  
  # Make plot within and between subjects effects of rank and bond strength, latent scale
  tar_target(p_bws_linpred,
             plot_bws_linpred(m_rank_bws, m_bond_bws, d_gesture)),
  
  # Make predictive plots on original scale
  tar_target(p_ws_ordpred,
             plot_ws_ordpred(m_rank_bws, m_bond_bws, d_gesture)),
  
  tar_target(p_bs_ordpred,
             plot_bs_ordpred(m_rank_bws, m_bond_bws, d_gesture)),
  
  # Standardized mean diff plots
  tar_target(p_bws_smd,
             plot_bws_smd(m_rank_bws, m_bond_bws, d_gesture))
  
)

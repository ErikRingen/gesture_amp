# First, source all files in dir to get bigger named functions
functions <- list.files(
  path = "R",
  pattern = "*.R",
  recursive = TRUE)

sapply(paste0("R/", functions[functions != "functions.R"]), source, .GlobalEnv)
rm(functions)

# Calculate posterior probability, generic
PP <- function(x) {
  if (mean(x) > 0) pp <- mean(x > 0)
  else pp <- mean(x < 0)
  
  return(round(pp, digits = 3))
}

# Generic regularizing priors
priors <- c(
  brms::prior(normal(0,1), class = b),
  brms::prior(exponential(1), class = sd)
)

fit_rank_bws <- function(data, chains = 8, cores = 8, iter = 1000, incl_dist = F) {
  
  m <- brm(
    amp ~ mean_rank_diff + rank_diff_z + (1 + mean_rank_diff + rank_diff_z|modality*species) + (1|interaction_id) + (1 + rank_diff_z|sender_id) + (1 + rank_diff_z|reciever_id) + (1|dyad_id),
    family = "cumulative",
    prior = priors,
    data = data,
    chains = chains,
    cores = cores,
    iter = iter,
    adapt_delta = 0.98,
    backend = "cmdstanr"
  )
  
  return(m)
}


fit_bond_bws <- function(data, chains = 8, cores = 8, iter = 1000, incl_dist = F) {
  
  m <- brm(
    amp ~ mean_bond_strength + bond_strength_z + (1 + mean_bond_strength + bond_strength_z|modality*species) + (1|interaction_id) + (1 + bond_strength_z|sender_id) + (1 + bond_strength_z|reciever_id) + (1|dyad_id),
    family = "cumulative",
    prior = priors,
    data = data,
    chains = chains,
    cores = cores,
    iter = iter,
    adapt_delta = 0.98,
    backend = "cmdstanr"
  )
  
  return(m)
}


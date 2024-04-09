library(tidyverse)
library(cmdstanr)

species_list <- c(`Barking Marsh Frog` = 13059,
                  `Eastern Sign-bearing Froglet` = 13131,
                  `Common Froglet` = 13134,
                  `Peron's Tree Frog` = 13204)

nspecies <- length(species_list)

stan_data <- readRDS("data/stan_data.rds")
cont_model_ms <- cmdstan_model("stan/continuous_model_ms.stan")

ni <- 500 #samples
nw <- 500 #warmups
nc <- 8 #chains

fit <- cont_model_ms$sample(data = stan_data,
                            chains = nc,
                            parallel_chains = nc,
                            show_messages = TRUE,
                            save_warmup = FALSE,
                            iter_sampling = ni,
                            iter_warmup = nw, adapt_delta = 0.95)

# Occupancy of sites not seen vs occupancy of sites seen
psi_summary <- fit$summary("psi")
psi_summary_sp <- list()
for(i in 1:nspecies) {
  var_names <- paste0("psi[",1:stan_data$nsites, ",", i, "]")
  psi_summary_sp[[i]] <- psi_summary %>%
    filter(variable %in% var_names) %>%
    mutate(species = names(species_list)[i],
           any_seen = stan_data$any_seen[,i])
}
psi_summary <- bind_rows(psi_summary_sp)

psi_summary %>%
  group_by(any_seen, species) %>%
  summarise(psi = mean(mean))

# Ovule-gated passenger RI empirical modelling template
# This is a starting scaffold, not a final analysis script.
# It assumes you have lineages, crosses, and a phylogeny in comparable formats.

library(tidyverse)
# library(brms)
# library(ape)
# library(phylolm)

lineages <- read_csv("data_templates/lineage_traits_template.csv")
crosses  <- read_csv("data_templates/crosses_template.csv")

# 1. Build K or a latent/proxy K.
# If pollen load is available, use it. If not, use pollen production times a
# sensitivity value psi. Run this for a grid of psi values.
psi <- 0.02
lineages <- lineages %>%
  mutate(
    pollen_eff_proxy = case_when(
      !is.na(pollen_load_mean) ~ pollen_load_mean,
      !is.na(pollen_production_per_flower) ~ psi * pollen_production_per_flower,
      TRUE ~ NA_real_
    ),
    K_proxy = pollen_eff_proxy / ovules_per_ovary,
    logK_proxy = log(K_proxy),
    logO = log(ovules_per_ovary)
  )

# 2. Join maternal and paternal traits onto each cross.
cross_dat <- crosses %>%
  left_join(lineages %>% select(species, maternal_ovules = ovules_per_ovary,
                               maternal_logO = logO,
                               maternal_logK = logK_proxy),
            by = c("maternal_species" = "species")) %>%
  left_join(lineages %>% select(species, paternal_ovules = ovules_per_ovary,
                               paternal_logO = logO,
                               paternal_logK = logK_proxy),
            by = c("paternal_species" = "species")) %>%
  mutate(
    seed_failures = ovules_exposed - seeds_matured,
    preRI_observed = 1 - seeds_matured / ovules_exposed
  )

# 3. Reciprocal-pair asymmetry table.
# This assumes both directions are present for a pair. It creates a pair_id
# independent of cross direction.
asym_dat <- cross_dat %>%
  mutate(
    sp_min = pmin(maternal_species, paternal_species),
    sp_max = pmax(maternal_species, paternal_species),
    pair_id = paste(sp_min, sp_max, sep = "__")
  ) %>%
  group_by(pair_id) %>%
  filter(n() == 2) %>%
  summarise(
    species_i = first(maternal_species),
    species_j = first(paternal_species),
    asym = first(preRI_observed) - last(preRI_observed),
    delta_logK = first(maternal_logK) - last(maternal_logK),
    delta_logO_proxy = last(maternal_logO) - first(maternal_logO),
    genetic_distance = mean(genetic_distance, na.rm = TRUE),
    .groups = "drop"
  )

# Basic non-phylogenetic diagnostic. This is only a first look.
lm_asym_K <- lm(asym ~ delta_logK + genetic_distance, data = asym_dat)
summary(lm_asym_K)

lm_asym_O <- lm(asym ~ delta_logO_proxy + genetic_distance, data = asym_dat)
summary(lm_asym_O)

# 4. Denominator-aware model for raw seed counts.
# With brms, a beta-binomial model can retain ovule denominators. Phylogenetic
# covariance can be added through gr(..., cov = A) once the tree and covariance
# matrix are prepared.

# Example skeleton only:
# tree <- read.tree("data_templates/phylogeny_template.nwk")
# A <- ape::vcv.phylo(tree, corr = TRUE)
#
# fit_pre <- brm(
#   seeds_matured | trials(ovules_exposed) ~
#     maternal_logK + genetic_distance + maternal_logK:genetic_distance +
#     (1 | gr(maternal_species, cov = A)) +
#     (1 | gr(paternal_species, cov = A)) +
#     (1 | maternal_species:paternal_species),
#   data = cross_dat,
#   family = beta_binomial2(),
#   data2 = list(A = A),
#   chains = 4, cores = 4
# )
# summary(fit_pre)

# 5. Sensitivity loop over psi.
psi_grid <- c(0.001, 0.003, 0.01, 0.02, 0.05, 0.10)
# Repeat the construction of K_proxy and the asymmetry model for each psi.

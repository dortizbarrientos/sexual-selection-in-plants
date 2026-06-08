# ==============================================================================
# FINAL PRODUCTION SCRIPT: PHYLOGENETIC MIXED MODEL (PMM) FOR DYADIC TRAITS
# Framework: Bayesian Skew-Normal Nested Multi-Membership Model via brms
# Target: Modeling Reproductive Isolation (RI) with bounded and hybrid vigour traits
# Daniel Ortiz-Barrientos
# June 8, 2026
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: REQUIRED LIBRARIES & CORE DEPENDENCIES
# ------------------------------------------------------------------------------
library(ape)            # Core phylogenetics tool for tree parsing and branch tracking
library(Matrix)         # Advanced matrix transformations for kinship mapping
library(brms)           # High-level interface to Stan for Bayesian Generalized Linear Mixed Models
library(tidyverse)      # Standard data management and manipulation pipeline (dplyr, ggplot2, etc.)
library(modelsummary)   # Package to construct automated publication tables directly from model objects
library(gt)             # Great Tables formatting engine for publishing quality outputs
library(marginaleffects)# Advanced posterior marginalization tools ignoring random effects

# ------------------------------------------------------------------------------
# SECTION 2: DATA PREPARATION & SISTER-CLADE TRACKING
# ------------------------------------------------------------------------------
# PRE-REQUISITE ASSUMPTION: 
# You have an empirical dataset 'my_dyad_data' containing columns: sp1, sp2, 
# genetic_distance, ovule_number_diff, ovule_category_mismatch, and reproductive_isolation.
# You also have your matching phylogenetic tree loaded as 'my_tree'.

# Extract the Most Recent Common Ancestor (MRCA) node index for every single species pair.
# This isolates the variation shared exclusively by deeply nested sister-clades.
my_dyad_data$mrca_node <- sapply(1:nrow(my_dyad_data), function(i) {
  getMRCA(my_tree, tip = c(my_dyad_data$sp1[i], my_dyad_data$sp2[i]))
})

# Format structural columns into proper factors for the Bayesian engine
my_dyad_data <- my_dyad_data %>%
  mutate(
    mrca_node = as.factor(paste0("node_", mrca_node)), # Categorical factor for node grouping
    sp1 = as.factor(sp1),                            # Factor for parent species 1
    sp2 = as.factor(sp2),                            # Factor for parent species 2
    ovule_category_mismatch = as.factor(ovule_category_mismatch) # Categorical factor (0 = same, 1 = mismatch)
  )

# Calculate the inverse phylogenetic covariance relationship matrix from your tree tips.
# This scales the branch lengths and maps tip distances back to standard correlation blocks.
inv_phylo_matrix <- inverseA(my_tree, scale = TRUE)$Ainv

# ------------------------------------------------------------------------------
# SECTION 3: MODEL SPECIFICATION & MULTI-CORE ESTIMATION
# ------------------------------------------------------------------------------
# Define the structural model formula. 
# Includes continuous and categorical fixed terms alongside nested random adjustments.
final_formula <- bf(
  reproductive_isolation ~ 
    genetic_distance + 
    ovule_number_diff + 
    ovule_category_mismatch + 
    (1 | mrca_node) +         # Random Effect 1: Baseline isolation offsets shared within deep sister-clades
    (1 | mm(sp1, sp2))        # Random Effect 2: Multi-membership tracking for cross-species parentage linked to phylo matrix
)

# Fit the Bayesian PMM via the Stan backend
final_pmm_model <- brm(
  formula = final_formula,
  data = my_dyad_data,
  data2 = list(A = inv_phylo_matrix), # Connects the inverse phylo correlation matrix directly to the multi-membership term
  family = skew_normal(),            # Rigorously maps bounded traits alongside negative hybrid vigour asymmetric tails
  prior = c(
    prior(normal(0, 2), class = "b"),           # Regularizing weakly informative priors for fixed effects
    prior(student_t(3, 0, 1), class = "sd"),     # Regularizing priors for variance parameters
    prior(normal(0, 4), class = "alpha")         # Weakly informative prior targeting distribution skewness
  ),
  chains = 4,                       # 4 Markov chains for standard publication-grade diagnostic verification
  cores = 4,                        # Leverages multi-core setup. Match this to your computer's high-spec cores
  iter = 3000,                      # Generous iteration parameters to ensure chain convergence over large matrix matrices
  warmup = 1000,                    # Drop initial 1000 noisy samples to isolate target posteriors
  control = list(adapt_delta = 0.95, max_treedepth = 15), # Heightened search settings to eliminate divergent transitions
  backend = "rstan"
)

# Output summary diagnostics, parameter estimates, Rhat diagnostics, and effective sample sizes (ESS)
summary(final_pmm_model)

# ------------------------------------------------------------------------------
# SECTION 4: PUBLICATION-READY SUMMARY TABLE GENERATION
# ------------------------------------------------------------------------------
# Mapping dictionary to format raw variable names into clean table row labels
var_labels <- c(
  "Intercept" = "Baseline Intercept (Constant)",
  "genetic_distance" = "Genetic Distance Gradient",
  "ovule_number_diff" = "Ovule Count Disparity (Continuous)",
  "ovule_category_mismatch1" = "Ovule Category Mismatch (Few vs. Many)",
  "alpha" = "Skewness (Hybrid Vigour Tail Asymmetry)"
)

# Extract and style model metrics into a clean HTML/LaTeX/gt object
publication_table <- modelsummary(
  final_pmm_model,
  metrics = c("estimate", "conf.low", "conf.high"),
  coef_map = var_labels,
  statistic = "[{conf.low}, {conf.high}]", # Neatly formats 95% Bayesian Credible Intervals inside standard brackets
  fmt = 3,                                  # Restricts numeric decimal values to a consistent 3 decimal places
  title = "Table 1: Phylogenetic Mixed Model (PMM) Estimates for Reproductive Isolation",
  notes = "Values inside brackets represent the 95% Bayesian Credible Intervals. Clade and Tip variances are implicitly calculated.",
  output = "gt"                             # Outputs a table object managed through the gt engine
)

# Apply formal academic styling layouts to the final table object
styled_table <- publication_table %>%
  tab_options(
    table.border.top.color = "black",       # Traditional double horizontal rule design boundaries
    table.border.bottom.color = "black",
    table_body.border.bottom.color = "black",
    heading.align = "left"                  # Align main headers cleanly to the left margins
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()       # Bold table column descriptive text labels
  )

# Display table inside the interactive viewer or export
print(styled_table)

# ------------------------------------------------------------------------------
# SECTION 5: ADVANCED POPULATION-SCALE MARGINAL VISUALIZATION
# ------------------------------------------------------------------------------
# Extract pure population-level expectations.
# Calculates the absolute effect of ovule disparities across three explicit genetic distance threshold slices.
marginal_predictions <- predictions(
  final_pmm_model,
  variables = "ovule_number_diff",
  datagrid(genetic_distance = c(0.1, 0.4, 0.8)), # Maps early, mid, and deep divergence evolutionary thresholds
  re_formula = NA                               # Essential parameter: sets all random effects to NA to compute true global population trends
)

# Build a publication-grade graphic tracking the interacting gradients
marginal_effects_plot <- ggplot(
  marginal_predictions, 
  aes(x = ovule_number_diff, y = estimate, color = as.factor(genetic_distance), fill = as.factor(genetic_distance))
) +
  # Generate shaded ribbon bands mapping out the 95% Bayesian Credible Intervals
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  # Overlay the primary continuous expected regression tracking paths
  geom_line(linewidth = 1.2) +
  # Custom aesthetic layout enhancements using color palettes optimized for scannability
  scale_color_brewer(palette = "Set1", name = "Genetic Distance Slice") +
  scale_fill_brewer(palette = "Set1", name = "Genetic Distance Slice") +
  labs(
    title = "Marginal Effect of Ovule Disparity on Reproductive Isolation",
    subtitle = "Calculated across evolutionary distance gradients using an asymmetric Skew-Normal framework.",
    x = "Phenotypic Disparity (Absolute Ovule Number Difference)",
    y = "Predicted Reproductive Isolation (RI)"
  ) +
  theme_classic(base_size = 13) +           # Clean theme framework with enhanced readable base font dimensions
  theme(
    legend.position = "bottom",              # Reposition legend to maximize graphic space width
    plot.title = element_text(face = "bold", size = 15),
    axis.line = element_line(color = "black")# Clean, classic black plot boundaries
  )

# Output final visualization map
print(marginal_effects_plot)

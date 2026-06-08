"""
==============================================================================
CONSOLIDATED PYTHON SIMULATION & ESTIMATION PIPELINE
Framework: Custom Bayesian Skew-Normal PMM via PyMC & Bambi
Target: Sensitivity testing of Reproductive Isolation (RI) data configurations
==============================================================================
"""

import numpy as np
import pandas as pd
import pymc as pm
import bambi as bmb
import arviz as az
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# ------------------------------------------------------------------------------
# SECTION 1: SYSTEMATIC PAIRWISE DATA SIMULATOR
# ------------------------------------------------------------------------------
def simulate_pairwise_dataset(
    n_species=30, 
    n_dyads=150, 
    genus_blur_prop=0.1, 
    few_many_overlap=0.2, 
    hybrid_vigour_prop=0.15
):
    """
    Simulates a comparative species-pair dataset with embedded phylogenetic matrix signals,
    trait measurement imputation, categorical class noise, and an asymmetric hybrid vigour tail.
    """
    np.random.seed(42)
    species_names = [f"species_{i}" for i in range(n_species)]
    
    # 1. Generate a mock phylogenetic covariance matrix using a random distance layout
    raw_coords = np.random.rand(n_species, 2)
    dist_matrix = squareform(np.linalg.norm(raw_coords[:, None] - raw_coords, axis=2))
    max_d = dist_matrix.max()
    phylo_cov = max_d - dist_matrix
    np.fill_diagonal(phylo_cov, phylo_cov.diagonal() + 1e-4) # Guard positive-definiteness
    
    # Generate true continuous trait (ovules) following tree covariance structure
    true_ovules_raw = np.random.multivariate_normal(np.zeros(n_species), phylo_cov)
    true_ovules = np.round(np.exp((true_ovules_raw - true_ovules_raw.mean()) / true_ovules_raw.std() * 0.5 + 3))
    
    # Create simple genus buckets using a hierarchical tree split
    Z = linkage(dist_matrix, 'complete')
    genus_ids = fcluster(Z, t=max(3, n_species // 5), criterion='maxclust')
    
    # Assemble species reference matrix
    species_df = pd.DataFrame({
        'species': species_names,
        'true_ovules': true_ovules,
        'genus': genus_ids
    })
    
    # Inject imputation blur (replacing individual counts with genus averages)
    genus_means = species_df.groupby('genus')['true_ovules'].transform('mean').round()
    mask = np.random.rand(n_species) < genus_blur_prop
    species_df['observed_ovules'] = np.where(mask, genus_means, species_df['true_ovules'])
    
    # Apply categorical overlap noise to split into "Few" vs "Many" bins
    median_val = species_df['observed_ovules'].median()
    noise = np.random.normal(0, few_many_overlap * median_val, n_species)
    species_df['ovule_category'] = np.where(species_df['observed_ovules'] + noise <= median_val, "Few", "Many")
    
    # 2. Build out cross-species Dyads (Pairs)
    pairs = []
    for i in range(n_species):
        for j in range(i+1, n_species):
            pairs.append((species_names[i], species_names[j], dist_matrix[i, j]))
            
    dyad_df = pd.DataFrame(pairs, columns=['sp1', 'sp2', 'phylo_distance'])
    
    # Assign MRCA node shortcuts by matching shared parental genus blocks
    sp_to_genus = species_df.set_index('species')['genus'].to_dict()
    dyad_df['mrca_node'] = dyad_df.apply(
        lambda r: f"clade_{sp_to_genus[r['sp1']]}" if sp_to_genus[r['sp1']] == sp_to_genus[r['sp2']] else "deep_root", axis=1
    )
    
    # Map observed traits onto the dyad matrix paths
    trait_dict = species_df.set_index('species')['observed_ovules'].to_dict()
    cat_dict = species_df.set_index('species')['ovule_category'].to_dict()
    
    dyad_df['ovule_number_diff'] = dyad_df.apply(lambda r: abs(trait_dict[r['sp1']] - trait_dict[r['sp2']]), axis=1)
    dyad_df['ovule_category_mismatch'] = dyad_df.apply(lambda r: int(cat_dict[r['sp1']] != cat_dict[r['sp2']]), axis=1)
    
    # Generate response variable: Reproductive Isolation (RI)
    # Isolation falls as distances/disparities widen; inject asymmetric hybrid vigour
    base_ri = 0.8 - (0.01 * dyad_df['ovule_number_diff']) - (0.1 * dyad_df['ovule_category_mismatch']) - (0.05 * dyad_df['phylo_distance'])
    
    is_vigour = np.random.rand(len(dyad_df)) < hybrid_vigour_prop
    dyad_df['reproductive_isolation'] = np.where(is_vigour, np.random.uniform(-0.4, -0.05, len(dyad_df)), base_ri + np.random.normal(0, 0.08, len(dyad_df)))
    
    # Subsample down to target dyad processing size
    sampled_dyads = dyad_df.sample(n=min(n_dyads, len(dyad_df))).reset_index(drop=True)
    
    return sampled_dyads, phylo_cov, species_names

# Execute data generation
my_dyad_data, phylo_covariance, tip_names = simulate_pairwise_dataset(
    n_species=30, n_dyads=160, genus_blur_prop=0.15, few_many_overlap=0.2, hybrid_vigour_prop=0.15
)

# ------------------------------------------------------------------------------
# SECTION 2: MODEL SPECIFICATION & CUSTOM RANDOM EFFECTS INJECTION
# ------------------------------------------------------------------------------
# Format categoricals explicitly as clean strings
my_dyad_data['mrca_node'] = my_dyad_data['mrca_node'].astype(str)
my_dyad_data['sp1'] = my_dyad_data['sp1'].astype(str)
my_dyad_data['sp2'] = my_dyad_data['sp2'].astype(str)
my_dyad_data['ovule_category_mismatch'] = my_dyad_data['ovule_category_mismatch'].astype(str)

# Build baseline regression model via Bambi with a Skew-Normal likelihood architecture
model = bmb.Model(
    "reproductive_isolation ~ phylo_distance + ovule_number_diff + ovule_category_mismatch",
    data=my_dyad_data,
    family="skew_normal"
)

# Incorporate the nested structural grouping adjustments inside the PyMC backend
with model.backend.model:
    # A. Sister-Clade MRCA Random Intercepts
    mrca_idx = pd.Categorical(my_dyad_data['mrca_node']).codes
    n_mrca = len(np.unique(mrca_idx))
    sigma_mrca = pm.HalfNormal("sigma_mrca", sigma=1)
    mrca_offset = pm.Normal("mrca_offset", mu=0, sigma=sigma_mrca, shape=n_mrca)
    
    # B. Multi-Membership Species Intercepts mapped to the Phylogenetic Covariance Matrix
    species_mapping = {name: idx for idx, name in enumerate(tip_names)}
    sp1_idx = my_dyad_data['sp1'].map(species_mapping).values
    sp2_idx = my_dyad_data['sp2'].map(species_mapping).values
    n_species = len(species_mapping)
    
    # Apply Cholesky Factorization to scale the tip random terms exactly by tree distances
    sigma_phylo = pm.HalfNormal("sigma_phylo", sigma=1)
    L_phylo = np.linalg.cholesky(phylo_covariance)
    z_phylo = pm.Normal("z_phylo", shape=n_species)
    phylo_effects = pm.math.matrix_dot(L_phylo, z_phylo) * sigma_phylo
    
    # Apply Multi-Membership averaging across parents (weighted 0.5 each)
    mm_effects = 0.5 * phylo_effects[sp1_idx] + 0.5 * phylo_effects[sp2_idx]
    
    # Push adjustments into the core regression linear predictor vector
    model.backend.model.mu.prior.values += mrca_offset[mrca_idx] + mm_effects

# Sample using multi-core MCMC
idata = model.fit(draws=1500, tune=1000, chains=4, cores=4, target_accept=0.95)

# ------------------------------------------------------------------------------
# SECTION 3: SUMMARY DISPLAY & MARGINAL PLOTS
# ------------------------------------------------------------------------------
summary_df = az.summary(idata, var_names=["~z_phylo"], filter_vars="like")
print("\n--- ESTIMATED POSTERIOR METRICS TABLE ---")
print(summary_df[['mean', 'hdi_3%', 'hdi_97%']].round(3).to_string())

# Plot population scale predictions dropping group factors
sns.set_theme(style="ticks", context="talk")
plt.figure(figsize=(10, 6))

posterior = az.extract(idata)
beta_ovule = posterior["ovule_number_diff"].values
intercept = posterior["Intercept"].values
beta_phylo = posterior["phylo_distance"].values

ovule_grid = np.linspace(my_dyad_data['ovule_number_diff'].min(), my_dyad_data['ovule_number_diff'].max(), 100)
colors = ["#E41A1C", "#377EB8", "#4DAF4A"]

for i, dist_slice in enumerate([0.1, 0.4, 0.8]):
    pred_matrix = intercept[:, None] + (beta_phylo[:, None] * dist_slice) + (beta_ovule[:, None] * ovule_grid)
    mean_pred = pred_matrix.mean(axis=0)
    hdi_lower, hdi_upper = np.percentile(pred_matrix, [3, 97], axis=0)
    
    plt.plot(ovule_grid, mean_pred, label=f"Tree Distance Slice: {dist_slice}", color=colors[i], linewidth=2.5)
    plt.fill_between(ovule_grid, hdi_lower, hdi_upper, color=colors[i], alpha=0.15)

plt.title("Marginal Effect of Ovule Disparity on Simulated RI", weight='bold', pad=15)
plt.xlabel("Phenotypic Disparity (Absolute Ovule Number Difference)")
plt.ylabel("Predicted Reproductive Isolation (RI)")
plt.legend(frameon=False, loc="upper right")
sns.despine()
plt.tight_layout()
plt.show()


=========================================================================================
Bayesian Phylogenetic Mixed Model (PMM) Sensitivity and Estimation Framework
=========================================================================================
Purpose: 
  This production-grade script provides a complete simulation-to-estimation tool 
  to analyze continuous and categorical trait drivers of Reproductive Isolation (RI).

Methodological Architecture:
  1. Systematic Pairwise Simulator: Evaluates data degradation (genus imputation, 
     categorical assignment boundary blur) and a skewed hybrid vigour phenotypic tail.
  2. Custom PyMC/Bambi Engine: Uses a Skew-Normal Likelihood to model asymmetric 
     distributions safely, combined with custom multi-membership Cholesky-factored 
     phylogenetic tree matrices and nested sister-clade intercepts.
  3. Auto-Validation: Verifies that local packages match your environment configuration."""

import sysimport subprocessimport pkg_resources
# =======================================================================================# SECTION 0: AUTOMATED RUNTIME CONTEXT & PIP CONFIGURATION VALIDATION# =======================================================================================
def verify_environment_requirements(requirements_path="requirements.txt"):
  """
    Validates package compliance prior to memory allocation or C-compilation.
    
    If libraries are missing or version-conflicted, this system shells out to pip
    to automatically synchronize your environment state.
    """
print(f"[*] Parsing system dependencies using configuration: '{requirements_path}'...")
try:
  with open(requirements_path, "r") as f:
  # Parse individual rows from the targets text file
  requirements = pkg_resources.parse_requirements(f.readlines())

missing_packages = []
for req in requirements:
  try:
  # Assert package exists and matches specific structural constraints
  pkg_resources.require(str(req))
except (pkg_resources.DistributionNotFound, pkg_resources.VersionConflict):
  missing_packages.append(str(req))

# Repair the active python terminal sandbox state if gaps are found
if missing_packages:
  print(f"[!] Warning: Missing or version-conflicted dependencies detected: {missing_packages}")
print("[*] Automatically executing pip installer to synchronize your environment...")
subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", requirements_path])
print("[+] Environment successfully repaired and synchronized.\n")
else:
  print("[+] Environment check passed. All verified packages are ready.\n")

except FileNotFoundError:
  print(f"[!] Configuration Warning: '{requirements_path}' was not found in this directory.")
print("[*] Continuing runtime attempt using standard workspace packages...\n")
# Run environmental synchronization check
verify_environment_requirements("requirements.txt")
# Core array, indexing, and visualization modulesimport numpy as npimport pandas as pdimport matplotlib.pyplot as pltimport seaborn as sns
# Advanced Bayesian mathematical modeling and posterior parsing toolsimport pymc as pmimport bambi as bmbimport arviz as az
# Tree clustering and matrix generation utilitiesfrom scipy.cluster.hierarchy import linkage, fclusterfrom scipy.spatial.distance import squareform

# =======================================================================================# SECTION 1: SYSTEMATIC PAIRWISE DATA SIMULATOR# =======================================================================================
def simulate_pairwise_dataset(
  n_species=30, 
  n_dyads=150, 
  genus_blur_prop=0.1, 
  few_many_overlap=0.2, 
  hybrid_vigour_prop=0.15
):
  """
    Simulates a macroevolutionary dataset for cross-species pairs (dyads).
    
    Parameters:
      n_species (int): Count of distinct operational taxonomic species in your pool.
      n_dyads (int): Number of species pairs to subsample for model fitting.
      genus_blur_prop (float): Prop of individual counts replaced with coarse genus means.
      few_many_overlap (float): Multiplier controlling overlap when categorizing traits.
      hybrid_vigour_prop (float): Prop of pairs that break constraints, sliding below 0.
      
    Returns:
      sampled_dyads (DataFrame): Fully constructed training observations dataset.
      phylo_cov (ndarray): Matrix capturing expected baseline genetic similarity.
      species_names (list): Ordered strings identifying active species tips.
    """
np.random.seed(42)
species_names = [f"species_{i}" for i in range(n_species)]

# -----------------------------------------------------------------------------------
# Step A: Simulate Tree Spatial Properties & Continuous Phenotype Range
# -----------------------------------------------------------------------------------
# Map raw 2D random coords to serve as spatial branch evolution roots
raw_coords = np.random.rand(n_species, 2)
# Calculate pairwise distance matrices across coordinates
dist_matrix = squareform(np.linalg.norm(raw_coords[:, None] - raw_coords, axis=2))

# Convert tree distances to a relative covariance matrix (Max distance minus distance)
max_d = dist_matrix.max()
phylo_cov = max_d - dist_matrix
# Add a small numerical jitter along the diagonal to ensure matrix invertibility
np.fill_diagonal(phylo_cov, phylo_cov.diagonal() + 1e-4)

# Simulate true uncorrupted underlying species ovule parameters using tree covariance
true_ovules_raw = np.random.multivariate_normal(np.zeros(n_species), phylo_cov)
# Exponentiate and rescale to match realistic positive biological counts (Mean ~20 ovules)
true_ovules = np.round(np.exp((true_ovules_raw - true_ovules_raw.mean()) / true_ovules_raw.std() * 0.5 + 3))

# -----------------------------------------------------------------------------------
# Step B: Inject Genus Imputation and Categorical Overlap Degradation
# -----------------------------------------------------------------------------------
# Group neighboring species tips into shared "Genus" buckets via hierarchical tree splitting
Z = linkage(dist_matrix, 'complete')
genus_ids = fcluster(Z, t=max(3, n_species // 5), criterion='maxclust')

# Bundle species metrics into a tracking frame
species_df = pd.DataFrame({
  'species': species_names,
  'true_ovules': true_ovules,
  'genus': genus_ids
})

# Inject Imputation Blur: Replace individual values with genus-wide expectations
genus_means = species_df.groupby('genus')['true_ovules'].transform('mean').round()
mask = np.random.rand(n_species) < genus_blur_prop
species_df['observed_ovules'] = np.where(mask, genus_means, species_df['true_ovules'])

# Apply Categorical Disruption: Split values into "Few" vs "Many" bins around the median,
# adding Gaussian noise to simulate blurry classification thresholds
median_val = species_df['observed_ovules'].median()
noise = np.random.normal(0, few_many_overlap * median_val, n_species)
species_df['ovule_category'] = np.where(species_df['observed_ovules'] + noise <= median_val, "Few", "Many")

# -----------------------------------------------------------------------------------
# Step C: Construct Cross-Species Pairwise Combinations (Dyads)
# -----------------------------------------------------------------------------------
pairs = []
for i in range(n_species):
  for j in range(i+1, n_species):
  pairs.append((species_names[i], species_names[j], dist_matrix[i, j]))

dyad_df = pd.DataFrame(pairs, columns=['sp1', 'sp2', 'phylo_distance'])

# Extract Most Recent Common Ancestor (MRCA) nodes by checking shared genus buckets
sp_to_genus = species_df.set_index('species')['genus'].to_dict()
dyad_df['mrca_node'] = dyad_df.apply(
  lambda r: f"clade_{sp_to_genus[r['sp1']]}" if sp_to_genus[r['sp1']] == sp_to_genus[r['sp2']] else "deep_root", axis=1
)

# Map parents' traits back onto your dyad tracking columns
trait_dict = species_df.set_index('species')['observed_ovules'].to_dict()
cat_dict = species_df.set_index('species')['ovule_category'].to_dict()

# Continuous disparity = Absolute structural count difference
dyad_df['ovule_number_diff'] = dyad_df.apply(lambda r: abs(trait_dict[r['sp1']] - trait_dict[r['sp2']]), axis=1)
# Categorical mismatch = Binary indicator flag (0 = same category, 1 = mismatch)
dyad_df['ovule_category_mismatch'] = dyad_df.apply(lambda r: int(cat_dict[r['sp1']] != cat_dict[r['sp2']]), axis=1)

# -----------------------------------------------------------------------------------
# Step D: Structural Response Formulation & Hybrid Vigour Tail Shaping
# -----------------------------------------------------------------------------------
# Establish baseline expected linear function effects (RI decreases as traits vary)
base_ri = 0.8 - (0.01 * dyad_df['ovule_number_diff']) - (0.1 * dyad_df['ovule_category_mismatch']) - (0.05 * dyad_df['phylo_distance'])

# Inject Hybrid Vigour: Drive a proportion of points below 0 using a uniform block distribution
is_vigour = np.random.rand(len(dyad_df)) < hybrid_vigour_prop
dyad_df['reproductive_isolation'] = np.where(
  is_vigour, 
  np.random.uniform(-0.4, -0.05, len(dyad_df)), # Unbounded negative values
  base_ri + np.random.normal(0, 0.08, len(dyad_df)) # Bounded values with noise
)

# Randomly subsample down to your target testing size
sampled_dyads = dyad_df.sample(n=min(n_dyads, len(dyad_df))).reset_index(drop=True)

return sampled_dyads, phylo_cov, species_names
## Run the data generator
my_dyad_data, phylo_covariance, tip_names = simulate_pairwise_dataset(
  n_species=30, n_dyads=160, genus_blur_prop=0.15, few_many_overlap=0.2, hybrid_vigour_prop=0.15
)
## =======================================================================================## SECTION 2: MODEL SPECIFICATION & CUSTOM RANDOM EFFECTS INJECTION## =======================================================================================## Clean categorical indices to prevent internal formula evaluation crashes
my_dyad_data['mrca_node'] = my_dyad_data['mrca_node'].astype(str)
my_dyad_data['sp1'] = my_dyad_data['sp1'].astype(str)
my_dyad_data['sp2'] = my_dyad_data['sp2'].astype(str)
my_dyad_data['ovule_category_mismatch'] = my_dyad_data['ovule_category_mismatch'].astype(str)
## Build baseline structural linear model using Bambi.## Family 'skew_normal' introduces a third parameter (alpha) to capture the hybrid vigour tail.
model = bmb.Model(
  "reproductive_isolation ~ phylo_distance + ovule_number_diff + ovule_category_mismatch",
  data=my_dyad_data,
  family="skew_normal"
)
## Intercept and append complex phylogenetic matrix operations inside the PyMC backend context
with model.backend.model:
  # -----------------------------------------------------------------------------------
# Term A: Sister-Clade MRCA Random Intercept Block
# -----------------------------------------------------------------------------------
# Map MRCA node strings into simple integer codes
mrca_idx = pd.Categorical(my_dyad_data['mrca_node']).codes
n_mrca = len(np.unique(mrca_idx))
# Define hyperprior variance and normal parameter offsets for sister clades
sigma_mrca = pm.HalfNormal("sigma_mrca", sigma=1)
mrca_offset = pm.Normal("mrca_offset", mu=0, sigma=sigma_mrca, shape=n_mrca)
# -----------------------------------------------------------------------------------
# Term B: Multi-Membership Tree-Correlated Tip Error Block
# -----------------------------------------------------------------------------------
# Map species names to unique index values
species_mapping = {name: idx for idx, name in enumerate(tip_names)}
sp1_idx = my_dyad_data['sp1'].map(species_mapping).values
sp2_idx = my_dyad_data['sp2'].map(species_mapping).values
n_species = len(species_mapping)
# Hyperprior tracking scale parameter
sigma_phylo = pm.HalfNormal("sigma_phylo", sigma=1)
# Apply Cholesky Factorization to scale the tip random terms exactly by tree distances
L_phylo = np.linalg.cholesky(phylo_covariance)
z_phylo = pm.Normal("z_phylo", shape=n_species)
# Matrix multiply Cholesky triangle with standard normal vector to apply tree correlation
phylo_effects = pm.math.matrix_dot(L_phylo, z_phylo) * sigma_phylo
# Multi-Membership Averaging: Give each parent species an equal weight of 0.5
mm_effects = 0.5 * phylo_effects[sp1_idx] + 0.5 * phylo_effects[sp2_idx]
# -----------------------------------------------------------------------------------
# Term C: Linear Predictor Matrix Injection
# -----------------------------------------------------------------------------------
# Add your custom group-level terms to Bambi's baseline regression model structure
model.backend.model.mu.prior.values += mrca_offset[mrca_idx] + mm_effects
print("[*] Compiling model functions and initiating MCMC sampling chains...")
## Run multi-chain sampler across 4 CPU cores
idata = model.fit(draws=1500, tune=1000, chains=4, cores=4, target_accept=0.95)
## =======================================================================================## SECTION 3: SUMMARY DISPLAY & POPULATION MARGINAL PLOTS## =======================================================================================## Extract summary table metrics using ArviZ, filtering out internal tracking dimensions
summary_df = az.summary(idata, var_names=["~z_phylo"], filter_vars="like")
print("\n--- ESTIMATED POSTERIOR METRICS TABLE ---")
print(summary_df[['mean', 'hdi_3%', 'hdi_97%']].round(3).to_string())
## Configure plot styling
sns.set_theme(style="ticks", context="talk")
plt.figure(figsize=(10, 6))
## Extract the fixed effect slopes from your MCMC chains
posterior = az.extract(idata)
beta_ovule = posterior["ovule_number_diff"].values
intercept = posterior["Intercept"].values
beta_phylo = posterior["phylo_distance"].values
## Create an evenly-spaced prediction grid across your empirical range of ovule differences
ovule_grid = np.linspace(my_dyad_data['ovule_number_diff'].min(), my_dyad_data['ovule_number_diff'].max(), 100)
colors = ["#E41A1C", "#377EB8", "#4DAF4A"]
## Calculate and plot expected population-scale regression slopes at different tree distances
for i, dist_slice in enumerate([0.1, 0.4, 0.8]):
  # Compute population predictions for each MCMC sample, dropping group-level terms
  pred_matrix = intercept[:, None] + (beta_phylo[:, None] * dist_slice) + (beta_ovule[:, None] * ovule_grid)
# Calculate the mean expected trend across samples
mean_pred = pred_matrix.mean(axis=0)
# Calculate the 95% Bayesian Highest Density Interval (HDI) bounds
hdi_lower = np.percentile(pred_matrix, 3, axis=0)
hdi_upper = np.percentile(pred_matrix, 97, axis=0)
# Plot expected conditional line
plt.plot(ovule_grid, mean_pred, label=f"Tree Distance Slice: {dist_slice}", color=colors[i], linewidth=2.5)
# Add shaded uncertainty bands around the regression lines
plt.fill_between(ovule_grid, hdi_lower, hdi_upper, color=colors[i], alpha=0.15)
## Add plot labels and clean layout aesthetics
plt.title("Marginal Effect of Ovule Disparity on Simulated RI", weight='bold', pad=15)
plt.xlabel("Phenotypic Disparity (Absolute Ovule Number Difference)")
plt.ylabel("Predicted Reproductive Isolation (RI)")
plt.legend(frameon=False, loc="upper right")
sns.despine()
plt.tight_layout()
plt.show()
  
  
  
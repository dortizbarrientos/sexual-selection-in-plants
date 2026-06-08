library(tidyverse)
library(ape)
library(nlme)
library(phytools)

table1 <- read.csv("Table1_2025_09_18.csv")
table2 <- read.csv("Table2_2025_12_15.csv")
table2b <- read.csv("Table2b_2025_11_25.csv")
gen_dist <- read.csv("Genetic_distances_results.csv")
gen_dist <- distinct(gen_dist)
#phy <- read.tree("All_PI_pairs_v3_2025_08.tre.txt")
phy <- read.tree("All_PI_pairs_Zuntini_2025_10.tre") # see below for code to produce this

# look at the data
plot(table(table2$life_history, table2$mating_system))

# clean up and summarize ovule number Lists
# remove duplicate/problematic values and find mean ovule number per species
species_ovule_numbers <- table2b %>% 
  filter(Remove!="yes") %>% # these are mostly cases where a value in cited in later papers - see googlesheet for details
  group_by(species) 

# plot ovule number distribution - note that some species are represented multiple times
hist(species_ovule_numbers$number, xlim=c(0,200), breaks=100000)
hist(species_ovule_numbers$number, xlim=c(0,50), breaks=1000000)

# find mean and median ovule number for each genus
# think about how to deal with genus-level values
genus_ovule_numbers <- species_ovule_numbers %>% 
  summarise(genus = first(genus), level = first(level), mean = mean(number, na.rm=T), count = n()) %>% # find average for each species
  group_by(genus) %>%
  summarise(genus_mean_ovule_number = mean(mean, na.rm=T), genus_median_ovule_number = median(mean, na.rm=T), count = n())

# plot genus-level ovule numbers
hist(log(genus_ovule_numbers$genus_mean_ovule_number))
hist(log(genus_ovule_numbers$genus_median_ovule_number))
plot(log(genus_ovule_numbers$genus_mean_ovule_number), log(genus_ovule_numbers$genus_median_ovule_number))

# add genus-level ovule number to table 2 
table2_plus <- left_join(table2, genus_ovule_numbers, by = "genus")

# compare genus-level averages to species-specific ovule numbers
# see googlesheet for a full explanation of confidence, but 1 = species-specific ovule number estimate, 5 = no data for any species in the genus
plot(log(table2_plus$ovule_number), log(table2_plus$genus_mean_ovule_number), col=table2_plus$confidence)
text(1:3, rep(10,3), c("1", "2", "3"), col=1:3)
plot(log(table2_plus$ovule_number), log(table2_plus$genus_median_ovule_number), col=table2_plus$confidence)
text(1:3, rep(10,3), c("1", "2", "3"), col=1:3)

# use genus-level ovule data to replace other values 
table2_plus <- table2_plus %>%
  # can edit below line to use mean or replace more or fewer ovule numbers
  mutate(ovule_number_plus = if_else(confidence > 3, genus_median_ovule_number, ovule_number))
  
# find genera with very little data
table2_plus %>% filter(count<3) %>% select(genus) %>% unique() # could look for more species in these genera
table2_plus %>% filter(is.na(count)) %>% select(genus) %>% unique()

# find family level values
family_ovule_numbers <- species_ovule_numbers %>% 
  group_by(family) %>%
  summarise(family_mean_ovule_number = mean(number, na.rm=T), count = n())

# add values for species with no data (see Missing_genera tab for more info)
table2_plus[table2_plus$genus == "Sinalliaria",]$ovule_number_plus <- 22
table2_plus[table2_plus$genus == "Achimenes",]$ovule_number_plus <- 1409.428571
table2_plus[table2_plus$genus == "Diplacus",]$ovule_number_plus <- 549.872222
table2_plus$ovule_number_plus[table2_plus$genus %in% c(
  "Cattleya", "Chiloglottis", "Dendrobium", "Habenaria", "Serapias", "Gymnadenia", "Satyrium"
)] <- 21036.52664

# look at distribution and generate categorical data
hist(log(table2_plus$ovule_number_plus), breaks = 50); abline(v=2, col="red", lty=2)
table2_plus <- table2_plus %>% mutate(ov_cat = case_when(log(ovule_number_plus) < 2 ~ "few", .default = "many"))

# select columns for final data about each species 
table2_final <- table2_plus %>% select(species, ovule_number_plus, ov_cat, life_history, mating_system)

# add table 2 and genetic distance data to table 1 
table1_plus <- left_join(table1, table2_final, by = c("species1" = "species")) %>% 
  left_join(., table2_final, by = c("species2" = "species")) %>%
  left_join(., gen_dist, by = c("species1" = "Species_1", "species2" = "Species_2"))
  
# summarize table 1 data for basic analyses 
table1_plus2 <- table1_plus %>%  mutate(mean_ov = (ovule_number_plus.x+ovule_number_plus.y)/2,
         min_ov = pmin(ovule_number_plus.x,ovule_number_plus.y),
         mating_system = case_when(mating_system.x=="outcrossing" & mating_system.y=="outcrossing" ~ "outcrossing",
                                   mating_system.x=="selfing" & mating_system.y=="selfing" ~ "selfing",
                                   .default = "mixed-mating"),
         life_history = case_when(life_history.x=="woody" & life_history.y=="woody" ~ "slow",
                                  life_history.x=="annual" & life_history.y=="annual" ~ "fast",
                                  .default = "intermediate"),
         ov_cat = case_when(ov_cat.x=="few" & ov_cat.y=="few" ~ "few",
                            ov_cat.x=="many" & ov_cat.y=="many" ~ "many",
                                  .default = "intermediate"),
         # added this category because original life histor category was giving weird resuls where "intermediate" had lower RI than both other types
         life_history2 = case_when(life_history.x=="woody" & life_history.y=="woody" ~ "woody",
                                   life_history.x=="woody" & life_history.y=="perennial" | life_history.x=="perennial" & life_history.y=="woody" ~ "WP",
                                   life_history.x=="perennial" & life_history.y=="perennial" ~ "perennial",
                                   life_history.x=="annual" & life_history.y=="perennial" | life_history.x=="perennial" & life_history.y=="annual" ~ "AP",
                                   life_history.x=="annual" & life_history.y=="annual" ~ "annual",
                                   .default = "mixed"))
                        
                              

# Look at relationships with genetic distance
hist(table1_plus2$K2P_Distance, breaks=25)
stripchart(table1_plus2$K2P_Distance ~ table1_plus2$mating_system, vertical=T, method="jitter")
stripchart(table1_plus2$K2P_Distance ~ table1_plus2$life_history, vertical=T, method="jitter")
stripchart(table1_plus2$K2P_Distance ~ table1_plus2$life_history2, vertical=T, method="jitter")

# replace all missing K2P with mean value to retain the 21 species pairs withour genetic distance 
# table1_plus2$K2P_Distance[is.na(table1_plus2$K2P_Distance)] <- mean(table1_plus2$K2P_Distance, na.rm=T)

# filter table 1 for basic analyses 
# note that the code we are using to run PGLS under BM cannot handle missing data 
table1_filt <- table1_plus2 %>% 
  filter(remove!="yes") %>% # this removes duplicates for pairs with multiple measures - keeps the best value (usually pollen competition), but we could try averaging these 
  filter(!is.na(PI_name)) %>% # this keeps only the phylogenetically independent pairs 
  filter(!is.na(mean_ov)) %>% # this remove pairs without ovule data (there should be none)
  filter(!is.na(K2P_Distance)) %>% # this removes pairs without genetic distance (currently 21 species pairs)
  #filter(K2P_Distance<0.2) %>%
  select(PI_name, RI, measure, min_ov, mean_ov, ov_cat, mating_system, life_history, life_history2, K2P_Distance)
names(table1_filt)[1] <- "species"
table1_filt$measure <- as.factor(table1_filt$measure)
table1_filt$mating_system <- as.factor(table1_filt$mating_system)
table1_filt$life_history <- as.factor(table1_filt$life_history)
table1_filt$ov_cat <- as.factor(table1_filt$ov_cat)

# plot relationships with RI
hist(table1_filt$RI)
cols <- colorRampPalette(c("blue", "red"))(100)
plot(table1_filt$K2P_Distance, table1_filt$RI, col=cols[cut(log(table1_filt$mean_ov), breaks = 100)]); text(0.25, c(0.95,0.9), c("many", "few"), col=c("red", "blue"))
plot(log(table1_filt$mean_ov), table1_filt$RI, col=cols[cut(table1_filt$K2P_Distance, breaks = 100)]); text(10, c(0.95,0.9), c("distant", "close"), col=c("red", "blue"))
stripchart(table1_filt$RI ~ table1_filt$mating_system, vertical=T, method="jitter")
stripchart(table1_filt$RI ~ table1_filt$life_history, vertical=T, method="jitter")
stripchart(table1_filt$RI ~ table1_filt$life_history2, vertical=T, method="jitter")
stripchart(table1_filt$RI ~ table1_filt$measure, vertical=T, method="jitter")
stripchart(table1_filt$RI ~ table1_filt$ov_cat, vertical=T, method="jitter")

# format data for PGLS
# drop tips from tre file that aren't in table1_filt
phy <- read.tree("All_PI_pairs_Zuntini_2025_10.tre") # see below for code to produce this
phy  <- drop.tip(phy, setdiff(phy$tip.label, table1_filt$species))
plot(phy)
# remove lines in table1_filt that aren't in tre file
table1_filt <- subset(table1_filt, species %in% phy$tip.label)
rownames(table1_filt) <- table1_filt$species

# run PGLS under BM

# method 1 - keep all covariates
# all covariates
fit_bm_4 <- gls(RI ~ K2P_Distance + mating_system + life_history + measure, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_4)
# full model with ovule number
fit_bm_5 <- gls(RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history + measure, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_5); anova(fit_bm_4, fit_bm_5)
# full model with ovule category
fit_bm_5b <- gls(RI ~ ov_cat + K2P_Distance + mating_system + life_history + measure, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_5b); anova(fit_bm_4, fit_bm_5b)

# method 2 - build up best model
# models with only one predictor - K2P is best based on AIC
fit_bm_1_gd <- gls(RI ~ K2P_Distance, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML") 
fit_bm_1_ms <- gls(RI ~ mating_system, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
fit_bm_1_lh <- gls(RI ~ life_history, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
fit_bm_1_m <- gls(RI ~ measure, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
fit_bm_1_ov <- gls(RI ~ log(mean_ov), data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML") 
anova(fit_bm_1_gd, fit_bm_1_ms, fit_bm_1_lh, fit_bm_1_m, fit_bm_1_ov)
# add predictors to K2P - measure is best and a marginally significant improvement
fit_bm_2_ms <- gls(RI ~ K2P_Distance + mating_system, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML") 
fit_bm_2_lh <- gls(RI ~ K2P_Distance + life_history, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML") 
fit_bm_2_m <- gls(RI ~ K2P_Distance + measure, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML") 
fit_bm_2_ov <- gls(RI ~ K2P_Distance + log(mean_ov), data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML") 
anova(fit_bm_1_gd, fit_bm_2_ms, fit_bm_2_lh, fit_bm_2_m, fit_bm_2_ov); anova(fit_bm_1_gd, fit_bm_2_m)
# add predictors to K2P + measure - ovule number is best but not a significant improvement
fit_bm_3_ms <- gls(RI ~ K2P_Distance + measure + mating_system, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML") 
fit_bm_3_lh <- gls(RI ~ K2P_Distance + measure + life_history, data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
fit_bm_3_ov <- gls(RI ~ K2P_Distance + measure + log(mean_ov), data = table1_filt, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
anova(fit_bm_2_m, fit_bm_3_ms, fit_bm_3_lh, fit_bm_3_ov); anova(fit_bm_2_m, fit_bm_3_ov)

# do PGLS using pollen competition only 
table1_pollen_comp <- table1_plus2 %>% 
  filter(measure=="pollen competition") %>%
  filter(!is.na(PI_name)) %>% # this keeps only the phylogenetically independent pairs 
  filter(!is.na(mean_ov)) %>% # this remove pairs without ovule data (there should be none)
  filter(!is.na(K2P_Distance)) %>% # this removes pairs without genetic distance (currently 21 species pairs)
  select(PI_name, RI, min_ov, mean_ov, mating_system, life_history, K2P_Distance)
names(table1_pollen_comp)[1] <- "species"
table1_pollen_comp$mating_system <- as.factor(table1_pollen_comp$mating_system)
table1_pollen_comp$life_history <- as.factor(table1_pollen_comp$life_history)

# format data for PGLS
phy <- read.tree("All_PI_pairs_Zuntini_2025_10.tre") # see below for code to produce this
phy  <- drop.tip(phy, setdiff(phy$tip.label, table1_pollen_comp$species))
# remove lines in table1_filt that aren't in tre file
table1_pollen_comp <- subset(table1_pollen_comp, species %in% phy$tip.label)
rownames(table1_pollen_comp) <- table1_pollen_comp$species

# method 1 - keep all covariates
# all covariates
fit_bm_pc3 <- gls(RI ~ K2P_Distance + mating_system + life_history, data = table1_pollen_comp, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_pc3)
# full model with ovule number
fit_bm_pc4 <- gls(RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history , data = table1_pollen_comp, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_pc4); anova(fit_bm_pc3, fit_bm_pc4)

plot(log(table1_pollen_comp$mean_ov), table1_pollen_comp$RI, col=cols[cut(table1_pollen_comp$K2P_Distance, breaks = 100)]); text(10, c(0.95,0.9), c("distant", "close"), col=c("red", "blue"))

# do PGLS using seed set only 
table1_seed_set <- table1_plus2 %>% 
  filter(measure=="seed set") %>%
  filter(!is.na(PI_name)) %>% # this keeps only the phylogenetically independent pairs 
  filter(!is.na(mean_ov)) %>% # this remove pairs without ovule data (there should be none)
  filter(!is.na(K2P_Distance)) %>% # this removes pairs without genetic distance (currently 21 species pairs)
  select(PI_name, RI, min_ov, mean_ov, mating_system, life_history, K2P_Distance)
names(table1_seed_set)[1] <- "species"
table1_seed_set$mating_system <- as.factor(table1_seed_set$mating_system)
table1_seed_set$life_history <- as.factor(table1_seed_set$life_history)

# format data for PGLS
phy <- read.tree("All_PI_pairs_Zuntini_2025_10.tre") # see below for code to produce this
phy  <- drop.tip(phy, setdiff(phy$tip.label, table1_seed_set$species))
# remove lines in table1_filt that aren't in tre file
table1_seed_set <- subset(table1_seed_set, species %in% phy$tip.label)
rownames(table1_seed_set) <- table1_seed_set$species

# method 1 - keep all covariates
# all covariates
fit_bm_ss3 <- gls(RI ~ K2P_Distance + mating_system + life_history, data = table1_seed_set, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_ss3)
# full model with ovule number
fit_bm_ss4 <- gls(RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history , data = table1_seed_set, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_ss4); anova(fit_bm_ss3, fit_bm_ss4)

plot(log(table1_seed_set$mean_ov), table1_seed_set$RI, col=cols[cut(table1_seed_set$K2P_Distance, breaks = 100)]); text(10, c(0.95,0.9), c("distant", "close"), col=c("red", "blue"))


# run other PGLS options - can't have branch lengths of 0
# model = "BM", "OUrandomRoot", "lambda", "delta", "kappa", etc.
library(phylolm)
phy <- read.tree("All_PI_pairs_Zuntini_2025_10.tre")
# code to change 0 branch lengths to 0.00001
hist(phy$edge.length)
phy$edge.length[phy$edge.length==0] <- 0.00001
plot(phy)

# model type matters a lot!
# *** we should look into this ***
summary(phylolm(RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history, data = table1_filt, phy = phy, model = "BM"))
summary(phylolm(RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history, data = table1_filt, phy = phy, model = "OUrandomRoot"))
summary(phylolm(RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history, data = table1_filt, phy = phy, model = "lambda"))
summary(phylolm(RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history, data = table1_filt, phy = phy, model = "delta"))
summary(phylolm(RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history, data = table1_filt, phy = phy, model = "kappa"))


# look at asymmetrical RI

# summarize table 1 data for asymmetry 
table1_plus3 <- table1_plus %>% mutate(RI_sp1_sp2 = 1-(2*interspecific_value_sp1_sp2/(interspecific_value_sp1_sp2+intraspecific_value_sp1)),
                                       RI_sp2_sp1 = 1-(2*interspecific_value_sp2_sp1/(interspecific_value_sp2_sp1+intraspecific_value_sp2)),
                                       RI_diff = RI_sp1_sp2-RI_sp2_sp1,
                                       ov_diff = ovule_number_plus.x-ovule_number_plus.y,
                                       mean_ov = (ovule_number_plus.x+ovule_number_plus.y)/2,
                                       mating_system = case_when(mating_system.x=="outcrossing" & mating_system.y=="outcrossing" ~ "O-O",
                                                                 mating_system.x=="outcrossing" & mating_system.y=="mixed-mating" ~ "O-M",
                                                                 mating_system.x=="outcrossing" & mating_system.y=="selfing" ~ "O-S",
                                                                 mating_system.x=="selfing" & mating_system.y=="selfing" ~ "S-S",
                                                                 mating_system.x=="selfing" & mating_system.y=="mixed-mating" ~ "S-M",
                                                                 mating_system.x=="selfing" & mating_system.y=="outcrossing" ~ "S-O",
                                                                 mating_system.x=="mixed-mating" & mating_system.y=="mixed-mating" ~ "M-M",
                                                                 mating_system.x=="mixed-mating" & mating_system.y=="outcrossing" ~ "M-O",
                                                                 mating_system.x=="mixed-mating" & mating_system.y=="selfing" ~ "M-S",
                                                                 .default = "other"),
                                       life_history = case_when(life_history.x=="woody" & life_history.y=="woody" ~ "W-W",
                                                                life_history.x=="woody" & life_history.y=="perennial" ~ "W-P",
                                                                life_history.x=="woody" & life_history.y=="annual" ~ "W-A",
                                                                life_history.x=="annual" & life_history.y=="annual" ~ "A-A",
                                                                life_history.x=="annual" & life_history.y=="perennial" ~ "A-P",
                                                                life_history.x=="annual" & life_history.y=="woody" ~ "A-W",
                                                                life_history.x=="perennial" & life_history.y=="perennial" ~ "P-P",
                                                                life_history.x=="perennial" & life_history.y=="woody" ~ "P-W",
                                                                life_history.x=="perennial" & life_history.y=="annual" ~ "P-A",
                                                                .default = "other"))

# filter data
table1_plus3 <- table1_plus3 %>% 
  filter(remove!="yes") %>% # this removes duplicates for pairs with multiple measures - keeps the best value (usually pollen competition), but we could try averaging these 
  filter(!is.na(PI_name)) %>% # this keeps only the phylogenetically independent pairs 
  filter(!is.na(K2P_Distance)) %>%
  filter(!is.na(RI_sp1_sp2)) %>%
  filter(!is.na(RI_sp2_sp1))

# look are the relationship between different directions of RI
hist(table1_plus3$RI_sp1_sp2)
hist(table1_plus3$RI_sp2_sp1)
plot(table1_plus3$RI_sp1_sp2, table1_plus3$RI_sp2_sp1)

# look at distributions of differences
hist(table1_plus3$RI_diff)
hist(table1_plus3$ov_diff)

# look at relationship between RI diff and various things
plot(table1_plus3$ov_diff, table1_plus3$RI_diff)
plot(log(table1_plus3$mean_ov), table1_plus3$RI_diff)
stripchart(table1_plus3$RI_diff ~ table1_plus3$mating_system, vertical=T, method="jitter")
stripchart(table1_plus3$RI_diff ~ table1_plus3$life_history, vertical=T, method="jitter")

# look at relationships between maternal ovule number and RI into mom
plot(log(table1_plus3$ovule_number_plus.x), table1_plus3$RI_sp1_sp2)
plot(log(table1_plus3$ovule_number_plus.y), table1_plus3$RI_sp2_sp1)

# format data for PGLS
# drop tips from tre file that aren't in table1_filt
phy <- read.tree("All_PI_pairs_Zuntini_2025_10.tre") # see below for code to produce this
phy  <- drop.tip(phy, setdiff(phy$tip.label, table1_plus3$PI_name))
plot(phy)
# remove lines in table1_filt that aren't in tre file
table1_plus3 <- subset(table1_plus3, PI_name %in% phy$tip.label)
rownames(table1_plus3) <- table1_plus3$PI_name

# run PGLS under BM for sp1-sp2
# all covariates
fit_bm_4 <- gls(RI_sp1_sp2 ~ K2P_Distance + mating_system + life_history + measure, data = table1_plus3, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
fit_bm_5 <- gls(RI_sp1_sp2 ~ log(ovule_number_plus.x) + K2P_Distance + mating_system + life_history + measure, data = table1_plus3, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_5); anova(fit_bm_4, fit_bm_5)

# run PGLS under BM for sp2-sp1
# all covariates
fit_bm_4 <- gls(RI_sp2_sp1 ~ K2P_Distance + mating_system + life_history + measure, data = table1_plus3, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
fit_bm_5 <- gls(RI_sp2_sp1 ~ log(ovule_number_plus.y) + K2P_Distance + mating_system + life_history + measure, data = table1_plus3, correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_5); anova(fit_bm_4, fit_bm_5)



### --- Make tre file with branch lengths! --- ###

Zuntini_trees <- read.tree("Zuntini_etal_2024.tre.txt")
Zuntini_tree <- Zuntini_trees[[1]]

# pull out tips within our focal genera 
Zuntini_tips <- print(Zuntini_tree$tip.label)
genera <- read.table("PI_pairs_genera.txt") %>% unique() %>% unlist()
matches <- Zuntini_tips[grepl(paste(genera, collapse = "|"), Zuntini_tips, ignore.case = TRUE)]
# write and edit csv to match original tip names to PI pair names 
write.csv(matches, "Zuntini_tip_matches.csv")
table1$PI_name %>% unique()
# Genera missing from Zuntini and the closest genus (in partheneses) in the tree according to Tree of Life
# Zuntini_tips[grepl("Lathraea", Zuntini_tips, ignore.case = TRUE)]
# Centaurium (Chironia), Erica (Calluna), Kalmia (Rhodothamnus), Mussaenda (Heinsia), Rhinanthus (Lathraea)

# prune tree to the tips we want to include
matches_list <- read.csv("Zuntini_tip_matches_edited.csv") %>% filter(complete.cases(PI_pair))
pruned_Zuntini_tree<-drop.tip(Zuntini_tree,Zuntini_tree$tip.label[-match(matches_list$zuntini_tip, Zuntini_tree$tip.label)])

# change tip names  
Zuntini_pruned_tips <- print(pruned_Zuntini_tree$tip.label) %>% as.data.frame()
temp <- data.frame(zuntini_tip = Zuntini_pruned_tips)
names(temp) <- "zuntini_tip"
temp2 <- left_join(temp, matches_list)
pruned_Zuntini_tree$tip.label <- temp2$PI_pair
plot(pruned_Zuntini_tree)

# add additional pairs as tips in the same genus
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Aquilegia02", where = which(pruned_Zuntini_tree$tip.label == "Aquilegia01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Aquilegia03", where = which(pruned_Zuntini_tree$tip.label == "Aquilegia01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Clarkia02", where = which(pruned_Zuntini_tree$tip.label == "Clarkia01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Clarkia03", where = which(pruned_Zuntini_tree$tip.label == "Clarkia01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Costus02", where = which(pruned_Zuntini_tree$tip.label == "Costus01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Costus03", where = which(pruned_Zuntini_tree$tip.label == "Costus01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Erythranthe02", where = which(pruned_Zuntini_tree$tip.label == "Erythranthe01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Erythranthe03", where = which(pruned_Zuntini_tree$tip.label == "Erythranthe01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Erythranthe04", where = which(pruned_Zuntini_tree$tip.label == "Erythranthe01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Ipomoea02", where = which(pruned_Zuntini_tree$tip.label == "Ipomoea01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Jaltomata02", where = which(pruned_Zuntini_tree$tip.label == "Jaltomata01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Jaltomata03", where = which(pruned_Zuntini_tree$tip.label == "Jaltomata01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Kalmia02", where = which(pruned_Zuntini_tree$tip.label == "Kalmia01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Leptosiphon02", where = which(pruned_Zuntini_tree$tip.label == "Leptosiphon01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Leucaena02", where = which(pruned_Zuntini_tree$tip.label == "Leucaena01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Leucaena03", where = which(pruned_Zuntini_tree$tip.label == "Leucaena01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Leucaena04", where = which(pruned_Zuntini_tree$tip.label == "Leucaena01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Leucaena05", where = which(pruned_Zuntini_tree$tip.label == "Leucaena01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Leucaena06", where = which(pruned_Zuntini_tree$tip.label == "Leucaena01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Nicotiana02", where = which(pruned_Zuntini_tree$tip.label == "Nicotiana01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Nicotiana03", where = which(pruned_Zuntini_tree$tip.label == "Nicotiana01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Phlox02", where = which(pruned_Zuntini_tree$tip.label == "Phlox01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Salvia02", where = which(pruned_Zuntini_tree$tip.label == "Salvia01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Schiedea02", where = which(pruned_Zuntini_tree$tip.label == "Schiedea01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Silene02", where = which(pruned_Zuntini_tree$tip.label == "Silene01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Solanum02", where = which(pruned_Zuntini_tree$tip.label == "Solanum01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Solanum03", where = which(pruned_Zuntini_tree$tip.label == "Solanum01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Solanum04", where = which(pruned_Zuntini_tree$tip.label == "Solanum01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Streptanthus02", where = which(pruned_Zuntini_tree$tip.label == "Streptanthus01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Streptanthus03", where = which(pruned_Zuntini_tree$tip.label == "Streptanthus01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Streptanthus04", where = which(pruned_Zuntini_tree$tip.label == "Streptanthus01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Streptanthus05", where = which(pruned_Zuntini_tree$tip.label == "Streptanthus01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Streptanthus06", where = which(pruned_Zuntini_tree$tip.label == "Streptanthus01"), position = 0, edge.length = 0)
pruned_Zuntini_tree <- bind.tip(pruned_Zuntini_tree, tip.label = "Streptanthus07", where = which(pruned_Zuntini_tree$tip.label == "Streptanthus01"), position = 0, edge.length = 0)
plot(pruned_Zuntini_tree)

# make sure the tips match the PI pairs
temp1 <- pruned_Zuntini_tree$tip.label
temp2 <- table1$PI_name %>% unique()
setdiff(temp1, temp2)
setdiff(temp2, temp1)

# write tre file
write.tree(pruned_Zuntini_tree, file = "All_PI_pairs_Zuntini_2025_10.tre")



###############################
### --- OLD CODE SCRAPS --- ###
###############################

### --- temp code to lift branch lengths from one tree to another

get_genus <- function(x) sub("^([A-Za-z]+).*", "\\1", x)

# Node ages in reference
bt_ref <- branching.times(Zuntini_tree)
gen_ref <- get_genus(Zuntini_tree$tip.label)
grp_ref <- split(Zuntini_tree$tip.label, gen_ref)

crown_age_ref <- sapply(grp_ref, function(tips) {
  if (length(tips) < 2) return(NA_real_)
  nd <- getMRCA(Zuntini_tree, tips)
  if (is.null(nd)) return(NA_real_)
  bt_ref[as.character(nd)]
})

# Build calibration table for target
gen_tgt <- get_genus(phy$tip.label)
cal <- data.frame(node = integer(), age.min = numeric(), age.max = numeric())

for (g in names(crown_age_ref)) {
  age <- crown_age_ref[[g]]
  if (is.na(age)) next
  tips_g <- phy$tip.label[gen_tgt == g]
  if (length(tips_g) < 2) next
  nd_tgt <- getMRCA(phy, tips_g)
  if (is.null(nd_tgt)) next
  cal <- rbind(cal, data.frame(node = nd_tgt, age.min = age, age.max = age))
}

# Rescale target tree with genus calibrations
phy_time <- chronos(phy, lambda = 1, calibrations = cal)


### --- Cleaning code

# use to identify duplicates and possible errors in the datasheet
df <- read.csv("Duplication_list.csv")
df_new <- df %>%
  group_by(species) %>%
  mutate(duplicate_flag = if_else(!is.na(Number) & n() > 1 & sum(!is.na(Number)) > 1, "duplicate", "unique")) %>%
  mutate(weird_values = if_else(max(!is.na(Number)) > min(!is.na(Number))*1.25, "check", "fine")) %>%
  ungroup()
table(df_new$duplicate_flag)
table(df_new$weird_values)
write.csv(df_new, file="temp.csv")








library(tidyverse)
library(ape)
library(nlme)
library(phytools)

table1 <- read.csv("Table1_2025_09_18.csv")
table2 <- read.csv("Table2_2025_09_18.csv")
table2b <- read.csv("Table2b_2025_11_25.csv")
gen_dist <- read.csv("Genetic_distances_results.csv")
gen_dist <- distinct(gen_dist)
#phy <- read.tree("All_PI_pairs_v3_2025_08.tre.txt")
phy <- read.tree("All_PI_pairs_Zuntini_2025_10.tre") # see below for code to produce this


#remove duplicate/problematic values and find mean & median ovule number per genus
  # think about how to deal with genus-level values
genus_ovule_numbers <- table2b %>% 
  filter(Remove!="yes") %>% 
  group_by(species) %>% 
  summarise(genus = first(genus), level = first(level), mean = mean(number, na.rm=T), count = n()) %>%
  group_by(genus) %>%
  summarise(genus_mean_ovule_number = mean(mean, na.rm=T), genus_median_ovule_number = median(mean, na.rm=T), count = n())

hist(log(genus_ovule_numbers$genus_mean_ovule_number))
hist(log(genus_ovule_numbers$genus_median_ovule_number))
plot(log(genus_ovule_numbers$genus_mean_ovule_number), log(genus_ovule_numbers$genus_median_ovule_number))

# add genus-level ovule number to table 2 
table2_plus <- left_join(table2, genus_ovule_numbers, by = "genus")

plot(table(table2_plus$life_history, table2_plus$mating_system))
plot(log(table2_plus$ovule_number), log(table2_plus$genus_mean_ovule_number), col=table2_plus$confidence)
plot(log(table2_plus$ovule_number), log(table2_plus$genus_median_ovule_number), col=table2_plus$confidence)

# use genus-level ovule data to replace other values 
table2_final <- table2_plus %>%
  # can edit below line to use mean or replace more or less ovule numbers
  mutate(ovule_number_plus = if_else(confidence > 3, genus_median_ovule_number, ovule_number)) %>%
  select(species, ovule_number_plus, life_history, mating_system)
  
# add table 2 and genetic distance data to table 1 
table1_plus <- left_join(table1, table2_final, by = c("species1" = "species")) %>% 
  left_join(., table2_final, by = c("species2" = "species")) %>%
  left_join(., gen_dist, by = c("species1" = "Species_1", "species2" = "Species_2"))
  
# summarize table 1 data
table1_plus2 <- table1_plus %>%  mutate(mean_ov = (ovule_number_plus.x+ovule_number_plus.y)/2,
         min_ov = pmin(ovule_number_plus.x,ovule_number_plus.y),
         mating_system = case_when(mating_system.x=="outcrossing" & mating_system.y=="outcrossing" ~ "outcrossing",
                                   mating_system.x=="selfing" & mating_system.y=="selfing" ~ "selfing",
                                   .default = "mixed-mating"),
         life_history = case_when(life_history.x=="woody" & life_history.y=="woody" ~ "slow",
                                  life_history.x=="annual" & life_history.y=="annual" ~ "fast",
                                  .default = "intermediate"))

hist(table1_plus2$K2P_Distance, breaks=25)
stripchart(table1_plus2$K2P_Distance ~ table1_plus2$mating_system, vertical=T, method="jitter")
stripchart(table1_plus2$K2P_Distance ~ table1_plus2$life_history, vertical=T, method="jitter")

# filter table 1
table1_filt <- table1_plus2 %>% 
  filter(remove!="yes") %>%
  filter(!is.na(PI_name)) %>% #133
  filter(!is.na(mean_ov)) %>% #125
  filter(!is.na(K2P_Distance)) %>% #104
  #filter(K2P_Distance<0.2) %>%
  select(PI_name, RI, measure, min_ov, mean_ov, mating_system, life_history, K2P_Distance)
names(table1_filt)[1] <- "species"
table1_filt$measure <- as.factor(table1_filt$measure)
table1_filt$mating_system <- as.factor(table1_filt$mating_system)
table1_filt$life_history <- as.factor(table1_filt$life_history)

table1_pollen_comp <- table1_filt %>% filter(measure=="pollen competition")

cols <- colorRampPalette(c("blue", "red"))(100)
plot(table1_filt$K2P_Distance, table1_filt$RI, col=cols[cut(log(table1_filt$mean_ov), breaks = 100)])
plot(log(table1_filt$mean_ov), table1_filt$RI, col=cols[cut(table1_filt$K2P_Distance, breaks = 100)])


# drop tips from tre file that aren't in table1_filt
phy  <- drop.tip(phy, setdiff(phy$tip.label, table1_filt$species))
plot(phy)

# remove lines in table1_filt that aren't in tre file
table1_filt <- subset(table1_filt, species %in% phy$tip.label)
rownames(table1_filt) <- table1_filt$species

# run PGLS under BM

# models with only one predictor
fit_bm_1_gd <- gls(
  RI ~ K2P_Distance, data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_1_gd)

fit_bm_1_ov <- gls(
  RI ~ log(mean_ov), data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_1_ov)

fit_bm_1_ms <- gls(
  RI ~ mating_system, data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_1_ms)

fit_bm_1_lh <- gls(
  RI ~ life_history, data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_1_lh)


# more complicated models
fit_bm_2 <- gls(
  RI ~ K2P_Distance + log(mean_ov), data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_2)

fit_bm_int <- gls(
  RI ~ K2P_Distance * log(mean_ov), data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_int)

anova(fit_bm_null, fit_bm_1_gd, fit_bm_2, fit_bm_3)

fit_bm_3 <- gls(
  RI ~ K2P_Distance + mating_system + life_history, data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_3)

fit_bm_4 <- gls(
  RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history, data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_4)

fit_bm_5 <- gls(
  RI ~ log(mean_ov) + K2P_Distance + mating_system + life_history + measure, data = table1_filt,
  correlation = corBrownian(form = ~ 1, phy = phy), method = "ML")
summary(fit_bm_5)

anova(fit_bm_3, fit_bm_4, fit_bm_5)



# run other PGLS options - can't have branch lengths of 0
library(phylolm)
phy <- read.tree("All_PI_pairs_Zuntini_2025_10.tre") # see below for code to produce this


# code to manipulate branch lengths
hist(phy$edge.length)
phy$edge.length[phy$edge.length==0] <- 0.00001
plot(phy)

# model = "BM", "OUrandomRoot", "lambda", "delta", "kappa", etc.
fit <- phylolm(RI ~ log(mean_ov) + K2P_Distance + mating_system, data = table1_filt, phy = phy, model = "lambda")
summary(fit)

fit2 <- phylolm(RI ~ K2P_Distance + mating_system, data = table1_filt, phy = phy, model = "lambda")
summary(fit)
plot(phy)


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








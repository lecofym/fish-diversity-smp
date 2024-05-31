# Download and load packages ####
# install.packages("pacman")
devtools::install_github("kkeenan02/diveRsity")
devtools::install_github("ahasverus/elbow")
devtools::install_github("cmartin/ggConvexHull")
pacman:: p_load(tidyverse, lme4, jtools, vegan, arm, ade4, diveRsity,
                cluster, clue, FD, ggrepel, reshape2, rcompanion,
                psych, car, mgcv, stats, emstreeR, ggConvexHull, Rcpp,
                geometry, rcdd, patchwork, plotly, devtools, ape,
                gtools, elbow)

# Clean console ####
rm(list = ls())
shell('cls')

#### Calculation of ecological indicators####
# Read the density dataset
D_dataset<-read.csv("Data/PNH_Density.csv", header = TRUE,
                    row.names = 1)
names(D_dataset)

# Select species
weight.D <- D_dataset[, c(9:107)]

# Simpson index
simpson <- matrix(nrow = nrow(weight.D), ncol = 1, byrow = T)
colnames(simpson) <- c('Simpson')
row.names(simpson) <- rownames(weight.D)
for (i in 1:nrow(weight.D)) {
  simpson[i, 1] <- diversity(weight.D[i, ], index = 'simpson')
}

simpson <- cbind(D_dataset[, c(3, 5, 7:8)], simpson)

# Read the traits dataset
FM <- read.csv("Data/PNH_Traits.csv", row.names = 1)
rownames(FM) <- gsub("\\.", "_", rownames(FM))
names(FM)

# Set two traits as factors
FM$Period_of_activity <- as.factor(FM$Period_of_activity)
FM$Diet <- as.factor(FM$Diet)

# Remove transects with less than 6 FEs
FE <- matrix(nrow = nrow(weight.D), ncol=2)
a <- for (m in 1:nrow(weight.D)){
  FE[m, ]<- dim(unique(FM[rownames(FM) %in%
                            colnames(weight.D[m, apply(weight.D[m, ],
                                                       2, sum)>0]), ]))
}
nbFE <- FE[, 1]
sum(nbFE < 6) # 193 transects with less than 6 FEs
weight <- data.frame (weight.D, nbFE)
weight.fe <- weight [!(weight$nbFE < 6), ]
sum(weight.fe$nbFE < 6)
dim(weight.fe) # 752 transects
names(weight.fe)
weight1 <- weight.fe[, -100] # remove nbFE column

# Pair the trait & weight datasets
FM <- FM[rownames(FM) %in% colnames(weight1),]
check_the_names<-data.frame(rownames(FM),colnames(weight1))
sum(check_the_names$rownames.FM.==check_the_names$colnames.weight1.)

# Convert weight data file to a matrix
weight1 <- as.matrix(weight1)

# Perform a PCoA
Dist <- daisy(FM, "gower")
pco <- pcoa(Dist)
eigenvalues<-pco$values

# Evaluation of the the quality of the functional space according to Maire et al. 2015
source("R/quality_funct_space.R")
qfs <- quality_funct_space(FM, traits_weights = NULL, nbdim = 14,
                           metric = "Gower", dendro = FALSE, plot = NA)
round(qfs$meanSD, 4)

# We selected 5 dimensions since the meanSD is < 0.005 and after this value, it increases
fd.coord <- qfs$details_funct_space$mat_coord[, 1:5]

# Explained variance
Gower<-qfs$details_funct_space$mat_dissim
fit <- cmdscale(Gower, eig = TRUE, k = 4)
Var.expl<- fit$eig[fit$eig >= 0] / sum(fit$eig[fit$eig > 0])

# Evaluation of the quality of the functional space according to Mouillot et al. 2021
inflection_points <- data.frame('Dimensions' = 1:length(Var.expl),
                               'Explained.variance' = cumsum(Var.expl))
perf <- elbow(data = inflection_points[, c('Dimensions',
                                            'Explained.variance')])
# We selected 5 dimensions since its value is similar to the elbow-based optimal dimensionality for this dataset (8 dimensions)
trait.D <- pco$vectors[, 1:5]
PCscores <- data.frame(trait.D)

# Save PCoA scores
write.csv(PCscores, file="Data/PNH_PCoA_scores.csv")

# Calculate ecological indicators
source("R/multidimFD.R")
Indice.den <- multidimFD(trait.D, weight1)
Indice.den <- data.frame(Indice.den)
names(Indice.den)

# Build a dataset with the weight and indices datasets
D_dataset1 <- D_dataset[rownames(D_dataset) %in% rownames(weight1), ]
check_the_names <- data.frame(rownames(D_dataset1), rownames(weight1))
sum(check_the_names$rownames.D_dataset1. == check_the_names$rownames.weight1)
PNH <- cbind(D_dataset1, Indice.den)
names(PNH)
PNH1 <- PNH[, c(1:8, 108:135,9:107)]

#### Temporal analyses of fish diversity based in LMMs ####
# S1 Text: Table B. Temporal analyses of fish diversity at the coast of Oaxaca, based in LMMs (Index ~ year + (1|site) + (1|Season))
NP <- subset(PNH1, PNH1$Human_use_level == "Non-protected")
MPA <- subset(PNH1, PNH1$Human_use_level == "MPA")

# LMMs for NP sites
NP.lmer.S <- lmer(Nb_sp ~ Year_cons + (1|Site) + (1|Season), data= NP)
summ(NP.lmer.S)

NP.lmer.D <- lmer(log(Tot_weight, base = 2) ~ Year_cons + (1|Site) +
                    (1|Season), data= NP)
summ(NP.lmer.D)

NP.lmer.Sim <- simpson %>% filter(Human_use_level == 'Non-protected') %>%
  lmer(Simpson ~ Year_cons + (1|Site) + (1|Season), data = .)
summ(NP.lmer.Sim)

NP.lmer.FR <- lmer(FRic ~ Year_cons + (1|Site) + (1|Season), data = NP)
summ(NP.lmer.FR)

NP.lmer.FD <- lmer(FDiv ~ Year_cons + (1|Site) + (1|Season), data = NP)
summ(NP.lmer.FD)

NP.lmer.FO <- lmer(FOri ~ Year_cons + (1|Site) + (1|Season), data = NP)
summ(NP.lmer.FO)

# Dependent variables were scaled to compare the estimate of each index
NP.lmer.S.sca <- lmer(scale(Nb_sp) ~ Year_cons + (1|Site) + (1|Season), data= NP)

NP.lmer.D.sca <- lmer(scale(log(Tot_weight, base = 2)) ~ Year_cons +
                        (1|Site)+ (1|Season), data= NP)

NP.lmer.Sim.sca <- simpson %>%
  filter(Human_use_level == 'Non-protected') %>%
  lmer(scale(Simpson) ~ Year_cons + (1|Site) + (1|Season), data = .)

NP.lmer.FR.sca <- lmer(scale(FRic) ~ Year_cons + (1|Site) + (1|Season), data= NP)

NP.lmer.FD.sca <-lmer(scale(FDiv) ~ Year_cons + (1|Site) + (1|Season), data= NP)

NP.lmer.FO.sca <- lmer(scale(FOri) ~ Year_cons + (1|Site) + (1|Season), data= NP)

# Figure 3. Standardized coefficients (mean ± 95% confidence interval) of LMMs for fish ecological indicators at NP (red) reefs of Oaxaca coast
NP.coefs.S <- as.data.frame(summary(NP.lmer.S.sca)$coefficients[-1, 1:2])

NP.coefs.D <- as.data.frame(summary(NP.lmer.D.sca)$coefficients[-1, 1:2])

NP.coefs.Sim <- as.data.frame(summary(NP.lmer.Sim.sca)$coefficients[-1, 1:2])

NP.coefs.FR <- as.data.frame(summary(NP.lmer.FR.sca)$coefficients[-1, 1:2])

NP.coefs.FD <- as.data.frame(summary(NP.lmer.FD.sca)$coefficients[-1, 1:2])

NP.coefs.FO <- as.data.frame(summary(NP.lmer.FO.sca)$coefficients[-1, 1:2])

a <- cbind(NP.coefs.S, NP.coefs.D, NP.coefs.Sim, NP.coefs.FR, NP.coefs.FD, NP.coefs.FO)

colnames(a) <- c("S", "d", "D", "FRic", "FDiv", "FOri")

coefs.data <- data.frame(t(a))
colnames(coefs.data) <- c("Estimate", "se")
coefs.data$Indices <- factor(rownames(coefs.data), ordered = T, levels = rev(c("S", "d", "D", "Shannon", "FRic", "FDiv", "FOri")))
str(coefs.data$Indices)

NP.fig.coef.model<- ggplot(coefs.data, aes(x = Indices, y = Estimate, fill = Indices)) +
  ggtitle("b) NP") + geom_hline(yintercept = 0, lty = 2, lwd = 1,
                                colour = "black") +
  geom_errorbar(aes(ymin = Estimate - 1.96*se, ymax = Estimate + 1.96*se),
                lwd = 1, colour = "black", width = 0) +
  geom_point(size = 13, pch = 21, stroke = 1) +
  scale_fill_manual(values = rev(c("coral2", "coral2", "gray80",
                                   "gray80", "gray80",
                                   "gray80", "gray80"))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 15, angle = 90),
        axis.title.x = element_text(size = 30, angle = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        text = element_text(family = "TT Arial"),
        plot.title = element_text(size = 30, face = "bold")) +
  coord_flip()
NP.fig.coef.model

# LMMs for MPA sites
MPA.lmer.S <- lmer(Nb_sp ~ Year_cons + (1|Site) + (1|Season), data= MPA)
summ(MPA.lmer.S)

MPA.lmer.D <- lmer(log(Tot_weight, base = 2) ~ Year_cons + (1|Site) +
                     (1|Season), data= MPA)
summ(MPA.lmer.D)

MPA.lmer.Sim <- simpson %>% filter(Human_use_level == 'MPA') %>%
  lmer(Simpson ~ Year_cons + (1|Site) + (1|Season), data = .)
summ(MPA.lmer.Sim)

MPA.lmer.FR <- lmer(FRic ~ Year_cons + (1|Site) + (1|Season), data = MPA)
summ(MPA.lmer.FR)

MPA.lmer.FD <-lmer(FDiv ~ Year_cons + (1|Site) + (1|Season), data = MPA)
summ(MPA.lmer.FD)

MPA.lmer.FO <- lmer(FOri ~ Year_cons + (1|Site) + (1|Season), data = MPA)
summ(MPA.lmer.FO)

# Dependent variables were scaled to compare the estimate of each index
MPA.lmer.S.sca <- lmer(scale(Nb_sp) ~ Year_cons + (1|Site) + (1|Season), data = MPA)

MPA.lmer.D.sca <- lmer(scale(log(Tot_weight, base = 2)) ~ Year_cons +
                         (1|Site) + (1|Season), data = MPA)

MPA.lmer.Sim.sca <- simpson %>% filter(Human_use_level == 'MPA') %>%
  lmer(scale(Simpson) ~ Year_cons + (1|Site) + (1|Season), data = .)

MPA.lmer.FR.sca <- lmer(scale(FRic) ~ Year_cons + (1|Site) + (1|Season), data = MPA)

MPA.lmer.FD.sca <-lmer(scale(FDiv) ~ Year_cons + (1|Site) + (1|Season), data = MPA)

MPA.lmer.FO.sca <- lmer(scale(FOri) ~ Year_cons + (1|Site) + (1|Season), data = MPA)

#Figure 3. Standardized coefficients (mean ± 95% confidence interval) of LMMs for fish ecological indicators at MPA (blue) reefs of Oaxaca coast
MPA.coefs.S <- as.data.frame(summary(MPA.lmer.S.sca)$coefficients[-1, 1:2])

MPA.coefs.D <- as.data.frame(summary(MPA.lmer.D.sca)$coefficients[-1, 1:2])

MPA.coefs.Sim <- as.data.frame(summary(MPA.lmer.Sim.sca)$coefficients[-1, 1:2])

MPA.coefs.FR <- as.data.frame(summary(MPA.lmer.FR.sca)$coefficients[-1, 1:2])

MPA.coefs.FD <- as.data.frame(summary(MPA.lmer.FD.sca)$coefficients[-1, 1:2])

MPA.coefs.FO <- as.data.frame(summary(MPA.lmer.FO.sca)$coefficients[-1, 1:2])

b <- cbind(MPA.coefs.S, MPA.coefs.D, MPA.coefs.Sim, MPA.coefs.FR, MPA.coefs.FD, MPA.coefs.FO)

colnames(b) <- c("S", "d", "D", "FRic", "FDiv", "FOri")

MPA.coefs.data <- data.frame(t(b))
colnames(MPA.coefs.data) <- c("Estimate", "se")
MPA.coefs.data$Indices <- factor(rownames(MPA.coefs.data), ordered = T, levels = rev(c("S", "d", "D", "FRic", "FDiv", "FOri")))

MPA.fig.coef.model<- ggplot(MPA.coefs.data, aes(x = Indices, y = Estimate, fill = Indices)) +
  ggtitle("a) MPA") + geom_hline(yintercept = 0, lty = 2, lwd = 1,
                                 colour = "black") +
  geom_errorbar(aes(ymin = Estimate - 1.96*se, ymax = Estimate + 1.96*se),
                lwd = 1, colour = "black", width = 0) +
  geom_point(size = 13, pch = 21, stroke = 1) +
  scale_fill_manual(values = rev(c("dodgerblue3", "dodgerblue3",
                                   "dodgerblue3", "dodgerblue3",
                                   "dodgerblue3", "dodgerblue3"))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 15, angle = 90),
        axis.title.x = element_text(size = 30, angle = 0),
        axis.text.y = element_text(size = 20, angle = 0),
        axis.title.y = element_blank(),
        legend.position = "none", text = element_text(family = "TT Arial"),
        plot.title = element_text(size = 30, face = "bold")) +
  scale_y_continuous (breaks = seq(-0.1, 0.2, 0.05)) +
  coord_flip()
MPA.fig.coef.model

Estimate.graph <- MPA.fig.coef.model + NP.fig.coef.model + plot_layout(nrow = 1, guides = "collect")

ggsave('Figs/Figure_3.tiff', plot = Estimate.graph,
       width = 2250, height = 2625, units = 'px', dpi = 320,
       compression = "lzw")


#### Spatial analyses of fish diversity ####
# As trait matrix we use PCoA scores calculated in the previous section
coord <- read.csv("Data/PNH_PCoA_scores.csv")
names(coord)
coord <- data.frame(coord[c(1:6)])
colnames(coord) <- c("name", "E1", "E2", "E3", "E4", "E5")

# As weight matrix we use the the average abundance per species at each human use level
assemblage.use <- read.csv("Data/PNH_spatial_analysis.csv")
coord_com <- assemblage.use %>%
  pivot_longer(cols =! Human_use_level) %>%
  full_join(coord, by = "name") %>% filter(value > 0) %>%
  mutate_at("Human_use_level", factor,
            levels = c("MPA", "Non-protected")) %>%
  group_by(Human_use_level) %>%
  mutate(centroid.E1 = mean(E1), centroid.E2 = mean(E2),
         centroid.E3 = mean(E3), centroid.E4 = mean(E4),
         centroid.E5 = mean(E5)) %>%
  group_by(Human_use_level,name) %>%
  mutate(dist.cen = sqrt(sum((c(E1, E2, E3, E4, E5) - c(centroid.E1, centroid.E2, centroid.E3, centroid.E4, centroid.E5))^2))) %>%
  ungroup()

# Abundance weighted centroid
trait_centroid<- coord_com %>%
  group_by(Human_use_level) %>% mutate(value = value/sum(value)) %>%
  mutate(E1 = E1 * value) %>% mutate(E2 = E2 * value) %>%
  mutate(E3 = E3 * value) %>% mutate(E4 = E4 * value) %>%
  mutate(E5 = E5 * value) %>% summarise(across(E1:E5, ~mean(.x)))

# Figures format
theme_set(theme_bw())
theme_update(legend.text = element_text(size = 28),
             legend.position = "bottom",  legend.box ="horizontal",
  legend.margin = margin(), panel.grid = element_blank(),
  axis.text = element_blank(), axis.ticks=element_blank(),
  plot.title = element_text(face = "bold", size = 28),
  axis.title = element_text(size = 28))
cpalette <- c("dodgerblue3", "coral2")

# Functional richness
(F_rich1 <- coord_com %>%
    ggplot(aes(x = E1, y = E2, colour = Human_use_level, fill = Human_use_level))+
    geom_convexhull(alpha = 0.3)+
    geom_point(shape = 19, size = 1) +
    labs(x = "PCoA1", y = "PCoA2", colour = "Level",
         title = "a) FRic")+
    scale_colour_manual(values = cpalette)+
    scale_fill_manual(values = cpalette)+
    coord_fixed(ylim = c(-0.5, 0.5), xlim = c(-0.55, 0.55)) +
    theme(plot.title = element_text(size = 20),
          axis.title = element_text(size = 18),
          legend.position = "none"))

(F_rich2<- coord_com %>%
    ggplot(aes(x = E3, y = E4, colour = Human_use_level, fill = Human_use_level))+
    geom_convexhull(alpha = 0.3)+ geom_point(shape = 19, size = 1) +
    labs(x = "PCoA3", y = "PCoA4", colour = "Level")+
    scale_colour_manual(values = cpalette)+
    scale_fill_manual(values = cpalette)+
    coord_fixed(ylim = c(-0.5, 0.5), xlim = c(-0.55, 0.55)) +
    theme(axis.title = element_text(size = 18),
          legend.position = "none"))

# Functional divergence
trait_max2 <- coord_com %>% group_by(Human_use_level) %>%
  top_n(n = 1, wt = value)


(F_div1 <- coord_com %>% group_by(Human_use_level) %>%
    mutate(centroid.E1 = mean(E1),
           centroid.E2 = mean(E2),
           centroid.E3 = mean(E3),
           centroid.E4 = mean(E4),
           centroid.E5 = mean(E5)) %>%
    group_by(Human_use_level, name) %>%
    mutate(dist.cen = sqrt(sum((c(E1, E2, E3, E4, E5) - c(centroid.E1, centroid.E2, centroid.E3, centroid.E4,centroid.E5))^2))) %>%
    ggplot(aes(x = E1, y = E2))+
    geom_point(aes(size = value, colour = Human_use_level), alpha = 0.3) +
    geom_point(data = trait_centroid,
               aes(x = E1, y = E2, fill = Human_use_level),
               shape = 23, colour = cpalette, size = 6)+
    stat_ellipse(data = subset(coord_com, coord_com$Human_use_level == "MPA"),
                 aes(x = E1, y = E2), type = "euclid",
                 level = trait_max2$dist.cen[1], colour = cpalette[1])+
    stat_ellipse(data = subset(coord_com, coord_com$Human_use_level == "Non-protected"),
                 aes(x = E1, y = E2), type = "euclid",
                 level = trait_max2$dist.cen[2], colour = cpalette[2])+
    scale_colour_manual(values = cpalette)+
    scale_fill_manual(values = cpalette)+
    scale_size(range = c(1, 7))+
    labs(x = "PCoA1", y = "PCoA2", colour = "Level",
         size = "Density", title = "b) FDiv")+
    guides(colour = guide_legend(order = 1, nrow = 1),
           size = guide_legend(order = 2, nrow = 1)) +
    coord_fixed(ylim = c(-0.5, 0.5), xlim = c(-0.55, 0.55))+
    theme(plot.title = element_text(size = 20),
          axis.title = element_text(size = 18), legend.position = "none"))

(F_div2 <- coord_com %>% group_by(Human_use_level) %>%
    mutate(centroid.E1 = mean(E1),
           centroid.E2 = mean(E2),
           centroid.E3 = mean(E3),
           centroid.E4 = mean(E4),
           centroid.E5 = mean(E5)) %>%
    group_by(Human_use_level, name) %>%
    mutate(dist.cen = sqrt(sum((c(E1, E2, E3, E4, E5) - c(centroid.E1, centroid.E2, centroid.E3, centroid.E4,centroid.E5))^2))) %>%
    ggplot(aes(x = E3, y = E4))+
    geom_point(aes(size = value, colour = Human_use_level), alpha = 0.3) +
    geom_point(data = trait_centroid,
               aes(x = E3, y = E4, fill = Human_use_level),
               shape = 23, colour = cpalette,size = 6)+
    stat_ellipse(data = subset(coord_com, coord_com$Human_use_level == "MPA"),
                 aes(x = E3, y = E4), type = "euclid",
                 level = trait_max2$dist.cen[1], colour = cpalette[1])+
    stat_ellipse(data = subset(coord_com, coord_com$Human_use_level == "Non-protected"),
                 aes(x = E3, y = E4), type = "euclid",
                 level = trait_max2$dist.cen[2], colour = cpalette[2])+
    scale_colour_manual(values = cpalette)+
    scale_fill_manual(values = cpalette)+
    scale_size(range = c(1, 7))+
    labs(x = "PCoA3", y = "PCoA4", colour = "Level",
         size = "Density")+
    guides(colour = guide_legend(order = 1, nrow = 1),
           size = guide_legend(order = 2, nrow = 1)) +
    coord_fixed(ylim = c(-0.5, 0.5), xlim = c(-0.55, 0.55))+
    theme(axis.title = element_text(size = 18),
          legend.position = "none"))

# Functional originality
coord <- read.csv("Data/PNH_PCoA_scores.csv", row.names = 1)
names (coord)
trait2 <- data.frame(coord[c(1:5)])
colnames(trait2) <-c("E1", "E2", "E3", "E4", "E5")
assemblage.use <- read.csv("Data/PNH_spatial_analysis.csv", row.names = 1)

entropy_data <- sapply(c("MPA","Non-protected"), function(j){
  subtrait <- trait2[colnames(assemblage.use)[which(assemblage.use[j, ] > 0)], ]
  #subtrait means only the the trait information of species present in Site j
  subcomm <- assemblage.use[j, which(assemblage.use[j, ] > 0)]
  #subcomm means only the the abundance information of species present in Site j
  comb_dat <- t(combn(nrow(subtrait), 2))
  #combinatory matrix to extract the correct values for each possible combination
  entropy_dat <- matrix(NA, nrow = choose(nrow(subtrait), 2), ncol = 5) #create empty matrix
  colnames(entropy_dat) <- c("E1", "E1_end", "E2", "E2_end", "den")
  for (i in 1:nrow(entropy_dat)){ #fill it with for loop
    entropy_dat[i, ] <- c(subtrait[comb_dat[i, 1], 1], #coordinate x of spp 1
                         subtrait[comb_dat[i, 2], 1],#coordinate x of spp 2
                         subtrait[comb_dat[i, 1], 2], #coordinate y of spp 1
                         subtrait[comb_dat[i, 2], 2], #coordinate y of spp 2
                         as.numeric(abs(subcomm[comb_dat[i, 1]] - subcomm[comb_dat[i, 2]])))
    #difference in abundance (absolute)
  }
  return(entropy_dat)
}, simplify = FALSE, USE.NAMES = TRUE) #return a list object with all pairwise for each Site

pair_data_ent <- map_df(entropy_data, ~ data.frame(.x), .id = "Human_use_level") #combine lists

entropy_data2 <- sapply(c("MPA","Non-protected"), function(j){
  subtrait <- trait2[colnames(assemblage.use)[which(assemblage.use[j, ] > 0)], ]
  #subtrait means only the the trait information of species present in Site j
  subcomm <- assemblage.use[j, which(assemblage.use[j, ] > 0)]
  #subcomm means only the the abundance information of species present in Site j
  comb_dat <- t(combn(nrow(subtrait), 2))
  #combinatory matrix to extract the correct values for each possible combination
  entropy_dat2 <- matrix(NA, nrow = choose(nrow(subtrait), 2), ncol = 5) #create empty matrix
  colnames(entropy_dat2) <- c("E3", "E3_end", "E4", "E4_end", "biom")
  for (i in 1:nrow(entropy_dat2))
    { #fill it with for loop
    entropy_dat2[i, ] <- c(subtrait[comb_dat[i, 1], 3], #coordinate x of spp 1
                          subtrait[comb_dat[i, 2], 3],#coordinate x of spp 2
                          subtrait[comb_dat[i, 1], 4], #coordinate y of spp 1
                          subtrait[comb_dat[i, 2], 4], #coordinate y of spp 2
                          as.numeric(abs(subcomm[comb_dat[i, 1]] - subcomm[comb_dat[i, 2]])))
    #difference in abundance (absolute)
    #####
  }
  return(entropy_dat2)
}, simplify = FALSE, USE.NAMES = TRUE) #return a list object with all pairwise for each Site

pair_data_ent <- map_df(entropy_data, ~ data.frame(.x), .id = "Human_use_level") #combine lists
ori_data <- sapply(c("MPA","Non-protected"), function(j){
  subtrait <- trait2[colnames(assemblage.use)[which(assemblage.use[j, ] > 0)], ]
  subcomm <- assemblage.use[j, which(assemblage.use[j, ] > 0)]
  comb_dat <- expand.grid(rownames(subtrait), rownames(subtrait))
  entropy_dat<-matrix(NA, nrow(comb_dat), ncol = 7)
  colnames(entropy_dat) <- c("Focal_spp_1", "Focal_spp_2", "E1", "E1_end", "E2", "E2_end", "dist")
  for (i in 1:nrow(comb_dat)) {
    entropy_dat[i, ] <- c(
      rownames(subtrait)[comb_dat[i, 1]],
      rownames(subtrait)[comb_dat[i, 2]],#I need to add the name of the species as an index for later
      subtrait[comb_dat[i, 1], 1],
      subtrait[comb_dat[i, 2], 1],
      subtrait[comb_dat[i, 1], 2],
      subtrait[comb_dat[i, 2], 2],
      as.numeric(dist(rbind(
        subtrait[comb_dat[i, 1], 1:2], subtrait[comb_dat[i, 2], 1:2]
      )))
    )
  }
  return(entropy_dat)
}, simplify = FALSE, USE.NAMES = TRUE)

ori_data<- lapply(ori_data, function(x){data.frame(x) %>%
    mutate_at(vars(-starts_with("Focal_spp")), as.numeric) %>%
    mutate_at(vars(starts_with("Focal_spp")), as.factor)})

pair_data_ori <- map_df(ori_data, ~ data.frame(.x), .id = "Human_use_level") %>%
  group_by(Human_use_level, Focal_spp_1) %>%  #we group by Site and species so we have only one pair of species per Site
  slice_min(dist, n = 2) %>%
  filter(dist > 0)

# Plot the originality
#```{r}
(F_ori1 <- coord_com %>% mutate(overall.cen.E1 = mean(E1),
                                overall.cen.E2 = mean(E2)) %>%
    ggplot(aes(x = E1, y = E2, colour = Human_use_level))+
    geom_point(aes(size = value), alpha = 0.7)+
    # ggrepel::geom_text_repel(aes(label=name))+
    geom_segment(data = pair_data_ori,
                 aes(xend = E1, x = E1_end, yend = E2, y = E2_end,
                     colour = Human_use_level), linewidth = 0.2)+
    scale_colour_manual(values = cpalette)+
    scale_fill_manual(values = cpalette)+
    scale_size(range = c(0, 7))+
    labs(x = "PCoA1", y = "PCoA2", colour = "Level", size = "Density",
         title = "c) FOri")+
    guides(colour = guide_legend(order = 1, nrow = 1),
           size = guide_legend(order = 2, nrow = 1))+
    coord_fixed(ylim = c(-0.5, 0.5), xlim = c(-0.55, 0.55))+
    theme(plot.title = element_text(size = 20),
          axis.title = element_text(size = 18),
          legend.position = "none"))

ori_data2 <- sapply(c("MPA", "Non-protected"), function(j){
  subtrait <- trait2[colnames(assemblage.use)[which(assemblage.use[j, ] > 0)], ]
  subcomm <- assemblage.use[j, which(assemblage.use[j, ] > 0)]
  comb_dat <- expand.grid(rownames(subtrait), rownames(subtrait))
  entropy_dat <- matrix(NA, nrow(comb_dat), ncol = 7)
  colnames(entropy_dat) <- c("Focal_spp_1", "Focal_spp_2", "E3", "E3_end", "E4", "E4_end", "dist")
  for (i in 1:nrow(comb_dat)) {
    entropy_dat[i, ] <- c(
      rownames(subtrait)[comb_dat[i, 1]],
      rownames(subtrait)[comb_dat[i, 2]],
      subtrait[comb_dat[i, 1], 1],
      subtrait[comb_dat[i, 2], 1],
      subtrait[comb_dat[i, 1], 2],
      subtrait[comb_dat[i, 2], 2],
      as.numeric(dist(rbind(
        subtrait[comb_dat[i, 1], 1:2], subtrait[comb_dat[i, 2], 1:2]
      )))
    )
  }
  return(entropy_dat)
}, simplify = FALSE, USE.NAMES = TRUE)

ori_data2 <- lapply(ori_data2, function(x){data.frame(x) %>%
    mutate_at(vars(-starts_with("Focal_spp")), as.numeric) %>%
    mutate_at(vars(starts_with("Focal_spp")), as.factor)})

pair_data_ori2 <- map_df(ori_data2, ~ data.frame(.x), .id = "Human_use_level") %>%
  group_by(Human_use_level, Focal_spp_1) %>%  #we group by Site and species so we have only one pair of species per Site
  slice_min(dist, n = 2) %>%
  filter(dist > 0)


#Plot the originality
#```{r}
(F_ori2 <- coord_com %>% mutate(overall.cen.E3 = mean(E3),
                                overall.cen.E4 = mean(E4)) %>%
    ggplot(aes(x = E3, y = E4, colour = Human_use_level))+
    geom_point(aes(size = value), alpha = 0.7)+
    # ggrepel::geom_text_repel(aes(label=name))+
    geom_segment(data = pair_data_ori2,
                 aes(xend = E3, x = E3_end, yend = E4, y = E4_end,
                     colour = Human_use_level), size = 0.2)+
    scale_colour_manual(values = cpalette)+
    scale_fill_manual(values = cpalette)+
    scale_size(range = c(0, 7))+
    labs(x = "PCoA3", y ="PCoA4", colour = "Level", size = "Density")+
    guides(colour = guide_legend(order = 1, nrow = 1),
           size = guide_legend(order = 2, nrow = 1))+
    coord_fixed(ylim = c(-0.5, 0.5), xlim = c(-0.55, 0.55))+
    theme(axis.title = element_text(size = 18),
          legend.position = "none"))

Spatial.graph <- F_rich1 + F_div1 + F_ori1 + F_rich2 + F_div2 + F_ori2 + plot_layout(nrow = 2, guides = "collect")

ggsave('Figs/Figure_2.tiff', plot = Spatial.graph,
       width = 2250, height = 2000, units = 'px', dpi = 320,
       compression = "lzw")

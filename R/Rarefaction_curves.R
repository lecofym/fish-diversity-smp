# Load packages ####
pacman:: p_load(tidyverse, vegan, extrafont)

# Clean console ####
rm(list = ls())
shell('cls')

#### Estimation of rarefaction curves ####
# Read the density dataset
D.dataset<-read.csv("Data/PNH_Density.csv", header = TRUE,
                    row.names = 1, stringsAsFactors = T)
levels(D.dataset$Human_use_level)

# Select data from MUMPAs and NP areas
MUMPA.data <- D.dataset %>% filter(Human_use_level == 'MPA') %>%
  .[, c(9:107)]

NP.data <- D.dataset %>% filter(Human_use_level == 'Non-protected') %>%
  .[, c(9:107)]


# Rarefaction curve estimation ####
MUMPA.curve <- decostand(MUMPA.data, method = 'pa') %>%
  specaccum(., method = 'rarefaction')
MUMPA.curve <- data.frame(Human.use = 'MPA',
                          Transect = MUMPA.curve$sites,
                          Richness = MUMPA.curve$richness,
                          SD = MUMPA.curve$sd)

NP.curve <- decostand(NP.data, method = 'pa') %>%
  specaccum(., method = 'rarefaction')
NP.curve <- data.frame(Human.use = 'NP',
                          Transect = NP.curve$sites,
                          Richness = NP.curve$richness,
                          SD = NP.curve$sd)

# Plotting and saving curves ###
plot.data <- rbind(MUMPA.curve, NP.curve)

(rare.curves <- ggplot(plot.data, aes(x = Transect, y = Richness,
                              color = Human.use, fill = Human.use)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(x = Transect, ymin = (Richness - 2*SD),
                  ymax = (Richness + 2*SD)), alpha = 0.2) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 10))+
  scale_fill_manual(values = c("dodgerblue3", "coral2")) +
  scale_color_manual(values = c("dodgerblue3", "coral2")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none",
        text = element_text(family = "Arial")))

ggsave("Figures and Tables/Rarefaction curves.tiff", plot = rare.curves,
       width = 2250, height = 1312, units = 'px', dpi = 320,
       compression = "lzw", bg = 'white')

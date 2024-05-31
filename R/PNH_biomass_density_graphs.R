# Download and load packages ####
devtools:: install_github("ricardo-bion/ggtech", dependencies = TRUE)
pacman:: p_load(tidyverse, ggtech, ggpubr)

# Clean console ####
rm(list = ls())
shell('cls')

# Load the biomass dataset
FBiomass <- read.csv("Data/PNH_FBiomass.csv", sep = ";", header=T)

FBiomass$Commercial = (log(FBiomass$Commercial)+1)
FBiomass$Non_commercial = (log(FBiomass$Non_commercial)+1)

# Set theme
theme_set(theme_tech(theme = "google")+
            theme(text=element_text(size=25)) +
            theme(legend.position="NA"))

# Graph: Biomass Commercial Fishes ####
(plot1 <- FBiomass %>%
    ggplot(aes(Year, Commercial, colour = Protection))+
    geom_point(aes(shape = Protection, size = 4, alpha = 0.7))+
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                level = 0.95, alpha = 0.5, aes(fill = Protection)) +
    scale_fill_manual(values = c("dodgerblue3", "coral2"))+
    scale_color_manual(values = c("dodgerblue3", "coral2")) +
    ggpubr::stat_cor(aes(color = Protection, size = 15,
                         label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                     size = 7) +
    scale_x_continuous(breaks = seq(2006, 2020, by = 2)) +
    scale_y_continuous(limits = c(3, 11)) +
    labs(x = "",
         y = expression("Biomass commercial fish (logX+1)")))

(plot2 <- plot1 + theme(text = element_text(size = 16),
                        axis.text = element_text(size = 14)))

ggsave("Figs/FBiomass_CommercialFish.tiff", plot2,
       width = 2250, height = 1312, units = 'px', dpi = 320,
       bg= "white", compression = "lzw")

# Graph: Biomass Non-Commercial ####
(plot3 <- FBiomass %>%
    ggplot(aes(Year, Non_commercial, colour = Protection))+
    geom_point(aes(shape = Protection, size = 4, alpha = 0.7))+
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
                alpha = 0.5, aes(fill = Protection)) +
    scale_fill_manual(values = c("dodgerblue3", "coral2"))+
    scale_color_manual(values = c("dodgerblue3", "coral2")) +
    ggpubr::stat_cor(aes(color = Protection, size = 15,
                         label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                     size = 7) +
    scale_x_continuous(breaks = seq(2006, 2020, by = 2)) +
    scale_y_continuous(limits = c(3, 11)) +
    labs(x = "",
         y = expression("Biomass non-commercial fish (logX+1)")))

(plot4 <- plot3 + theme(text = element_text(size = 16),
                        axis.text = element_text(size = 14)))

ggsave("Figs/FBiomass_Non-commercialFish.tiff", plot4,
       width = 2250, height = 1312, units = 'px', dpi = 320,
       bg= "white", compression = "lzw")

# ------------------------------------------------------ Density ####

FDensity <- read.csv("Data/PNH_FDensity.csv", header = T)

# Density Commercial Fishes ####
(plot10 <- FDensity %>%
    ggplot(aes(Year, Commercial, colour = Protection))+
    geom_point(aes(shape = Protection, size = 4, alpha = 0.7))+
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
                alpha = 0.5, aes(fill = Protection)) +
    scale_fill_manual(values = c("dodgerblue3", "coral2"))+
    scale_color_manual(values = c("dodgerblue3", "coral2")) +
    ggpubr::stat_cor(aes(color = Protection, size = 15,
                         label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                     size = 7) +
    scale_x_continuous(breaks = seq(2006, 2020, by = 2))+
    scale_y_continuous(limits = c(0, 0.025)) +
    labs(x = "",
         y = expression("Density commercial fish")))

(plot11 <- plot10 + theme(text = element_text(size = 16),
                          axis.text = element_text(size = 14)))

ggsave("Figs/FDensity_CommercialFish.tiff", plot11,
       width = 2250, height = 1312, units = 'px', dpi = 320,
       bg= "white", compression = "lzw")

# Figure: Density Non Commercial Fishes ####
(plot20 <- FDensity %>%
    ggplot(aes(Year, Non_commercial, colour = Protection))+
    geom_point(aes(shape = Protection, size = 4, alpha = 0.7))+
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
                alpha = 0.5, aes(fill = Protection)) +
    scale_fill_manual(values = c("dodgerblue3", "coral2"))+
    scale_color_manual(values = c("dodgerblue3", "coral2")) +
    ggpubr::stat_cor(aes(color = Protection, size = 15,
                         label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                     size = 7) +
    scale_x_continuous(breaks = seq(2006, 2020, by = 2))+
    labs(x = "",
         y = expression("Density non-commercial fish")))

(plot21 <- plot20 + theme(text = element_text(size = 16),
                          axis.text = element_text(size = 14)))

ggsave("Figs/FDensity_Non-commercialFish.tiff", plot21,
       width = 2250, height = 1312, units = 'px', dpi = 320,
       bg= "white", compression = "lzw")

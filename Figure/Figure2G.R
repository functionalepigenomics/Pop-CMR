PD<-read.table("input",header=T)
library(ggplot2)
library(webr)
library(dplyr)
# library(tidyr)
PieDonut(PD, aes(Survived, Class, count=Freq),explode = 2)
# ggsave(plot = g, width = 6, height = 6, dpi = 300, filename = "pca_aml0.1.pdf")
# lexicon <- data.frame("Level1" = c(rep("Flavour", 11), rep("Appearance", 4)),
#                       "Level2" = c(rep("Misc", 6), rep("Pungent", 5), rep("Colour", 4)),
#                       "Level3" = c("Fresh", "Refreshing", "Soapy", "Minty", "Nutty", "Milky", "Peppery", "Sharp", "Horseradish", "Mustard hot", "Spicy", "Colourful"," Fresh Green", "Dark Green", "Bright Green")
# )
# 
# lexicon %>%
#   mutate(top_level = Level1) %>%
#   pivot_longer(1:2) %>%
#   group_by(name, value) %>%
#   mutate(width = n()) %>%
#   unique() %>%
#   arrange(name) %>%
#   group_by(name) %>%
#   mutate(ymid = as.numeric(sub("\\D+", "", name)),
#          ymax = ymid + 0.5, ymin = ymid - 0.5,
#          xmin = c(0, head(cumsum(width), -1)),
#          xmax = cumsum(width),
#          xmid = (xmax + xmin) / 2) %>%
#   ggplot(aes(xmid, ymid, fill = top_level)) +
#   geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
#                 alpha = name, color = top_level)) +
#   geomtextpath::geom_textpath(aes(y = ymid + 0.25, label = value, 
#                                   group = value)) +
#   scale_alpha_manual(values = c(1, 0.3, 0.1)) +
#   scale_fill_manual(values = c("#cd9900", "#00817e")) +
#   scale_colour_manual(values = c("#cd9900", "#00817e")) +
#   scale_y_continuous(limits = c(-0.5, 3.6)) +
#   coord_polar() +
#   theme_void() +
#   theme(legend.position = "none")
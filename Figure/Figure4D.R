# Load ggplot2
library(ggplot2)

# Create Data
data <- data.frame(
  group=LETTERS[1:2],
  value=c(0.342535787321063,1-0.342535787321063)
)

# Basic piechart
g<-ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + # remove background, grid, numeric labels
  scale_fill_manual(values =c("#feb24c","#f03b20"))+
  theme(legend.position="none")
ggsave("sas.pdf",width = 2, height = 2)
# Load ggplot2
library(ggplot2)

a<-c(21.98,7.42,7.14,6.59,4.12)
# Create Data
count.data <- data.frame(
  group=c("Metabolism process","Biosynthesis process","Binding","Cellular component","Cell organization and biogenesis","Others"),
  # group=c("Metabolism process","Biosynthesis process","Cellular component","nucleobase-containing compound metabolic process","Intracellular anatomical structure","Others"),
  prop=c(a,100-sum(a))
)

count.data <- count.data %>%
  arrange(desc(group)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
# count.data$group<-factor(count.data$group, levels = rev(levels(count.data$group)))

# Basic piechart
g<-ggplot(count.data, aes(x = "", y = prop, fill = group)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "black")+
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(6)) +
  theme_void()
ggsave("pie_5a.pdf",width = 6, height = 6)


### ========
a<-c(29.57,12.61,7.83,6.96,6.09)
# Create Data
count.data <- data.frame(
  group=c("Metabolism process","Biosynthesis process","Cellular component","nucleobase-containing compound metabolic process","Intracellular anatomical structure","Others"),
  prop=c(a,100-sum(a))
)

count.data <- count.data %>%
  arrange(desc(group)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
# count.data$group<-factor(count.data$group, levels = rev(levels(count.data$group)))

# Basic piechart
g<-ggplot(count.data, aes(x = "", y = prop, fill = group)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "black")+
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(6)) +
  theme_void()
ggsave("pie_5b.pdf",width = 6, height = 6)

### ========
a<-c(30.81,13.64,9.09,7.07,6.06)
# Create Data
count.data <- data.frame(
  group=c("Metabolism process","Biosynthesis process",
          "nucleobase-containing compound metabolic process","Cellular component","Intracellular anatomical structure","Others"),
  prop=c(a,100-sum(a))
)

count.data <- count.data %>%
  arrange(desc(group)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
# count.data$group<-factor(count.data$group, levels = rev(levels(count.data$group)))

# Basic piechart
g<-ggplot(count.data, aes(x = "", y = prop, fill = group)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "black")+
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(6)) +
  theme_void()
ggsave("pie_6.pdf",width = 6, height = 6)
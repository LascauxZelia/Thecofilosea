################################################################################
## Analyzes on genome size
################################################################################
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)
library(gridExtra)

## Reading genome sizes and number of genes
gs <- read_tsv("legioTableFormatted.tab")
ss <- read_tsv("summaryStats.tab")
ss$ShortName <- gsub(" ", "_", ss$Node)
## Harmonizing names
gs$ShortName[!(gs$ShortName %in% ss$ShortName)] ## OK

## Joining tables
df <- gs %>% right_join(ss, by = "ShortName")
df$Group <- df$Group.y

## Reorder levels
group_order <- c("Deep_anc", "Legionella", "Legionella_anc", 
                 "Coxiella", "Coxiella_anc", 
                 "Aquicella", "Aquicella_anc", 
                 "Berkiella", "Berkiella_anc",
                 "Piscirickettsia", 
                 "Gamma", "Gamma_anc",
                 "Fangia/Francisella", "Fangia_Francisella_anc", 
                 "Outgroup", "Outgroup_anc")
df$Group <- factor(df$Group, levels = group_order)

dfl <- df %>% filter(Group %in% c("Aquicella", "Berkiella", "Coxiella", 
                                  "Legionella"))

## Linear regression
lm1 <- lm(Length ~ Present - 1, data = dfl)
summary(lm1)

## Checking genome size per group
xlab <- c(1:5)
gg_size_vs_nprot <- ggplot(dfl, aes(Present, Length, colour = Group)) +
  geom_point() + 
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, col = "black", 
              formula = y ~ x - 1, size = 0.5) +
  #facet_wrap(vars(Group)) + 
  scale_y_continuous(labels = xlab, breaks = xlab*1e6) +
  labs(x = "Number of protein families",
       y = "Genome size (Mb)", 
       title = "A") +
  annotate("text", x = c(500, 500), y = c(4, 4.5)*1e6, 
           label = c(paste("italic(R) ^ 2 ==", 
                           signif(summary(lm1)$adj.r.squared, 4)), 
                     paste("slope == ", signif(lm1$coefficients[1], 4))),
           parse = TRUE, hjust = 0) +
  theme_minimal()
gg_size_vs_nprot

## Linear regression
lm1 <- lm(Length ~ Present - 1, data = dfl)
summary(lm1)

## Predict ancestral genome sizes from protein family number 
df$Length[114:224] <- predict.lm(lm1, ss[114:224,"Present"])

## Violin plot of sizes
## Setting "" instead of NA in labels
df$Label[is.na(df$Label)] <- ""
df
## Label only Deep_anc
pos <- position_jitter(width = 0.3, seed = 12345)
cex <- ifelse(df$Label[df$Group %in% group_order[1:9]] == "", 0.5, 2)
pch <- ifelse(df$Label[df$Group %in% group_order[1:9]] == "", 1, 16)

gg_ancSizes <- ggplot(subset(df, Group %in% group_order[1:9]), 
                      aes(Group, Length, colour = Group, label = Label)) +
  geom_jitter(position = pos, show.legend = FALSE, cex = cex, pch = pch) + 
  #geom_jitter(position = pos, data=subset(df, Label != ""), cex = 2, pch = 16) +
  geom_boxplot(data=subset(df, Group %in% group_order[2:9]), alpha = 0.5, 
               outlier.shape = NA, show.legend = FALSE) +
  scale_x_discrete(limits = rev) + 
  scale_y_continuous(labels = xlab, breaks = xlab*1e6) +
  coord_flip() +
  geom_text_repel(position = pos, force = 5, max.iter = 100000, 
                  min.segment.length = 0, show.legend = FALSE, cex=3) +
  labs(x = NULL, y = "Genome Size (Mb)", colour = NULL) +
  theme_minimal() +
  geom_hline(yintercept = 1650000, linetype = "dotted", color = "black") +
  annotate("text", y = 1600000, x = 1.6, label = "Ca. Pokemonas kadabra (1.65 Mbp)", size=2.9, angle=90, fontface = "italic") +
  geom_hline(yintercept = 2370000, linetype = "dotted", color = "black") +
  annotate("text", y = 2320000, x = 1.6, label = "Ca. Fiscibacter Pecunius (2.37 Mbp)", size=2.9, angle=90, fontface = "italic")
gg_ancSizes

gg_ancetralSizesProt <- grid.arrange(gg_size_vs_nprot, gg_ancSizes, nrow = 1)
ggsave("AncestralGenomeSizes.pdf", plot = gg_ancSizes, 
       device = "pdf", width = 28, height = 18.5, unit = "cm")

## Calculate stats
df %>%
  group_by(Group) %>%
  select(Present, Length) %>%
  summarise(
    Present = mean(Present, na.rm = TRUE),
    Length = mean(Length, na.rm = TRUE)
  )

# # A tibble: 16 x 3
# Group                  Present   Length
# <fct>                    <dbl>    <dbl>
#   1 Deep_anc                 1795. 2338119.
# 2 Legionella               2665. 3418420.
# 3 Legionella_anc           2650. 3452557.
# 4 Coxiella                 1194. 1678195.
# 5 Coxiella_anc             1384. 1803492.
# 6 Aquicella                1436. 1921180 
# 7 Aquicella_anc            1610. 2097634.
# 8 Berkiella                2259. 3403412 
# 9 Berkiella_anc            2173. 2831289.
# 10 Piscirickettsia          1657  3450064 
# 11 Gamma                    2841. 4857046.
# 12 Gamma_anc                2656. 3460535.
# 13 Fangia/Francisella       1520. 2419759 
# 14 Fangia_Francisella_anc   1521. 1982085.
# 15 Outgroup                 2206. 3930047.
# 16 Outgroup_anc             2049. 2668796.

## Show genome size in some ancestors
df[df$Label != "", c("ShortName", "Length", "Present", "Label")]
# A tibble: 10 x 4
# ShortName   Length Present Label              
# <chr>        <dbl>   <dbl> <chr>              
#   1 55        2989615.   2295. Legionella (55)    
# 2 56        2345571.   1800. Legionella (56)    
# 3 57        2241649.   1721. Legionellaceae (57)
# 4 78        2048019.   1572. Coxiella (78)      
# 5 93        2111858.   1621. Aquicella (93)     
# 6 94        2090746.   1605. Coxiellaceae (94)  
# 7 95        2280987.   1751. LLCA ss (95)       
# 8 96        2194267.   1684. LLCA sl (96)       
# 9 97        2167027.   1663. LLCA Pisci (97)    
# 10 107       2710195.   2080. LFL (107)     

df %>% 
  filter(Group %in% c("Aquicella", "Coxiella")) %>%
  summarise(
    Present = mean(Present, na.rm = TRUE),
    Length = mean(Length, na.rm = TRUE)
  )

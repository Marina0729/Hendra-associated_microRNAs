.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(edgeR)
library(limma)
library(tidyverse)

microRNA_counts <- read_csv("Data/Redlands_counts.csv")

d0_d3 <- microRNA_counts %>% 
  select(gene, s1, s2, s3, s7, s8, s11, s12, s13)


redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv")
  
  
metadata_numeric <- redlands_horse_metadata %>% 
  filter(sample != "s6") %>%
  mutate(day = sub("d","", condition)) %>% 
  select(-condition) %>% 
  mutate(day = as.numeric(day)) %>%
  filter(day != 5, day !=7)
  


horse_counts <- DGEList(counts = d0_d3[, -1], genes = d0_d3[, 1], samples = metadata_numeric)

horse_counts

# filter with filterByExpr

design_matrix <- model.matrix(~ day + animal, data = metadata_numeric)

gene_filter <- filterByExpr(horse_counts, min.count = 6, min.total.count = 80, design = design_matrix)

horse_filtered <- horse_counts[gene_filter, , keep.lib.sizes = FALSE]

horse_filtered

# normalise TMM

horse_norm <- calcNormFactors(horse_filtered)


# normalise voom because library sizes are quite variable
horse_voom <- voom(horse_norm, plot = T, design = design_matrix)


plotMDS(horse_voom, col = as.numeric(as.factor(horse_voom$targets$animal)))

horse_voom$E
horse_voom$targets


# fit a linear model for each gene
horse_lm <- lmFit(horse_voom)

# calculate test statistics

horse_stats <- eBayes(horse_lm)

topTable(horse_stats, number = 20)


plot_gene <- function(gene_name) {
  horse_voom$targets %>% 
    as_tibble %>% 
    mutate(expn = horse_voom$E[horse_voom$genes$gene == gene_name, ]) %>% 
    ggplot(aes(x = day, y = expn, colour = animal)) +
    geom_point() +
    labs(title = gene_name)
}

topTable(horse_stats, coef = 2:5)

plot_gene("5_64196")
plot_gene("eca-miR-381")
plot_gene("eca-miR-127")
plot_gene("eca-miR-379")
plot_gene("eca-miR-143")
plot_gene("18_31552")
plot_gene("20_38845")
plot_gene("eca-miR-215")
plot_gene("8_76957")
plot_gene("eca-miR-146a") # negative regulator of inflammation 
plot_gene("eca-miR-223") #negative reulator of innate immunty 
plot_gene("eca-miR-10a")
plot_gene("eca-miR-21")
plot_gene("eca-miR-24")
plot_gene("eca-miR-145")
plot_gene("20_38845")

#doing some overall plotting with expression data 
microRNA_expn <- horse_voom$E %>% 
  as.tibble %>%
  bind_cols(as.tibble(horse_voom$genes)) %>% 
  gather(sample, expression, -gene) %>% 
  left_join(redlands_horse_metadata, by = "sample") 
  
ggplot(microRNA_expn, aes(x = condition, y = expression )) +
  geom_boxplot()

ggplot(microRNA_expn, aes(x = animal, y = expression )) +
  geom_boxplot()

tail(microRNA_expn)


#spread and convert to a matrix for PCA analysis 

scaled_microRNAs <- microRNA_expn %>%
  spread(gene, expression) %>%
  select(-animal, -condition) %>% 
  column_to_rownames("sample") %>% 
  scale()

scaled_microRNAs[1:5, 1:5]
plot(scaled_microRNAs)

pca_microRNAs <- prcomp(scaled_microRNAs)

summary(pca_microRNAs)
plot(pca_microRNAs)

PC1_PC2 <- pca_microRNAs$x %>% 
  as_tibble(rownames = "sample") %>%
  gather(PC, expression, -sample) %>% 
  left_join(redlands_horse_metadata, by = "sample") %>%
  spread(PC, expression) %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_text(aes(label = condition, color = animal))
  
  
PC2_PC3 <- pca_microRNAs$x %>% 
  as_tibble(rownames = "sample") %>%
  gather(PC, expression, -sample) %>% 
  left_join(redlands_horse_metadata, by = "sample") %>%
  spread(PC, expression) %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_text(aes(label = condition, color = animal))

PC3_PC4 <- pca_microRNAs$x %>% 
  as_tibble(rownames = "sample") %>%
  gather(PC, expression, -sample) %>% 
  left_join(redlands_horse_metadata, by = "sample") %>%
  spread(PC, expression) %>% 
  ggplot(aes(x = PC3, y = PC4)) +
  geom_text(aes(label = condition, color = animal))







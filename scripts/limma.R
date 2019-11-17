.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(edgeR)
library(limma)
library(tidyverse)

microRNA_counts <- read_csv("Data/Redlands_counts.csv") %>%                 #name it microRNA_counts
  select(-s6)

redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv") %>%  #Name it redlands_horse_metadata
  filter(sample != "s6")

horse_counts <- DGEList(counts = microRNA_counts[, -1], genes = microRNA_counts[, 1], samples = redlands_horse_metadata)
horse_counts

# filter with filterByExpr

design_matrix <- model.matrix(~ condition + animal, data = redlands_horse_metadata)

gene_filter <- filterByExpr(horse_counts, min.count = 1, min.total.count = 20, design = design_matrix)

horse_filtered <- horse_counts[gene_filter, , keep.lib.sizes = FALSE]

horse_filtered

# normalise TMM

horse_norm <- calcNormFactors(horse_filtered)

plot(horse_norm)

# normalise voom because library sizes are quite variable
horse_voom <- voom(horse_norm, plot = T, design = design_matrix)

horse_voom


plotMDS(horse_voom, col = as.numeric(as.factor(horse_voom$targets$animal)))

horse_voom$E



# fit a linear model for each gene
horse_lm <- lmFit(horse_voom)

# calculate test statistics

horse_stats <- eBayes(horse_lm)

topTable(horse_stats)

plot_gene <- function(gene_name) {
  horse_voom$targets %>% 
    as_tibble %>% 
    mutate(expn = horse_voom$E[horse_voom$genes$gene == gene_name, ]) %>% 
    ggplot(aes(x = condition, y = expn, colour = animal)) +
    geom_point() +
    labs(title = gene_name)
}
topTable(horse_stats, coef = 2:5)

plot_gene("eca-miR-30c")

plot_gene("18_31552")

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







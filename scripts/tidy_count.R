.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(tidyverse)
library(cowplot)

read_csv("data/Redlands_counts.csv")                                    #Read in the data 
microRNA_counts <- read_csv("Data/Redlands_counts.csv")                 #name it microRNA_counts

read_csv("data/redlands_horse_metadata.csv")                            #Read in the metadata
redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv") #Name it redlands_horse_metadata

#Install edgeR package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)




microRNA_counts

redlands_horse_metadata


#convert condition column to numeric 

redlands_horse_metadata_tidy <- redlands_horse_metadata %>% 
  mutate(day = sub("d","", condition )) %>% 
  mutate(day = as.numeric(day)) %>% 
  select(-condition)

#tidy the microRNA_count data and join to metadata 
microRNA_counts_tidy <- microRNA_counts %>%
  gather(sample, counts, -gene) %>%
  left_join(redlands_horse_metadata_tidy, by = "sample") %>%
  select(-sample) %>% 
  group_by(gene)
  summarise(mean_counts = mean(counts))
  
#Filtering to remove lowly expressed genes
#Genes with very low counts across all libraries provide little evidence for differential expression and they 
#interfere with some of the statistical approximations that are used later in the pipeline. They also add to 
#the multiple testing burden when estimating false discovery rates, reducing power to detect differentially 
#expressed genes. These genes should be filtered out prior to further analysis.
#Want to create a minimum counts threshold present in at least three replicates at each timepoint. Would ususally
#use CPM or RPKM to account for the different library sizes between samples. But have no information about library 
#size so will do filtering based on counts. 

#what is the average counts at each time point as these represent replicates?  

d0 <- microRNA_counts_tidy %>% 
  filter(day == 0) %>%
  filter(counts > 100)

d1 <- microRNA_counts_tidy %>% 
  filter(day == 1 ) %>%
  filter(counts > 100)

d3 <- microRNA_counts_tidy %>% 
  filter(day == 3 ) %>%
  filter(counts > 1000)

d5 <- microRNA_counts_tidy %>% 
  filter(day == 5 ) %>%
  filter(counts > 1000)

d7 <- microRNA_counts_tidy %>% 
  filter(day == 7) %>%
  filter(counts > 1000)

days0_1 <- full_join(d0, d1)

days3_5 <- full_join(d3, d5)

days0_5 <- full_join(days0_1, days3_5)


days0_7 <- full_join(days0_5, d7)

days0_7

expressedmicroRNA_counts <- days0_7

expressedmicroRNA_counts

#making some plots
plot_counts_day <- ggplot(data = expressedmicroRNA_counts) +
  geom_point( mapping = aes(x = day,
                            y = counts)) + 
  facet_wrap(~ gene, nrow = 8) +
  scale_y_log10() +
  
  theme(axis.title = element_text(size = 2), 
        axis.text.x = element_text(size = 2), 
        axis.text.y = element_text(size = 5))





miR103 <- expressedmicroRNA_counts %>% 
  filter(gene == "eca-miR-103")

plot_miR103_counts_day <- ggplot(data = miR103) +
  geom_line( mapping = aes(x = day,
                            y = counts, 
                            color = animal))

miR16 <- expressedmicroRNA_counts %>% 
  filter(gene == "eca-miR-16")

plot_miR16_counts_day <- ggplot(data = miR16) +
  geom_line( mapping = aes(x = day,
                           y = counts, 
                           color = animal))



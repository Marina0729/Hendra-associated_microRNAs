.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(tidyverse)
library(cowplot)

read_csv("data/Redlands_counts.csv")                                    #Read in the data 
microRNA_counts <- read_csv("Data/Redlands_counts.csv")                 #name it microRNA_counts

read_csv("data/redlands_horse_metadata.csv")                            #Read in the metadata
redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv") #Name it redlands_horse_metadata


microRNA_counts

redlands_horse_metadata

#Filtering to remove lowly expressed genes
#Genes with very low counts across all libraries provide little evidence for differential expression and they 
#interfere with some of the statistical approximations that are used later in the pipeline. They also add to 
#the multiple testing burden when estimating false discovery rates, reducing power to detect differentially 
#expressed genes. These genes should be filtered out prior to further analysis.

#Install edgeR package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")



Day0counts <- microRNA_counts["s1", "s6", "s11"]

microRNA_counts %>% 
  select()





#tidy the microRNA_count data and join to metadata 
microRNA_counts_tidy <- microRNA_counts %>%                             #assign to new variable
  gather(sample, counts, -gene) %>%                                     #gather to get sample names into a column
  left_join(redlands_horse_metadata, by = "sample") %>%                 #join by sample names
  rename(day = condition) %>%                                           #rename column 
  select(-sample)                                                       #remove sample column
  


#making some plots
plot_counts_day <- ggplot(data = redlands_counts_tidy) +
  geom_point( mapping = aes(x = day,
                            y = counts, 
                            color = animal))





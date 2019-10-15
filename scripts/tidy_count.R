.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(tidyverse)
library(cowplot)
library(modelr)


read_csv("data/Redlands_counts.csv")                                    #Read in the data 
microRNA_counts <- read_csv("Data/Redlands_counts.csv")                 #name it microRNA_counts

read_csv("data/redlands_horse_metadata.csv")                            #Read in the metadata
redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv") #Name it redlands_horse_metadata

microRNA_counts

redlands_horse_metadata

#long winded way of integrating metadata with count data to create unique genes in first column
counts <- microRNA_counts %>% 
  rename(h1d0 = s1) %>% 
  rename(h1d1 = s2) %>% 
  rename(h1d3 = s3) %>% 
  rename(h1d5 = s4) %>% 
  rename(h1d7 = s5) %>% 
  rename(h2d0 = s6) %>% 
  rename(h2d1 = s7) %>% 
  rename(h2d3 = s8) %>% 
  rename(h2d5 = s9) %>% 
  rename(h2d7 = s10) %>% 
  rename(h3d0 = s11) %>% 
  rename(h3d1 = s12) %>% 
  rename(h3d3 = s13) %>% 
  rename(h3d5 = s14) %>% 
  rename(h3d7 = s15) %>% 
  select(-h1d0, -h2d0, -h3d0) %>% 
  filter(h1d1 !=0, h1d3 !=0, h1d5 !=0, h1d7 !=0, h2d1 !=0, h2d3 !=0, h2d5 !=0, h2d7 !=0, h3d1 !=0, h3d3 !=0, h3d5 !=0, h3d7 !=0)

#how about joining the two data frames to create columns with days 
#must first convert condition column to numeric 

redlands_horse_metadata_tidy <- redlands_horse_metadata %>% 
  mutate(day = sub("d","", condition )) %>% 
  mutate(day = as.numeric(day)) %>% 
  select(-condition)

#tidy the microRNA_count data and join to metadata 
microRNA_counts_tidy <- microRNA_counts %>%
  gather(sample, counts, -gene) %>%
  left_join(redlands_horse_metadata_tidy, by = "sample") %>%
  select(-sample)
  
#Add this pipe if want to calculate mean counts per day  
  arrange(gene) %>% 
  group_by(gene, day) %>% 
  summarise(avg_counts = mean(counts), sd_counts = sd(counts)) %>% 
  mutate(cv = sd_counts / avg_counts) %>% 
  mutate(diff = avg_counts - lag(avg_counts, default = 0)) %>% 
  filter(gene = )

microRNA_counts_tidy %>% 
  spread()


#Modelling for each microRNA 
miR103 <- microRNA_counts_tidy %>% 
  filter(gene == "eca-miR-103")

miR103

ggplot(miR103, aes(day, counts)) +
  geom_point() +
  geom_smooth(model = glm)

miR103_mod <- glm(counts ~ day(), data = miR103)
coef(miR103_mod)
summary()
#(Intercept)    day 
#4958.5041    462.2175 

library(splines)
mod1 <- lm(counts ~ ns(day, 1), data = miR103)
coef(miR103_mod)


miR103 <- microRNA_counts_tidy %>% 
  filter(gene == "eca-miR-103")

miR103

ggplot(miR103, aes(day, counts)) +
  geom_point() +
  geom_smooth(model = lm)

miR103_mod <- lm(counts ~ day, data = miR103)
coef(miR103_mod)
#(Intercept)    day 
#7942.7833   -106.2167







plot_counts_day <- ggplot(data = microRNA_counts_tidy) +
  geom_point( mapping = aes(x = day,
                            y = counts)) + 
  facet_wrap(~ gene, nrow = 8) +
  scale_y_log10() +
  theme(axis.title = element_text(size = 2), 
        axis.text.x = element_text(size = 2), 
        axis.text.y = element_text(size = 5))


  
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



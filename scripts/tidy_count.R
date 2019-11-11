.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(tidyverse)
library(cowplot)
library(modelr)
library(splines)
library(broom)
library(edgeR)
library(limma)
library(RColorBrewer)
library(gplots)


read_csv("data/Redlands_counts.csv")                                    #Read in the data 
microRNA_counts <- read_csv("Data/Redlands_counts.csv")                 #name it microRNA_counts

read_csv("data/redlands_horse_metadata.csv")                            #Read in the metadata
redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv") #Name it redlands_horse_metadata

microRNA_counts

#Quick look at means and variance across libraries 
plot1 <- microRNA_counts %>% 
  gather(sample, counts, -gene) %>%
  mutate(sample = sub("s","", sample )) %>% 
  ggplot(microRNA_counts, mapping = aes(x = sample, y= counts, group = sample)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_x_discrete(limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  theme_bw() +
  labs(
    title = "Distribution of counts across libraries"
  )

redlands_horse_metadata


redlands_horse_metadata_tidy <- redlands_horse_metadata %>% 
  mutate(day = sub("d","", condition )) %>%                   #make a new column with no "d" from day column
  mutate(day = as.numeric(day)) %>%                           #convert new column to numeric
  select(-condition)                                          #remove old "condition" column


#join to metadata
microRNA_counts_tidy <- microRNA_counts %>%                   #no sequencing data came back for sample 6, removed from df
  gather(sample, counts, -gene) %>%                           #gathered all values into sample (key) and counts (value) columns leaving the gene untouched
  left_join(redlands_horse_metadata_tidy, by = "sample") %>%  #join the two tidy data frames by "sample"
  rename(library = sample) %>% 
  mutate(library = sub("s","", library))

#Plots to describe counts across libraries and samples
ggplot(microRNA_counts_tidy, aes( x = library, y= counts, color = animal, alpha = day)) +
  geom_jitter(size = 1) +
  geom_boxplot(outlier.shape = NA, aes(color = NA), alpha = 0.2, show.legend = FALSE) +
  scale_y_log10()+
  scale_x_discrete(limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  theme_bw() +
  labs(
    title = "Total miRNA counts"
  )

ggplot(microRNA_counts_tidy, aes( x = library, y= counts, fill = day, color = animal)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_classic() +
  scale_fill_viridis_c()+
  scale_x_discrete(limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  labs(
    title = "Total miRNA counts"
  )


ggplot(microRNA_counts_tidy, aes( x = gene , y = animal, fill= counts)) + 
  geom_tile() +
  facet_wrap( ~ day)
#not very informative  

#How many counts in each library?
microRNA_counts_tidy %>% 
  group_by(library) %>% 
  summarise(counts_sum = sum(counts)) %>% 
  ggplot(aes(y = counts_sum, x = library)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
  

#Does variance differ at different count ranges? 
plot1_above10000 <- microRNA_counts %>% 
    gather(sample, counts, -gene) %>% 
    mutate(sample = sub("s","", sample )) %>% 
    filter(counts > 10000) %>% 
    ggplot(microRNA_counts, mapping = aes(x = sample, y= counts, group = sample)) +
    geom_boxplot() +
    scale_y_log10() +
    scale_x_discrete(limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
    theme_bw() +
    labs(
      title = "Distribution of counts above 10,000 across libraries"
    )

plot1_above10000hist <- microRNA_counts %>% 
  gather(sample, counts, -gene) %>% 
  mutate(sample = sub("s","", sample )) %>% 
  filter(counts > 10000) %>% 
  ggplot(aes(x= counts)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample ) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Distribution of counts above 10,000 across libraries"
  )

ggplot(aes(x = sample_value)) + geom_histogram() + facet_wrap(~ sample_id) +
  labs(title = "Height")

plot1_1000to10000 <- microRNA_counts %>% 
  gather(sample, counts, -gene) %>% 
  mutate(sample = sub("s","", sample )) %>% 
  filter(counts > 1000 | counts < 10000) %>% 
  ggplot(microRNA_counts, mapping = aes(x = sample, y= counts, group = sample)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_x_discrete(limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  theme_bw() +
  labs(
    title = "Distribution of counts 1000 to 10,000 across libraries"
  )
 
plot1_100to1000 <- microRNA_counts %>% 
  gather(sample, counts, -gene) %>% 
  mutate(sample = sub("s","", sample )) %>% 
  filter(counts < 1000 & counts > 100) %>% 
  ggplot(microRNA_counts, mapping = aes(x = sample, y= counts, group = sample)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_x_discrete(limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  theme_bw() +
  labs(
    title = "Distribution of counts between 100 and 1000 across libraries"
  )

plot1_below100 <- microRNA_counts %>% 
  gather(sample, counts, -gene) %>% 
  mutate(sample = sub("s","", sample )) %>% 
  filter(counts < 100) %>% 
  ggplot(microRNA_counts, mapping = aes(x = sample, y= counts, group = sample)) +
  geom_boxplot() +
  scale_x_discrete(limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  theme_bw() +
  labs(
    title = "Distribution of counts below 100 across libraries"
  )

plot_grid(plot1_below100, plot1_100to1000, plot1_above1000)


#Following RNAseq tutorial 
#Create a matrix with only counts and gene names as column names so can 
#use cpm() function

# Remove first column from seqdata
countdata <- microRNA_counts[,-(1)]

# Store GeneID as rownames
rownames(countdata) <- microRNA_counts$gene

View(countdata)
head(countdata)

#remove s6 and calculate CPM
countdata_nos6 <- countdata %>% 
  select(-s6)

myCPM <- cpm(countdata_nos6)

#Have a look at the output
head(myCPM)
tail(myCPM)

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
#need to convert to tibble and tidy for plotting

microRNA_cpm <- as_tibble(myCPM, rownames = "gene")

tail(microRNA_cpm)


microRNA_counts_nos6 <- microRNA_counts %>% 
  select(-s6)

tail(microRNA_counts_nos6)

microRNA_cpm_tidy <- microRNA_cpm %>%
  gather(sample, cpm, -gene)                           
  

combined_cpm_counts <- bind_cols(microRNA_cpm_tidy, microRNA_counts_nos6_tidy) %>% 
  select(-gene1, -sample1) %>%
  left_join(redlands_horse_metadata_tidy, by = "sample") %>% 
  rename(library = sample) %>% 
  mutate(library = sub("s","", library))
  
plot_cpm_counts <- ggplot(combined_cpm_counts, aes( x = counts, y = cpm)) +
  geom_point(alpha = 0.2) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 3)) +
  geom_hline(yintercept = 0.5, color = "red", size = 0.5) +
  facet_wrap( ~ library)

ggsave(filename = "results/cpm_vs_counts.png", plot = plot_cpm_counts)


# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs

head(thresh)

# Summary of how many TRUEs there are in each row
table(rowSums(thresh))

# we would like to keep genes that have at least 3 TRUES in each row of thresh
keep <- rowSums(thresh) >= 3

keep
as_tibble(keep, rownames = "gene")

#filter combined cpm and counts by keep tibble 



# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata_nos6[keep,]
summary(keep)
counts.keep

#dimensions of tibble. 572 genes have at least three TRUES 
dim(counts.keep)


  
#plot data to look at distributions for each gene 
gene_summaries <- microRNA_counts_tidy %>% 
  gather(sample, counts, -gene) %>% 
  left_join(redlands_horse_metadata_tidy, by ="sample") %>% 
  select(-sample) %>%
  group_by(gene) %>%
  summarize(
    avg_count = mean(counts),
    avg_log_count = mean(log(1 + counts))
    )


ggplot(gene_summaries, aes(avg_log_count)) +
  geom_density()


#average counts across the three horses and filtering out lowly expressed microRNAs
mean_counts <- microRNA_counts_tidy %>% 
  group_by(gene, day) %>% 
  summarise(avg_count = mean(counts)) %>% 
  filter(avg_count > 10)

#plot the data
ggplot(data = mean_counts, mapping = aes(x = day, y = avg_count, group = gene)) +
  geom_line(alpha = 0.2) +
  scale_y_log10() +
  facet_wrap(~ gene)

#A lot of miRs seem to go down after day 5, perhaps we should do a regression on the day 0 to 5 only

#Too many microRNAs to plot! want to determine if the relationship between day and counts is linear for each microRNA
#then can filter out those without a linear relationship

get_lm_se <- function(microRNA_counts_tidy){            #create a variable where we get standard error from lm()
  fit <- lm(counts ~ day, data = microRNA_counts_tidy)  #fit a linear model to y=counts, x=day using microRNA_counts_tidy
  data.frame(term = names(fit$coefficients),            #structure the data frame
             slope = fit$coefficients,                  #return slope coefficient
             se = summary(fit)$coefficient[,2],         #return standard error coefficient
             rsq = summary(fit)$r.squared)              #return the r squared 
}                                                       #close nested function

summary(fit)

mod_microRNAs <- microRNA_counts_tidy %>%                             
  group_by(gene, animal) %>%                            #apply model by gene
  do(get_lm_se(.))                                      #do anything function


#remove intercept rows, remove term column and filter out genes with at least one slope = 0 and rsq more than 0.9 
mod_microRNAs_slopes <- mod_microRNAs %>% 
  filter(term == "day") %>% 
  select(-term) %>%
  group_by(gene) %>% 
  filter(!any(slope == 0)) %>% 
  filter(rsq > 0.9) %>% 
  filter(slope > 1| slope < -1)


ggplot(data = mod_microRNAs_slopes, mapping = aes(x = gene, 
                                                  y = slope,
                                                  color = animal)) + 
  geom_point() + 
  geom_linerange(aes(ymin = slope - se, ymax = slope + se))


#Want to be able to filter the microRNA_counts tidy data frame using the gene column in this mod_microRNAs_slopes data frame
filter(gene %in% c(mod_microRNAs_slopes$gene))

#plotting for some microRNAs of interest

#miR-221
miR1403p <- microRNA_counts_tidy %>% 
  filter(gene == "eca-miR-140-3p")

miR1403p

plot_miR1403p <- ggplot(miR1403p, aes(x = day, y = counts, color = animal)) +
  geom_point() +
  scale_y_log10() +
  stat_summary(aes(y = counts, group = 1), fun.y = mean, colour="red", geom ="line", group = 1) +
  labs(
  title = "miR1403p",                      
  x = "Day post infection",             
  y = "Counts")                          

ggsave(filename = "results/miR1403p.png", plot = plot_miR1403p, width = 15, height = 15, dpi = 600, units = "cm")

#horse 3 time point 7 is consitently lower in counts across all microRNAs 


#getting some stats on the modelling
miR103_mod <- glm(counts ~ day(), data = miR103)
summary(miR103_mod)

mod1 <- lm(counts ~ ns(day, 1), data = miR103)
summary(miR103_mod)


  
#Filtering to remove lowly expressed genes

d0 <- microRNA_counts_tidy %>% 
  filter(day == 0) %>%
  filter(counts > 10)

d1 <- microRNA_counts_tidy %>% 
  filter(day == 1 ) %>%
  filter(counts > 10)

d3 <- microRNA_counts_tidy %>% 
  filter(day == 3 ) %>%
  filter(counts > 10)

d5 <- microRNA_counts_tidy %>% 
  filter(day == 5 ) %>%
  filter(counts > 10)

d7 <- microRNA_counts_tidy %>% 
  filter(day == 7) %>%
  filter(counts > 10)

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
  scale_y_log10()




.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(tidyverse)
library(cowplot)
library(modelr)
library(splines)
library(broom)
library(reprex)

install.packages("reprex")

read_csv("data/Redlands_counts.csv")                                    #Read in the data 
microRNA_counts <- read_csv("Data/Redlands_counts.csv")                 #name it microRNA_counts

read_csv("data/redlands_horse_metadata.csv")                            #Read in the metadata
redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv") #Name it redlands_horse_metadata

microRNA_counts

redlands_horse_metadata


# joining the two data frames to create columns with days 
# first convert condition column to numeric 

redlands_horse_metadata_tidy <- redlands_horse_metadata %>% 
  mutate(day = sub("d","", condition )) %>%                   #make a new column with no "d" from day column
  mutate(day = as.numeric(day)) %>%                           #convert new column to numeric
  select(-condition)                                          #remove old "condition" column

#eliminate s6 and join to metadata then filter out day 7 for linear regression
microRNA_counts_tidy <- microRNA_counts %>%
  select(-s6) %>%                                             #no sequencing data came back for sample 6, removed from df
  gather(sample, counts, -gene) %>%                           #gathered all values into sample (key) and counts (value) columns leaving the gene untouched
  left_join(redlands_horse_metadata_tidy, by = "sample") %>%  #join the two tidy data frames by "sample"
  select(-sample) %>%
  filter(day != 7)

microRNA_counts_tidy %>% 
  filter(gene == %in%$mod_microRNA_slopes)
    

#average counts across the three horses and filtering out lowly expressed microRNAs
mean_counts <- microRNA_counts_tidy %>% 
  group_by(gene, day) %>% 
  summarise(avg_count = mean(counts)) %>% 
  filter(avg_count > 10)

#plot the data
ggplot(data = mean_counts, mapping = aes(x = day, y = avg_count, group = gene)) +
  geom_line(alpha = 0.2) +
  scale_y_log10()

#Al lot of miRs seem to go down after day 5, perhaps we should do a regression on the day 0 to 5 only

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


slopes <-  mod_microRNAs_slopes %>% 
  select(gene, slope, rsq)


slopes 

lm_genes <- slopes %>% 
  select(gene)

lm_genes


df1 <-  matrix(c("eca-miR-486", -3, 0.9352, "eca-miR-16", -3, 0.9436), nrow=2, ncol=3, byrow = TRUE)
df1

microRNA_counts_tidy %>% 
  filter(gene == lm_genes)

(y <- 1:4)
mean(y)
  
bigtibble <- matrix(2584, 3425, 4352, 356, 352, 453, "h1", "h2", "h3", "h1", "h2", "h3", 0, 0, 0, 1, 1, 1, nrow = 6, ncol = 4, byrow = FALSE, dimnames = list(c("mir1","mir1", "mi1","mir1",   c("gene","counts", "animal", "day")))
bigtibble                    

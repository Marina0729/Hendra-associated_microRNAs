.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(tidyverse)
library(cowplot)
library(modelr)
library(splines)
library(broom)

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

#tidy the microRNA_count data and join to metadata 
microRNA_counts_tidy <- microRNA_counts %>%
  select(-s6) %>%                                             #no sequencing data came back for sample 6, removed from df
  gather(sample, counts, -gene) %>%                           #gathered all values into sample (key) and counts (value) columns leaving the gene untouched
  left_join(redlands_horse_metadata_tidy, by = "sample") %>%  #join the two tidy data frames by "sample"
  select(-sample)                                             #remove the sample column
  
#want to determine if the relationship between day and counts is linear for each microRNA

get_lm_se <- function(microRNA_counts_tidy){            #create a variable where we get standard error from lm()
  fit <- lm(counts ~ day, data = microRNA_counts_tidy)  #fit a linear model to y=counts, x=day using microRNA_counts_tidy
  data.frame(term = names(fit$coefficients),            #structure the data frame
             slope = fit$coefficients,                  #return slope coefficient
             se = summary(fit)$coefficient[,2])         #return standard error coefficient 
}                                                       #close nested function

microRNA_counts_tidy %>%                             
  group_by(gene) %>%                                    #apply model by gene
  do(get_lm_se(.))                                      #do anything function



    
    summarise(avg_counts = mean(counts), sd_counts = sd(counts)) %>% 
  mutate(cv = sd_counts / avg_counts) %>% 
  mutate(diff = avg_counts - lag(avg_counts, default = 0)) %>% 
  filter(gene = )

microRNA_counts_tidy %>% 
  spread()


#Modelling for some microRNAs of interest

#miR103
miR103 <- microRNA_counts_tidy %>% 
  filter(gene == "eca-miR-103")

miR103

ggplot(miR103, aes(day, counts)) +
  geom_point() +
  geom_smooth(model = glm) +
  scale_y_log10()

miR103_mod <- glm(counts ~ day(), data = miR103)
summary(miR103_mod)

coefficient <- coef(miR103_mod)
coefficient
slope <- coefficent  



mod1 <- lm(counts ~ ns(day, 1), data = miR103)
coef(miR103_mod)

#miR16

miR16 <- microRNA_counts_tidy %>% 
  filter(gene == "eca-miR-16")

miR16

ggplot(miR16, aes(day, counts)) +
  geom_point() +
  geom_smooth(model = lm)

miR16_mod <- lm(counts ~ day, data = miR16)
coef(miR16_mod)
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



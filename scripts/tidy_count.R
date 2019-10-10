.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(tidyverse)
library(cowplot)

read_csv("data/Redlands_counts.csv")
Redlands_counts <- read_csv("Data/Redlands_counts.csv")

read_csv("data/redlands_horse_metadata.csv")
redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv")


Redlands_counts

redlands_horse_metadata

#tidy the count data 
redlands_counts_tidy <- Redlands_counts %>%              #assign to new variable
  gather(sample, counts, -gene) %>%                      #gather to get sample names into a column
  left_join(redlands_horse_metadata, by = "sample") %>%  #join by sample names
  rename(day = condition) %>%                            #rename column 
  select(-sample) %>%                                    #remove sample column
  mutate()


#making some plots
plot_counts_day <- ggplot(data = redlands_counts_tidy) +
  geom_point( mapping = aes(x = day,
                            y = counts, 
                            color = animal))





.libPaths(c("C:/Users/ale097/Data School/Packages"))
library(tidyverse)
library(cowplot)

read_csv("data/Redlands_counts.csv")
Redlands_counts <- read_csv("Data/Redlands_counts.csv")

read_csv("data/redlands_horse_metadata.csv")
redlands_horse_metadata <- read_csv("Data/redlands_horse_metadata.csv")


Redlands_counts

redlands_horse_metadata

counts_morethan100 <- Redlands_counts %>% 
  mutate(sum_counts = s1 +s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14 + s15) %>% 
  filter(sum_counts >= 100)



plot_all <- ggplot(data = counts_morethan100,
                        mapping = aes(x = , y = Temp_min, group = )) + geom_point(alpha = 0.2)


?rowSums

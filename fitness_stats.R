
## install packages
#install.packages("tidyverse")

## "We tested for significant hybrid breakdown by comparing 
## fitness of the F2 hybrid with its maternal population, 
## using a Mann–Whitney U-test and an alpha of 0.05 in R 
## 2.15.1 (R Development Core Team). We tested for significant
## recovery in lines for which the mean reached or passed the 
## reference parental fitness, adjusting the P-value when 
## multiple comparisons occur at the same time."

#As you can see from Figure 2, where I plotted means +-1SE:
#  1- F2s are not significantly less fit than the respective parental;
#2- only some hybrid lines (probably SC4, SC7 and SC8) and significantly more fit than the F2 with the same mitochondria.
#But we need some p-values to please reviewers.

#Here is the detailed list of what we need:
#  - is there significant hybrid breakdown, in each mitochondrial background?
#  comparisons:
#  SC vs SC-F2;
#SD vs SD-F2
#simple Mann–Whitney U-test, without correction of p-value since it's only one comparison per pair.

## load packages
library(readr)
library(dplyr)
library(stringr)


setwd("/Users/tessa/dib-lab/2020-fitness-stats")

# read in data
fecundity_data <- readr::read_csv("fitness_data_fecundity.csv", col_types = cols(class = col_factor(NULL)))
survivorship_data <- read_csv("fitness_data_survivorship.csv",  col_types = cols(class = col_factor(NULL)))

# sanity check data
barplot(table(fecundity_data$class))
barplot(table(survivorship_data$class))
table(survivorship_data$class)


## Survivorship data ##

# Is there significant hybrid breakdown, in each mitochondrial background?
## SC vs SC-F2

### survivorship
sc_scf2 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SC", "SCF2"))

sc_scf2_wc <- wilcox.test(survivorship~ class,  
                   data = sc_scf2, 
                   exact = FALSE)

### fecundity
f_sc_scf2 <- fecundity_data %>% 
  select(class, fecundity) %>% 
  filter(class %in% c("SC", "SCF2"))

f_sc_scf2_wc <- wilcox.test(fecundity~ class,  
                            data = f_sc_scf2, 
                            exact = FALSE)

## SD vs SD-F2

## survivorship
sd_sdf2 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SD", "SDF2"))

sd_sdf2_wc <- wilcox.test(survivorship~ class,  
                          data = sd_sdf2, 
                          exact = FALSE)

# fecundity
f_sd_sdf2 <- fecundity_data %>% 
  select(class, fecundity) %>% 
  filter(class %in% c("SD", "SDF2"))

f_sd_sdf2_wc <- wilcox.test(fecundity~ class,  
                            data = f_sd_sdf2, 
                            exact = FALSE)


## Recovery from hybrid breakdown (in each mtDNA background)
# SCF2 vs SC5, SC6, SC7, SC8

scf2_sc4 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SCF2", "SC4"))
scf2_sc4_wc <- wilcox.test(survivorship~ class,  
                          data = scf2_sc4, 
                          exact = FALSE)
scf2_sc4_wcPadj <- p.adjust(scf2_sc4_wc$p.value, method="bonferroni", n=4)

scf2_sc6 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SCF2", "SC6"))
scf2_sc6_wc <- wilcox.test(survivorship~ class,  
                           data = scf2_sc6, 
                           exact = FALSE)
scf2_sc6_wcPadj <- p.adjust(scf2_sc6_wc$p.value, method="bonferroni", n=4)

scf2_sc7 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SCF2", "SC7"))
scf2_sc7_wc <- wilcox.test(survivorship~ class,  
                           data = scf2_sc7, 
                           exact = FALSE)
scf2_sc7_wcPadj <- p.adjust(scf2_sc7_wc$p.value, method="bonferroni", n=4)

scf2_sc8 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SCF2", "SC8"))
scf2_sc8_wc <- wilcox.test(survivorship~ class,  
                           data = scf2_sc8, 
                           exact = FALSE)
scf2_sc8_wcPadj <- p.adjust(scf2_sc8_wc$p.value, method="bonferroni", n=4)


#  SDF2 vs SD3, SD4, SD7
sdf2_sd3 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SDF2", "SD3"))
sdf2_sd3_wc <- wilcox.test(survivorship~ class,  
                           data = sdf2_sd3, 
                           exact = FALSE)
sdf2_sd3_wcPadj <- p.adjust(sdf2_sd3_wc$p.value, method="bonferroni", n=3)

sdf2_sd4 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SDF2", "SD4"))
sdf2_sd4_wc <- wilcox.test(survivorship~ class,  
                           data = sdf2_sd4, 
                           exact = FALSE)
sdf2_sd4_wcPadj <- p.adjust(sdf2_sd4_wc$p.value, method="bonferroni", n=3)

sdf2_sd7 <- survivorship_data %>% 
  select(class, survivorship) %>% 
  filter(class %in% c("SDF2", "SD7"))
sdf2_sd7_wc <- wilcox.test(survivorship~ class,  
                           data = sdf2_sd7, 
                           exact = FALSE)
sdf2_sd3_wcPadj <- p.adjust(sdf2_sd3_wc$p.value, method="bonferroni", n=3)






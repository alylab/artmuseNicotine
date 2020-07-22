## setup ----

require(tidyverse)
require(magrittr)

source(here::here("R", "paths.R"))
load(paste(stats_dir, "raw.rda", sep = "/"))

## construct PCA on nicotine dependence covariates ----

pca_nicotine <- demos %>%
  select(subj_num, ppm_on, ppm_off, cigs_per_day_est, ftnd, years_smoke) %>% 
  mutate(ppm_diff = ppm_on - ppm_off) %>% 
  select(-subj_num) %>% 
  filter(complete.cases(.)) %>% 
  prcomp(center = TRUE, scale = TRUE)

save(pca_nicotine, file = paste(stats_dir, "pca.rda", sep = "/"))

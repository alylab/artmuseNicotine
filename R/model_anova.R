## setup ----

require(tidyverse)
require(magrittr)
require(rlang)

source(here::here("R", "paths.R"))
load(paste(stats_dir, "raw.rda", sep = "/"))

## set up repeated measures anova ----

aprime_anova <- sdt_metrics %>% 
  # All valid trials, so this should be true
  mutate(probe = cue) %>% 
  aov(aprime ~ exptCond * probe * on_smoking + Error(subj_num), data = .)

save(aprime_anova, file = paste(stats_dir, "models_anova.rda", sep = "/"))

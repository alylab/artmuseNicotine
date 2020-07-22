## setup ----

require(tidyverse)
require(magrittr)
options(mc.cores = parallel::detectCores())
future::plan(future::multiprocess)

source(here::here("R", "paths.R"))
load(paste(stats_dir, "raw.rda", sep = "/"))

## resample the bootstrap iterations ----

boots_by_subj <- raw %>%
  filter(valid == 1) %>%
  select(-c(valid, session_num, session_order, cue)) %>%
  mutate(corResp = recode(corResp, `0` = "fa", `1` = "hit"),
         on_smoking = recode(on_smoking, `0` = "off", `1` = "on"),
         exptCond = recode(exptCond, `0` = "control", `1` = "relational")) %>%
  nest(trials = -c(subj_num, ppm, years_smoke, on_smoking, exptCond, probe, corResp)) %>%
  nest(conditions = -subj_num) %>%
  nest(subjects = everything()) %>%
  # the each arg sets the number of iterations, my friends!
  # note that here, 1:n() is 1:1
  slice(rep(1:n(), each = 500)) %>%
  mutate(iteration = 1:n(),
         # resample outer level: subjects
         subjects = furrr::future_map(subjects, ~sample_frac(.x, size = 1, replace = TRUE) %>%
                                        mutate(subj_num = 1:nrow(.)),
                                      .progress = TRUE)) %>%
  unnest(subjects) %>%
  unnest(conditions) %>%
  # CHANGE N BOOTSTRAPS HERE
  # resample inner level: trials (honoring condition structure)
  mutate(boots = furrr::future_map(trials, ~rsample::bootstraps(.x, times = 1) %>%
                                     mutate(rate = map_dbl(splits, ~mean(as.data.frame(.x)$resp))) %>% 
                                     select(rate),
                                   .progress = TRUE)) %>%
  select(-trials) %>%
  unnest(boots) %>%
  pivot_wider(names_from = corResp, values_from = rate, names_prefix = "rate_") %>%
  mutate(aprime = sdt_aprime(rate_hit, rate_fa),
         aprime = coalesce(aprime, 0.5)) %>%
  nest(data = -iteration)

## make aprime_by_ppm (regression on the RAW data) ----

aprime_by_ppm <- raw %>%
  filter(valid == 1) %>%
  mutate(corResp = recode(corResp, `0` = "fa", `1` = "hit"),
         on_smoking = recode(on_smoking, `0` = "off", `1` = "on"),
         exptCond = recode(exptCond, `0` = "control", `1` = "relational")) %>%
  group_by(subj_num, ppm, years_smoke, on_smoking, exptCond, probe, corResp) %>%
  summarize(rate = mean(resp)) %>%
  ungroup() %>% 
  pivot_wider(names_from = corResp, values_from = rate, names_prefix = "rate_") %>%
  mutate(aprime = sdt_aprime(rate_hit, rate_fa)) %>% 
  pivot_longer(cols = c(rate_hit, rate_fa, aprime),
               names_to = "metric_type",
               values_to = "value") %>% 
  pivot_wider(names_from = on_smoking, values_from = c(value, ppm)) %>% 
  mutate(ppm_diff = ppm_on - ppm_off,
         value_diff = value_on - value_off,
         # luckily for us, centering at 0.5 is relatively interpretable
         # for hit rate, false alarm rate, and A'
         # because 0.5 is chance in a balanced 2-choice task
         value_off_c = value_off - 0.5) %>%
  arrange(metric_type, exptCond, probe, value_off, value_diff) %>%
  group_by(metric_type, exptCond, probe) %>%
  # For graphing, if you want to graph on-off improvement by baseline level
  # for sussing out regression to the mean
  mutate(rank_value_off = 1:n()) %>%
  nest(data = -c(metric_type, exptCond, probe)) %>%
  mutate(resid_value_diff = map(data, ~lm(value_diff ~ value_off_c, data = .x) %>%
                                   broom::augment() %>%
                                   select(value_diff_resid = .resid)),
         resid_ppm_diff = map(data, ~lm(ppm_diff ~ value_off_c, data = .x) %>%
                                broom::augment() %>%
                                select(ppm_diff_resid = .resid)),
         data = pmap(list(data, resid_value_diff, resid_ppm_diff),
                     function(a, b, c) {bind_cols(a, b, c)}),
         model_main = map(data, ~lm(value_diff ~ value_off_c + ppm_diff, data = .x)),
         model_covar = map(data, ~lm(value_diff ~ value_off_c + ppm_diff + years_smoke, data = .x)),
         model_diffonly = map(data, ~lm(value_diff ~ ppm_diff, data = .x)),
         augs = map(model_main, ~.x %>%
                      broom::augment() %>%
                      select(value_diff_fit = .fitted)),
         resid_augs = map(data, ~lm(value_diff_resid ~ ppm_diff_resid, data = .x) %>%
                            broom::augment() %>%
                            select(value_diff_fit_resid = .fitted)),
         data = map2(data, augs, ~bind_cols(.x, .y)),
         data = map2(data, resid_augs, ~bind_cols(.x, .y))) %>%
  select(-starts_with("resid"), -augs)

## make aprime_by_ppm_boot ----

aprime_by_ppm_boot <- boots_by_subj %>%
  unnest(data) %>%
  pivot_longer(cols = c(rate_hit, rate_fa, aprime),
               names_to = "metric_type",
               values_to = "value") %>% 
  pivot_wider(names_from = on_smoking, values_from = c(value, ppm)) %>% 
  mutate(ppm_diff = ppm_on - ppm_off,
         value_diff = value_on - value_off,
         # luckily for us, centering at 0.5 is relatively interpretable
         # for hit rate, false alarm rate, and A'
         # because 0.5 is chance in a balanced 2-choice task
         value_off_c = value_off - 0.5) %>%
  nest(data = -c(iteration, metric_type, exptCond, probe)) %>%
  mutate(resid_value_diff = map(data, ~lm(value_diff ~ value_off_c, data = .x) %>% broom::augment() %>% select(value_diff_resid = .resid)),
         resid_ppm_diff = map(data, ~lm(ppm_diff ~ value_off_c, data = .x) %>% broom::augment() %>% select(ppm_diff_resid = .resid)),
         data = pmap(list(data, resid_value_diff, resid_ppm_diff), function(a, b, c) {bind_cols(a, b, c)}),
         model_resid = map(data, ~lm(value_diff_resid ~ ppm_diff_resid, data = .x)),
         coefs_main = map(data, ~lm(value_diff ~ value_off_c + ppm_diff, data = .x) %>%
                       broom::tidy()),
         coefs_covar = map(data, ~lm(value_diff ~ value_off_c + ppm_diff + years_smoke, data = .x) %>% 
                             broom::tidy()),
         coefs_diffonly = map(data, ~lm(value_diff ~ ppm_diff, data = .x) %>%
                             broom::tidy())) %>%
  left_join(aprime_by_ppm %>%
              select(metric_type,
                     exptCond,
                     probe,
                     data,
                     coefs_main = model_main,
                     coefs_covar = model_covar,
                     coefs_diffonly = model_diffonly) %>%
              mutate_at(vars(starts_with("coefs")),
                        ~map(., ~broom::tidy(.x))),
            by = c("metric_type", "exptCond", "probe"),
            suffix = c("_boot", "_raw")) %>%
  mutate(predicted_resid = map2(data_raw, model_resid,
                                ~.x %>%
                                  select(ppm_diff_resid) %>%
                                  predict(.y, newdata = .)),
         data_raw = map2(data_raw, predicted_resid,
                         ~.x %>%
                           mutate(obs = 1:nrow(.), predicted = .y))) %>%
  select(-ends_with("resid"), -data_boot)

## finish up ----

save(boots_by_subj, file = paste(stats_dir, "boots.rda", sep = "/"))

save(aprime_by_ppm, aprime_by_ppm_boot, file = paste(stats_dir, "aprime_by_ppm.rda", sep = "/"))

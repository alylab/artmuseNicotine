## setup ----

require(tidyverse)

source(here::here("R", "paths.R"))

snodgrass <- function (num, denom) return((num+.5)/(denom+1))
snodgrass_vec <- function (vec) return((sum(vec) + 0.5) / (length(vec) + 1))
sdt_pr <- function (hit, fa) return(hit - fa)
sdt_br <- function(hit, fa) return(fa / (1 - (hit - fa)))
sdt_aprime <- function (hit, fa) return(.5 + (sign(hit-fa)*((hit-fa)^2 + abs(hit-fa))/(4*pmax(hit,fa)-(4*hit*fa))))
sdt_b2prime <- function (hit, fa)  return(sign(hit-fa)*(hit*(1-hit) - fa*(1-fa))/(hit*(1-hit) + fa*(1-fa)))
sdt_dprime <- function (hit, fa) return(qnorm(hit) - qnorm(fa))
sdt_c <- function (hit, fa) return(-.5 * (qnorm(hit) + qnorm(fa)))

## read in and preprocess trialwise data from raw ----

demos = read_csv(paste(data_dir, "demo_smoke_habits.csv", sep = "/")) %>%
  mutate(cigs_per_day_est = 10 * cigs_per_day_est) %>%
  mutate_if(is.numeric, as.integer) %>%
  mutate(cigs_per_day_est = .1 * cigs_per_day_est,
         ppm_on = if_else(session_order == 1, ppm_s1, ppm_s2),
         ppm_off = if_else(session_order == 1, ppm_s2, ppm_s1))

# Read in all .txt files in the raw data folder
raw = tibble(filename = list.files(data_dir_task, pattern = ".txt", full.names = TRUE)) %>%
  mutate(data = map(filename, ~.x %>%
                      read_delim(delim = "\t",
                                 skip = 9,
                                 col_names = FALSE) %>%
                      # Skipping the "actual" header rows for easier reading-in,
                      # then manually renaming cols
                      select(trial = X1,
                             stimOnsetTime = X2,
                             exptCond = X4,
                             stimNum = X6,
                             cue = X7,
                             probe = X8,
                             valid = X9,
                             resp = X10,
                             corResp = X11,
                             acc = X12,
                             rt = X13) %>%
                      # Coerce cols to integer where appropriate to reduce memory size
                      mutate_at(c("trial", "exptCond", "valid", "resp", "corResp", "acc"),
                                as.integer) %>%
                      # recode this col to more informative labels for now
                      mutate_at(c("cue", "probe"), ~recode(., `1` = "art", `2` = "room")))) %>%
  # pull out just the subject and session numbers from the whole file name
  mutate(subj_num = as.integer(str_sub(filename, start = -10L, end = -8L)),
         session_num = str_sub(filename, start = -5L, end = -5L),
         session_num = if_else(session_num == "2", 2L, 1L)) %>%
  select(-filename) %>%
  # bind on desired demos
  # create "on_smoking" based on "session_order"
  # which indicates which session (1 or 2) was the smoking session
  left_join(demos %>%
              select(subj_num, session_order, ppm_on, ppm_off) %>%
              pivot_longer(names_to = "on_smoking", values_to = "ppm", cols = starts_with("ppm")) %>% 
              mutate(on_smoking = recode(on_smoking, ppm_on = 1L, ppm_off = 0L),
                     # move session_num to 0 and 1, I gotta plan
                     session_order_reverse = if_else(session_order == 1L, 2L, 1L),
                     session_num = if_else(on_smoking == 1, session_order, session_order_reverse),
                     ppm_z = c(scale(ppm))),
            by = c("subj_num", "session_num")) %>%
  select(-session_order_reverse) %>%
  # Scrub all subjects who don't appear in demographics data
  filter(!is.na(session_order)) %>%
  unnest(data)

## unmodeled summary statistics from raw ----

sdt_metrics <- raw %>%
  filter(valid == 1) %>%
  mutate(corResp = recode(corResp, `0` = "fa", `1` = "hit"),
         on_smoking = recode(on_smoking, `0` = "off", `1` = "on"),
         exptCond = recode(exptCond, `0` = "control", `1` = "relational")) %>%
  group_by(subj_num, on_smoking, exptCond, cue, corResp) %>%
  summarize(rate = mean(resp), snod = snodgrass_vec(resp), n_trials = n()) %>%
  pivot_wider(names_from = corResp, values_from = c(rate, snod)) %>%
  mutate(aprime = sdt_aprime(rate_hit, rate_fa),
         aprime_snod = sdt_aprime(snod_hit, snod_fa),
         dprime = sdt_dprime(snod_hit, snod_fa),
         pr = sdt_pr(rate_hit, rate_fa))

## finish it up ----

save(raw, sdt_metrics, file = paste(stats_dir, "raw.rda", sep = "/"))

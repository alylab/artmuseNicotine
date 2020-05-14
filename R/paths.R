# path strings for loading and saving data
data_dir <- here::here("ignore", "materials", "data")
data_dir_task <- here::here("ignore", "materials", "data", "task_data")
stats_dir <- here::here("ignore", "data_R")

# signal detection metric helper functions that occur in more than one script
snodgrass <- function (num, denom) return((num+.5)/(denom+1))
snodgrass_vec <- function (vec) return((sum(vec) + 0.5) / (length(vec) + 1))
sdt_pr <- function (hit, fa) return(hit - fa)
sdt_br <- function(hit, fa) return(fa / (1 - (hit - fa)))
sdt_aprime <- function (hit, fa) return(.5 + (sign(hit-fa)*((hit-fa)^2 + abs(hit-fa))/(4*pmax(hit,fa)-(4*hit*fa))))
sdt_b2prime <- function (hit, fa)  return(sign(hit-fa)*(hit*(1-hit) - fa*(1-fa))/(hit*(1-hit) + fa*(1-fa)))
sdt_dprime <- function (hit, fa) return(qnorm(hit) - qnorm(fa))
sdt_c <- function (hit, fa) return(-.5 * (qnorm(hit) + qnorm(fa)))

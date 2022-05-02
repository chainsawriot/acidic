require(tidyverse)
require(ggplot2)
require(furrr)
require(acidic)
plan(multisession)

ss <- c(500, 1000, 10000, 30000)
es <- c(0, .05, .1, .15, .2, .3, .4, .5, .75)
z <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9)
Q <- 0.7
rep_n <- 1000
tol <- 0.01

scenarios <- tibble::tibble(expand.grid(ss, es, z))
colnames(scenarios) <- c("ss", "es", "z")
res <- list()

for (i in sample(seq_len(nrow(scenarios)))) {
    print(i)
    res[[i]] <- furrr::future_map_dfr(rep(scenarios$es[i], rep_n), ~ simulate_phi(., n = scenarios$ss[i], z = scenarios$z[i], acidic:::gen_acc, tol = tol, Q = Q), .progress = TRUE, .options = furrr_options(seed = 123))
    saveRDS(res, "acc_res.RDS")
}

scenarios$res <- res
saveRDS(scenarios, "acc_scenarios.RDS")

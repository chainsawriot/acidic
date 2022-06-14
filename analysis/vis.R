require(tidyverse)
require(ggplot2)
require(acidic)
require(furrr)
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
    res[[i]] <- furrr::future_map_dfr(rep(scenarios$es[i], rep_n), ~ simulate_phi(., n = scenarios$ss[i], z = scenarios$z[i], acidic::gen_f1, tol = tol, Q = Q), .progress = TRUE, .options = furrr_options(seed = 123))
    saveRDS(res, "f1_res.RDS")
}

scenarios$res <- res
saveRDS(scenarios, "f1_scenarios.RDS")


scenarios <- readRDS("f1_scenarios.RDS")

scenarios %>% mutate(mdn = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.5)), pct25 = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.025)), pct975 = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.975))) %>% filter(pct25 < 0) %>% print(n = 1000)

scenarios <- readRDS("acc_scenarios.RDS")
scenarios %>% mutate(mdn = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.5, na.rm = TRUE)), pct25 = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.025, na.rm = TRUE)), pct975 = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.975, na.rm = TRUE))) %>% ggplot(aes(y = mdn, ymin = pct25, ymax = pct975, x = z)) + geom_line() + geom_ribbon(alpha = 0.5) + facet_grid(ss~es) + theme_minimal() + xlab("F1") + ylab("Observed effect size")

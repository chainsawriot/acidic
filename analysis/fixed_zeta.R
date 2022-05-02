require(acidic)
require(furrr)

rho <- .2

s <- c(0, .25, .5, .75, .9)

scenarios <- tibble::tibble(expand.grid(rho, s, s))

colnames(scenarios) <- c("es", "sen", "spec")
res <- list()
Q <- 0.7
tol <- 0.01
rep_n <- 1000
res <- list()

for (i in sample(seq_len(nrow(scenarios)))) {
    print(i)
    res[[i]] <- furrr::future_map_dfr(rep(scenarios$es[i], rep_n), ~ simulate_phi(., n = 10000, tp = scenarios$sen[i], tn = scenarios$spec[i], tol = tol, Q = Q, z = NULL, fun = NULL), .progress = TRUE, .options = furrr_options(seed = 123))
    saveRDS(res, "fixed_res.RDS")
}

scenarios$res <- res
saveRDS(scenarios, "fixed_scenarios.RDS")


require(tidyverse)
scenarios <- readRDS("fixed_scenarios.RDS")
scenarios %>% mutate(mdn = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.5, na.rm = TRUE)), pct25 = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.025, na.rm = TRUE)), pct975 = purrr::map_dbl(res, ~quantile(.$d2_phi, probs = 0.975, na.rm = TRUE))) %>% ggplot(aes(y = mdn, ymin = pct25, ymax = pct975, x = spec)) + geom_line() + geom_ribbon(alpha = 0.5) + facet_grid(.~sen) + theme_minimal() + xlab("Specificity") + ylab("Observed effect size (R)") + theme(axis.text.x = element_text(angle = 90)) + ylim(-.2, +.2)

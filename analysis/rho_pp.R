require(acidic)

total_n <- 5000
rho <- .1

cal_pp <- function(rho, total_n = 5000, tol = .01) {
    x <- NA
    while(length(x) == 1) {
        x <- estimate_phi(rho, n = total_n, tol = tol, return_data = FALSE)
    }
    pp <- (x[1] + x[3]) / total_n
    return(pp)
}


rho <- seq(0, 1, by = .02)
set.seed(1212121)
res <- purrr::map(rho, ~replicate(300, cal_pp(.)))
saveRDS(purrr::map2_dfr(rho, res, function(x, y) data.frame(rho = x, pp = y)), "pp_sim.RDS")

pp_sim <- readRDS("pp_sim.RDS")

require(ggplot2)
require(magrittr)
pp_sim %>% ggplot(aes(x = rho, y = pp)) + geom_point(alpha = 0.3) + stat_smooth(method = "loess", formula = "y ~ x") + xlab("True effect size") + ylab("Prevalence") + theme_minimal() -> sim_fig
ggsave("sim_fig.pdf", sim_fig, width = 4, height = 4)

## require(tidyverse)
## require(ggplot2)
## source("emulator.R")
## require(furrr)
## plan(multisession)

## rep_n <- 1000
## es <- .1111
## tol <- 0.01
## n <- 831

simulate_perfect <- function(phi, n, tol = 0.01, max_iters = 100) {
    data <- NA
    while(!is.data.frame(data)) {
        data <- estimate_phi(phi, n, tol = tol, max_iters = max_iters)
    }
    return(cor.test(as.numeric(data$x), as.numeric(data$y))$p.value)
}

simulate_distort <- function(phi, n, tp, tn, Q = 0.7, tol = 0.03) {
    data <- NA
    while(!is.data.frame(data)) {
        data <- estimate_phi(phi, n, tol = tol)
    }
    dx1 <- distort_gt(data$x, tp = tp, tn = tn)
    dx2 <- distort_gt(dx1, Q = Q)
    return(cor.test(as.numeric(dx2), as.numeric(data$y))$p.value)
}

## simulate_perfect(0.1111, 852, tol = 0.01)

## expected_p <- replicate(rep_n, simulate_perfect(phi = es, n = n, tol = tol))
## beta <- sum(expected_p > 0.05) / rep_n
## (ori_power <- 1 - beta)

## set.seed(12121)
## obs_p <- replicate(rep_n, simulate_distort(es, n = 852, tp = .831, tn = .721, Q = .928, tol = tol))
## beta <- sum(obs_p > 0.05) / rep_n
## (now_power <- 1 - beta)

.gen_p <- function(prob, n, rep_n) {
    if (n == Inf) {
        return(rep(prob, rep_n))
    } else {
        return(rbinom(rep_n, n, prob) / n)
    }
}

cal_bigphi <- function(rho, n, tpr, tnr, Q, alpha = 0.05, tol = 0.01, rep_n = 2000, pn = Inf, nn = Inf, Qn = Inf, seed = sample(1:2^15, 1), .progress = TRUE) {
    set.seed(seed)
    esv <- rep(rho, rep_n)
    tpv <- .gen_p(tpr, pn, rep_n)
    tnv <- .gen_p(tnr, nn, rep_n)
    Qv <- .gen_p(Q, Qn, rep_n)
    args <- list(phi = esv, tp = tpv, tn = tnv, Q = Qv, tol = tol, n = n)
    obs_p <- furrr::future_pmap_dbl(args, simulate_distort, .options = furrr::furrr_options(seed = seed), .progress = .progress)
    beta <- sum(obs_p > alpha) / rep_n
    return(1 - beta)
}

## cal_bigphi(rho = 0.1, n = 15000, tp = 0.8, tn = 0.6, Q = 0.8)

## cal_bigphi(rho = 0.1, n = 15000, tp = 0.6, tn = 0.6, Q = 0.6)

## cal_bigphi(rho = 0.1, n = 15000, tp = 0.6, tn = 0.6, Q = 0.8, pn = 300, nn = 600, Qn = 100, rep_n = 2000)

## cal_bigphi(rho = 0.1111, n = 852, tp = 0.831, tn = 0.721, Q = 0.928, seed = 123, rep_n = 3000)

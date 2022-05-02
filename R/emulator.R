## Simulation level parameters

## Model level parameters
## TP
## TN

## tp <- 0.3
## tn <- 0.7

cal_f1 <- function(tp, tn) {
    (2 * tp) / (2 * (tp + (1 - tp) + (1 - tn)))
}

gen_ppv <- function(f1, tp) {
    (f1 * tp) / (2 * tp - f1)
}

## f1 <- 0.6857143
## tp <- 0.75
## ppv <- gen_ppv(f1, tp)
## pp <- (120 + 40) / (120 + 40 + 70 + 170)
## 1 - (tp * pp * (1 - ppv)) / (ppv * (1-pp))

#' generate a random combination of tp and tn, given f1. The combination generated is always >= 0 and <= 1.
#' @export
gen_f1 <- function(f1, pp, tp = runif(1)) {
    tn <- -1
    while(tn < 0 | tn > 1) {
        ppv <- gen_ppv(f1, tp)
        tn <- 1 - (tp * pp * (1 - ppv)) / (ppv * (1 - pp))
        if (tn < 0 | tn > 1) {
            tp <- runif(1)
        }
    }
    return(tibble::tibble(tp = tp, tn = tn))
}

#' generate a random combination of tp and tn, given f1. The combination generated is always >= 0 and <= 1.
## gen_f1 <- function(f1, tp = runif(1)) {
##     tn <- -1
##     while(tn < 0 | tn > 1) {
##         tn <- tp + 2 - ((2 * tp) / f1)
##         if (tn < 0 | tn > 1) {
##             tp <- runif(1)
##         }
##     }
##     return(tibble::tibble(tp = tp, tn = tn))
## }

#' @export
gen_acc <- function(acc, pp, tp = runif(1)) {
    tn <- -1
    while(tn < 0 | tn > 1) {
        tn <- (acc - (pp * tp)) / (1 - pp)
        if (tn < 0 | tn > 1) {
            tp <- runif(1)
        }
    }
    return(tibble::tibble(tp = tp, tn = tn))
}

.emulate_ml <- function(ground_truth, tp = 1, tn = 1) {
    if (ground_truth) {
        if (runif(1) > tp) {
            return(!ground_truth)
        }
    } else {
        if (runif(1) > tn) {
            return(!ground_truth)
        }
    }
    return(ground_truth)
}

.emulate_coding <- function(ground_truth, Q = 1) {
    if (runif(1) > Q) {
        return(sample(c(TRUE, FALSE), size = 1, replace = TRUE))
    } else {
        return(ground_truth)
    }
}

gen_x <- function(prob, size = 1) {
    if (length(prob) == 1) {
        sample(c(TRUE, FALSE), size = size, prob = c(prob, (1-prob)), replace = TRUE)
    } else {
        purrr::map_lgl(prob, ~sample(c(TRUE, FALSE), size = 1, replace = TRUE, prob = c(., 1-.)))
    }
}

#' @examples
#' mod <- gen_acc(.9)
#' ground_truth <- gen_x(0.3, size = 100)
#' table(ground_truth, distort_gt(ground_truth, tp = mod$tp, tn = mod$tn))  
distort_gt <- function(x, tp = 1, tn = 1, Q = NULL) {
    if (!is.null(Q)) {
        return(purrr::map_lgl(x, .emulate_coding, Q = Q))
    }
    purrr::map_lgl(x, .emulate_ml, tp = tp, tn = tn)
}

## Vectorize
cal_p <- function(x, beta1, beta0 = 0) {
    vt <- beta0 + (beta1 * x)
    1 / (1 + exp(-vt))
}

## simulate_f1 <- function(f1, beta1, beta0, prev, n) {
##     x <- gen_x(prev, n)
##     p <- cal_p(x, beta1, beta0)
##     y <- gen_x(p)
##     mod <- gen_f1(f1)
##     dy <- distort_gt(y, tp = mod$tp, tn = mod$tn)
##     c1 <- coef(glm(y~x, family = binomial))
##     c2 <- coef(glm(dy~x, family = binomial))
##     tibble::tibble(f1 = f1, beta1 = beta1, beta0 = beta0, prev = prev, n = n, obeta1 = c1[2], obeta0 = c1[1], dbeta1 = c2[2], dbeta0 = c2[1])
## }

#' @param z an accuarcy metric
#' @param fun a function that generate TP and TN
#' @param beta1 the preset beta1 of the logistic regression model
#' @param beta0 the preset beta0 of the logistic regression model
#' @param prev the prevalence of the positive case
#' @param n the size of the generated data
simulate_pred <- function(z, fun, x, Q = 0.7, tp = NULL, tn = NULL) {
    pp <- sum(x$x) / nrow(x)
    if (is.null(tp) & is.null(tn)) {
        mod <- fun(z, pp)
        tp <- mod$tp
        tn <- mod$tn
    }
    dx1 <- distort_gt(x$x, tp = tp, tn = tn)
    dx2 <- distort_gt(dx1, Q = Q)
    d1_phi <- cal_phixy(dx1, x$y)
    d2_phi <- cal_phixy(dx2, x$y)
    o_phi <- cal_phixy(x$x, x$y)
    tibble::tibble(z = z, o_phi = o_phi, d1_phi = d1_phi, d2_phi = d2_phi)
}

## simulate_f1 <- function(f1, beta1, beta0, prev, n) {
##     simulate_pred(z = f1, fun = gen_f1, beta1 = beta1, beta0 = beta0, prev = prev, n = n)
## }

## res <- dplyr::bind_rows(purrr::rerun(1000, simulate_pred(0.99, beta1 = 1, beta0 = 0, prev = 0.3, fun = gen_f1, n = 3000)))


## a <- 6
## b <- 2
## c <- 1
## d <- 3
## n <- 12

cal_phi <- function(a, b, c, d) {
    nom <- (a * d - b * c)
    ## The as.numeric is to prevent integer overflow
    denom <- sqrt(as.numeric((a + c)) * (a + b) * (d + c) * (d + b))
    nom / denom
}

cal_phixy <- function(x, y) {
    va <- sum(x & y)
    vb <- sum(x & !y)
    vc <- sum(!x & y)
    vd <- sum(!x & !y)
    cal_phi(va, vb, vc, vd)
}

## parameters
## n <- 1000
## phi <- 0.4
## tol <- 0.03

gen_data <- function(va, vb, vc, vd) {
    dfa <- data.frame(x = rep(TRUE, va), y = rep(TRUE, va))
    dfb <- data.frame(x = rep(TRUE, vb), y = rep(FALSE, vb))
    dfc <- data.frame(x = rep(FALSE, vc), y = rep(TRUE, vc))
    dfd <- data.frame(x = rep(FALSE, vd), y = rep(FALSE, vd))
    df <- rbind(dfa, dfb, dfc, dfd)
    df <- df[sample(seq_len(nrow(df))), ]
    rownames(df) <- NULL
    return(df)
}

#' @export
estimate_phi <- function(phi, n, tol = 0.03, max_iters = 100, return_data = TRUE) {
    prev_cd <- runif(1, 0.1, 0.5)
    prev_c_in_cd <- runif(1, 0.01, 0.5)
    total_cd <- floor(n * prev_cd)
    vc <- floor(total_cd * prev_c_in_cd)
    vd <- total_cd - vc
    total_ab <- n - total_cd
    ## initial guess
    prev_a_in_ab_min <- prev_c_in_cd ## it must be higher than this.
    prev_a_in_ab_max <- 1
    diff_phi <- -100
    iter <- 0
    while(abs(diff_phi) > tol) {
        prev_a_in_ab <- mean(c(prev_a_in_ab_max, prev_a_in_ab_min))
        ##print(prev_a_in_ab)
        va <- floor(total_ab * prev_a_in_ab)
        vb <- total_ab - va
        cur_phi <- cal_phi(va, vb, vc, vd)
        diff_phi <- cur_phi - phi
        ##print(diff_phi)
        ##print(iter)
        if (abs(diff_phi) > tol) {
            if (diff_phi < 0) {
                prev_a_in_ab_min <- prev_a_in_ab
            } else {
                prev_a_in_ab_max <- prev_a_in_ab
            }
        }
        iter <- iter + 1
        if (iter > max_iters) {
            return(NA)
        }
    }
    if (return_data) {
        return(gen_data(va, vb, vc, vd))
    } else {
        return(c(va, vb, vc, vd))
    }
}

#' @export
simulate_phi <- function(phi, n, z, fun, Q = 0.7, tol = 0.03, max_iters = 100, tp = NULL, tn = NULL) {
    data <- NA
    while(!is.data.frame(data)) {
        data <- estimate_phi(phi, n, tol = tol, max_iters = max_iters)
    }
    simulate_pred(z = z, fun = fun, x = data, Q = Q, tp = tp, tn = tn)
}


## for (i in 1:5) {
##     print(i)
##     res[[i]] <- furrr::future_map_dfr(rep(scenarios$es[i], 1000), ~ simulate_phi(., n = scenarios$ss[i], z = scenarios$z[i], gen_f1, tol = 0.01), .progress = TRUE, .options = furrr_options(seed = 123))
##     saveRDS(res, "f1_res.RDS")
## }

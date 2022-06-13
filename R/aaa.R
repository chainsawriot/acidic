## Style guide: tpr, tnr (true positive rate, true negative rate)
## pp (prevalence of positive), ppv (positive predictive value)

## A B
## C D

## To prevent problems of using "c", use "va", "vb", "vc", and "vd"

cal_f1 <- function(tpr, tnr, pp) {
    ppv <- (tpr * pp) / (tnr * pp + (1 - tnr) * (1 - pp))
    2 * (tpr + ppv) / (tpr * ppv)
}

cal_acc <- function(tpr, tnr, pp) {
    pp * tpr + (1 - pp) * tnr
}

.gen_ppv <- function(f1, tpr) {
    (f1 * tpr) / (2 * tpr - f1)
}

#' generate a random combination of tp and tn, given f1. The combination generated is always >= 0 and <= 1.
#' @export
gen_f1 <- function(f1, pp, tpr = runif(1)) {
    tnr <- -1
    while(tnr < 0 | tnr > 1) {
        ppv <- .gen_ppv(f1, tpr)
        tnr <- 1 - (tpr * pp * (1 - ppv)) / (ppv * (1 - pp))
        if (tnr < 0 | tnr > 1) {
            tpr <- runif(1)
        }
    }
    return(tibble::tibble(tpr = tpr, tnr = tnr, pp = pp, f1 = f1))
}

#' @export
gen_acc <- function(acc, pp, tpr = runif(1)) {
    tnr <- -1
    while(tnr < 0 | tnr > 1) {
        tnr <- (acc - (pp * tpr)) / (1 - pp)
        if (tnr < 0 | tnr > 1) {
            tpr <- runif(1)
        }
    }
    return(tibble::tibble(tpr = tpr, tnr = tnr, pp = pp, acc = acc))
}

.emulate_ml <- function(ground_truth, tpr = 1, tnr = 1) {
    if (ground_truth) {
        if (runif(1) > tpr) {
            return(!ground_truth)
        }
    } else {
        if (runif(1) > tnr) {
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
#' mod <- gen_acc(.9, .3)
#' ground_truth <- gen_x(0.3, size = 100)
#' table(ground_truth, distort_gt(ground_truth, tp = mod$tp, tn = mod$tn))  
distort_gt <- function(x, tpr = 1, tnr = 1, Q = NULL) {
    if (!is.null(Q)) {
        return(purrr::map_lgl(x, .emulate_coding, Q = Q))
    }
    purrr::map_lgl(x, .emulate_ml, tpr = tpr, tnr = tnr)
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
#' @param x a data.frame of generated dataset
#' @param Q coder agreement
#' @param tpr true positive rate
#' @param tnr true negatve rate
simulate_pred <- function(z, fun, x, Q = 0.7, tpr = NULL, tnr = NULL) {
    pp <- sum(x$x) / nrow(x)
    if (is.null(tpr) & is.null(tnr)) {
        mod <- fun(z, pp)
        tpr <- mod$tpr
        tnr <- mod$tnr
    }
    dx1 <- distort_gt(x$x, tpr = tpr, tnr = tnr)
    dx2 <- distort_gt(dx1, Q = Q)
    d1_phi <- cal_phixy(dx1, x$y)
    d2_phi <- cal_phixy(dx2, x$y)
    o_phi <- cal_phixy(x$x, x$y)
    tibble::tibble(z = z, o_phi = o_phi, d1_phi = d1_phi, d2_phi = d2_phi)
}

cal_phi <- function(va, vb, vc, vd) {
    nom <- (va * vd - vb * vc)
    ## The as.numeric is to prevent integer overflow
    denom <- sqrt(as.numeric((va + vc)) * (va + vb) * (vd + vc) * (vd + vb))
    nom / denom
}

cal_phixy <- function(x, y) {
    va <- sum(x & y)
    vb <- sum(x & !y)
    vc <- sum(!x & y)
    vd <- sum(!x & !y)
    cal_phi(va, vb, vc, vd)
}

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
simulate_phi <- function(phi, n, z, fun, Q = 0.7, tol = 0.03, max_iters = 100, tpr = NULL, tnr = NULL) {
    data <- NA
    while(!is.data.frame(data)) {
        data <- estimate_phi(phi, n, tol = tol, max_iters = max_iters)
    }
    simulate_pred(z = z, fun = fun, x = data, Q = Q, tpr = tpr, tnr = tnr)
}

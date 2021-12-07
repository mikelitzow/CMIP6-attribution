## Model checking utilities
## Michael Malick

## diagnostic checks ---------------------------------------
rhat_highest <- function(stanfit, k = 4, pars) {
    rhat <- get_rhat(stanfit, pars = pars)
    rhat.max <- rev(sort(rhat)[(length(rhat) - k):length(rhat)])
    return(rhat.max)
}

neff_lowest <- function(stanfit, k = 4, pars) {
    neff <- get_neff(stanfit, pars = pars)
    neff.min <- sort(neff)[1:k]
    return(neff.min)
}

pairs_lowest <- function(stanfit, k = 4, pars) {
    n <- get_neff(stanfit, pars = pars)
    n.min <- names(sort(n))[1:k]
    pairs(stanfit, pars = n.min)
}

get_rhat <- function(stanfit, pars) {
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    summary(stanfit, pars = pars)$summary[ , "Rhat"]
}

get_neff <- function(stanfit, pars) {
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    summary(stanfit, pars = pars)$summary[ , "n_eff"]
}

total_draws <- function(stanfit) {
    ## N chains * N draws  -- post warmup
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    dim(stanfit)[1] * dim(stanfit)[2]
}


## trace_plot ----------------------------------------------
trace_plot <- function(stanfit, ...) {
    ## Plot trace for each parameter
    ##
    ## stanfit = stanfit object
    ## ... = passed to rstan::extract
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }

    draws     <- extract(stanfit, permuted = FALSE, ...)
    warmup    <- stanfit@stan_args[[1]]$warmup
    thin      <- stanfit@stan_args[[1]]$thin
    nchains   <- dim(stanfit)[2]
    niter     <- dim(stanfit)[1]
    colors    <- 2:(nchains + 1)
    par_names <- names(draws[1, 1, ])
    xvals     <- (warmup + 1):(niter + warmup)

    for(i in seq_along(par_names)) {
        matplot(xvals, draws[, , i],
                type = "l",
                lty = 1,
                col = colors,
                xlab = "Iteration",
                ylab = "Value",
                main = par_names[i],
                axes = FALSE)
        graphics::box(col = "grey50")
        graphics::axis(1, lwd = 0, col = "grey50", lwd.ticks = 1)
        graphics::axis(2, lwd = 0, col = "grey50", las = 1, lwd.ticks = 1)
        txt <- paste("warmup =", warmup,
                     "   niter =", niter,
                     "   nchain =", nchains,
                     "   thin =", thin)
        graphics::mtext(txt, side = 3, cex = 0.7)
    }

}

if(FALSE) {

    ecode <- '
      parameters {
        real<lower=0> y[2];
      }
      model {
        y ~ exponential(1);
      }
    '
    fit <- stan(model_code = ecode, iter = 1000, chains = 2, thin = 3)
    trace_plot(fit)

}


## brms_R2_plot --------------------------------------------
brms_R2_plot <- function(fit, breaks = 50) {
    ## Plot posterior for R^2
    ##
    ## fit = a brmsfit object
    ## breaks = see ?hist

    if(!is(fit, "brmsfit"))
        stop("Input not a 'brmsfit' object")
    if(length(fit$criteria$bayes_R2) == 0)
        stop("R2 not saved in model object")

    r2 <- as.vector(fit$criteria$bayes_R2)
    x <- hist(r2,
              breaks = breaks,
         main = paste("Median R^2 =", round(median(r2), 2)),
         xlab = "R^2 value",
         axes = FALSE,
         pch = 16,
         col = "grey75",
         border = "white")
    axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey65")
    axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey65")
    box(col = "grey70")
    abline(v = median(r2), col = "red3", lwd = 2)
    abline(v = quantile(r2, probs = 0.025), col = "steelblue", lty = 2, lwd = 2)
    abline(v = quantile(r2, probs = 0.975), col = "steelblue", lty = 2, lwd = 2)
    rug(r2, col  = "wheat3")

}


## plot_y_yhat ---------------------------------------------
plot_y_yhat <- function(y, yhat,
                        ylab = "Response",
                        xlab = "Fitted value",
                        title = NULL) {
    ## Plot response ~ fitted-values
    ##
    ## y = observed response
    ## yhat = fitted values

    lim <- range(c(y, yhat))
    plot(yhat, y, pch = 19, col = "grey40",
         xlim = lim, ylim = lim,
         ylab = ylab,
         xlab = xlab,
         axes = FALSE,
         main = title,
         panel.first = grid(lty = 1, col = "grey90"))
    abline(0, 1, col = "red3", lty = 2, lwd = 2)
    rug(yhat, side = 1, col = "grey50")
    rug(y, side = 2, col = "grey50")
    box(col = "grey50")
    axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey50")
    axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey50")
}

if(FALSE) {

    y <- rnorm(100)
    x <- runif(100)
    fit <- lm(y ~ x)
    plot_y_yhat(y, fitted.values(fit))

}


## plot_resid_yhat -----------------------------------------
plot_resid_yhat <- function(resid, yhat,
                            ylab = "Residual",
                            xlab = "Fitted value",
                            loess = TRUE,
                            title = NULL) {
    ## Plot residuals ~ fitted-values
    ##
    ## resid = model residuals
    ## yhat = fitted values
    ## loess = logical, should a LOESS model fit be plotted

    plot(yhat, resid, pch = 19, col = "grey40",
         ylab = ylab,
         xlab = xlab,
         axes = FALSE,
         main = title,
         panel.first = grid(lty = 1, col = "grey90"))
    abline(h = 0, col = "grey50", lty = 2, lwd = 2)
    if(loess) {
        lw <- loess.smooth(yhat, resid)
        lines(lw$x, lw$y, col = "red3", lwd = 3)
    }
    rug(resid, side = 2, col = "grey50")
    rug(yhat, side = 1, col = "grey50")
    box(col = "grey50")
    axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey50")
    axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey50")
}

if(FALSE) {

    y <- rnorm(100)
    x <- runif(100)
    fit <- lm(y ~ x)
    plot_resid_yhat(residuals(fit), fitted.values(fit))
    plot_resid_yhat(residuals(fit), fitted.values(fit), loess = FALSE)

}


## brms_diag_plot ------------------------------------------
brms_diag_plot <- function(fit, coda = FALSE, yhat.loess = TRUE) {
    ## Wrapper for brms diagnostic graphic functions
    ##
    ## fit = a brmsfit object
    ## coda = logical, should coda MCMC diagnostic plots be created
    ## yhat.loess = logical, should a LOESS model fit be plotted for yhat diag

    if(!is(fit, "brmsfit"))
        stop("Input not a 'brmsfit' object")

    mcmc <- rstan::As.mcmc.list(fit$fit)
    np   <- brms::nuts_params(fit)
    lp   <- brms::log_posterior(fit)
    neff <- summary(fit$fit)$summary[ , "n_eff"]
    rhat <- summary(fit$fit)$summary[ , "Rhat"]
    max.neff <- coda::niter(mcmc) * coda::nchain(mcmc)

    print(bayesplot::mcmc_neff(bayesplot::neff_ratio(fit$fit)) +
          ggtitle(paste("Min Neff =", round(min(neff), 1), "/", max.neff)))
    print(bayesplot::mcmc_rhat(bayesplot::rhat(fit$fit)) +
          ggtitle(paste("Max Rhat =", round(max(rhat), 3))))
    print(bayesplot::mcmc_nuts_divergence(np, lp))

    if(length(fit$criteria$loo) != 0) {
        plot(fit$criteria$loo, diagnostic = "k")
    }

    if(length(fit$criteria$bayes_R2) != 0) {
        brms_R2_plot(fit)
    }

    form <- as.formula(fit$formula)
    y <- fit$data[[all.vars(form)[1]]]
    fit.v <- fitted(fit, re_formula = NULL)
    res <- residuals(fit)

    plot_y_yhat(y, fit.v[ , "Estimate"],
                title = "Response vs. fitted values",
                ylab = "Response",
                xlab = "Mean fitted value")

    plot_resid_yhat(res[ , "Estimate"], fit.v[ , "Estimate"],
                    title = "Realized residuals vs. fitted values",
                    loess = yhat.loess,
                    xlab = "Average fitted y",
                    ylab = "Average realized residual")

    print(brms::pp_check(fit, type = "dens_overlay", nsamples = 50))
    print(brms::pp_check(fit, type = "error_scatter_avg", nsamples = 50))
    print(brms::pp_check(fit, type = "boxplot", nsamples = 10))
    print(brms::pp_check(fit, type = "scatter_avg", nsamples = NULL))
    print(brms::pp_check(fit, type = "stat_2d", nsamples = NULL))

    if(coda) {
        codatools::coda_diag(mcmc)
    }
}

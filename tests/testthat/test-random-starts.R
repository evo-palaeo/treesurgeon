test_that("multiple random starts return the best finite fit", {
    model <- make_constraint_test_model()

    set.seed(2026)
    single_fit <- treesurgeon::fit_timescale_ml(
        model,
        method = "Nelder-Mead",
        control = list(maxit = 2000, reltol = 1e-8)
    )

    set.seed(2026)
    multi_fit <- treesurgeon::fit_timescale_ml(
        model,
        method = "Nelder-Mead",
        control = list(maxit = 2000, reltol = 1e-8),
        n_starts = 3
    )

    expect_s3_class(multi_fit, "timescale_ml_fit")
    expect_equal(multi_fit$n_starts, 3L)
    expect_true(multi_fit$best_start %in% seq_len(3))
    expect_true(is.finite(multi_fit$logLik))
    expect_gte(multi_fit$logLik, single_fit$logLik - 1e-8)
})


test_that("n_starts must be a positive integer", {
    model <- make_constraint_test_model()
    bad_values <- list(0, -1, 1.5, NA_real_, Inf, c(1, 2), numeric(0))

    for (value in bad_values) {
        expect_error(
            treesurgeon::fit_timescale_ml(model, n_starts = value),
            "n_starts"
        )
    }
})


test_that("multiple starts reach the best optimum found by an independent grid search", {
    model <- make_constraint_test_model(sigma_min = 0.01, root_max_age = 20)
    control <- list(maxit = 5000, reltol = 1e-10)

    # Grid over the feasible region. Each node is placed at a fraction of its
    # feasible interval, so every grid point respects the hard age bounds.
    fractions <- seq(0.01, 0.99, length.out = 60)
    grid <- expand.grid(rep(list(fractions), model$n_node))
    grid_values <- apply(grid, 1, function(u) {
        treesurgeon:::objective(stats::qlogis(u), model)
    })
    expect_true(any(is.finite(grid_values)))

    # Polish the best grid points to obtain a reference optimum.
    best_rows <- order(grid_values)[seq_len(5)]
    reference <- vapply(best_rows, function(i) {
        stats::optim(
            stats::qlogis(unlist(grid[i, ])),
            treesurgeon:::objective,
            model = model,
            method = "Nelder-Mead",
            control = control
        )$value
    }, numeric(1))
    reference_loglik <- -min(reference)

    for (seed in c(1, 2, 3)) {
        set.seed(seed)
        fit <- suppressMessages(treesurgeon::fit_timescale_ml(
            model,
            method = "Nelder-Mead",
            control = control,
            n_starts = 10L
        ))

        expect_gte(fit$logLik, reference_loglik - 1e-6)
    }
})


test_that("a non-converged fit warns and a converged fit does not", {
    model <- make_constraint_test_model()

    expect_warning(
        suppressMessages(treesurgeon::fit_timescale_ml(
            model,
            method = "Nelder-Mead",
            control = list(maxit = 2)
        )),
        "did not converge"
    )

    fit <- withCallingHandlers(
        treesurgeon::fit_timescale_ml(
            model,
            method = "Nelder-Mead",
            control = list(maxit = 5000, reltol = 1e-10)
        ),
        message = function(m) invokeRestart("muffleMessage")
    )
    expect_identical(as.integer(fit$convergence), 0L)
})

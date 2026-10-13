test_that("defaults use a small variance floor but leave the root unbounded", {
    tree <- make_constraint_test_tree()
    model <- treesurgeon::timescale_ml(
        tree = tree,
        tip_ages = rep(0, ape::Ntip(tree)),
        calibrations = make_constraint_test_calibrations(tree)
    )

    expect_equal(model$sigma_min, 0.01)
    expect_null(model$root_max_age)

    # Without a root maximum, theta[1] retains the exponential transform:
    # root age = minimum permissible root age + exp(theta[1]).
    theta <- rep(0, model$n_node)
    ages <- treesurgeon:::theta_to_node_ages(theta, model)
    expect_equal(
        unname(ages[model$root]),
        unname(model$minimum_ages[model$root] + 1),
        tolerance = 1e-12
    )
})


test_that("sigma_min = NULL explicitly disables the variance floor", {
    model <- make_constraint_test_model(sigma_min = NULL)
    expect_null(model$sigma_min)
})


test_that("valid constraint arguments are stored on the model", {
    model <- make_constraint_test_model(sigma_min = 0.05, root_max_age = 20)

    expect_equal(model$sigma_min, 0.05)
    expect_equal(model$root_max_age, 20)
})


test_that("sigma_min must be one finite positive number or NULL", {
    bad_values <- list(0, -0.1, NA_real_, Inf, -Inf, c(0.1, 0.2), numeric(0))

    for (value in bad_values) {
        expect_error(
            make_constraint_test_model(sigma_min = value),
            "sigma_min"
        )
    }
})


test_that("root_max_age must be one finite number or NULL", {
    bad_values <- list(NA_real_, Inf, -Inf, c(10, 20), numeric(0))

    for (value in bad_values) {
        expect_error(
            make_constraint_test_model(root_max_age = value),
            "root_max_age"
        )
    }
})


test_that("root_max_age must exceed the minimum permissible root age", {
    # The calibration fixture gives the non-root internal node a hard minimum
    # age of 1, which propagates to a minimum permissible root age of 1.
    expect_error(
        make_constraint_test_model(root_max_age = 1),
        "minimum permissible root age"
    )
    expect_error(
        make_constraint_test_model(root_max_age = 0.9),
        "minimum permissible root age"
    )
})


test_that("a bounded root transform always returns ages inside its interval", {
    model <- make_constraint_test_model(root_max_age = 20)
    root <- model$root
    root_min <- model$minimum_ages[root]

    for (root_theta in c(-5, 0, 5)) {
        theta <- rep(0, model$n_node)
        theta[1] <- root_theta
        ages <- treesurgeon:::theta_to_node_ages(theta, model)

        expect_gt(ages[root], root_min)
        expect_lt(ages[root], model$root_max_age)
        expect_true(all(ages[model$parent] > ages[model$child]))
    }
})


test_that("node-age and parameter transforms round-trip with a root maximum", {
    model <- make_constraint_test_model(root_max_age = 20)
    ages0 <- treesurgeon:::initial_node_ages(model)
    theta <- treesurgeon:::node_ages_to_theta(ages0, model)
    ages1 <- treesurgeon:::theta_to_node_ages(theta, model)

    expect_equal(ages1, ages0, tolerance = 1e-10)
})


test_that("the variance floor is inactive when the empirical variance exceeds it", {
    state <- list(log_rates = c(-1, 0, 1))

    unconstrained <- treesurgeon:::estimate_lognormal_parameters(state)
    constrained <- treesurgeon:::estimate_lognormal_parameters(
        state,
        sigma_min = 0.1
    )

    expect_equal(constrained$sigma2_mle, 2 / 3, tolerance = 1e-12)
    expect_equal(constrained$sigma2, unconstrained$sigma2, tolerance = 1e-12)
    expect_false(constrained$floor_active)

    expect_equal(
        treesurgeon:::logLik_lognormal_profile(state, sigma_min = 0.1),
        treesurgeon:::logLik_lognormal_profile(state),
        tolerance = 1e-12
    )
})


test_that("the variance floor uses the full normal likelihood when active", {
    log_rates <- c(-0.01, 0, 0.01)
    state <- list(log_rates = log_rates)
    sigma_min <- 0.1

    pars <- treesurgeon:::estimate_lognormal_parameters(
        state,
        sigma_min = sigma_min
    )

    expected_sigma2_mle <- mean((log_rates - mean(log_rates))^2)
    expected_sigma2 <- sigma_min^2
    expected_loglik <-
        -0.5 * length(log_rates) * log(2 * pi * expected_sigma2) -
        sum((log_rates - mean(log_rates))^2) / (2 * expected_sigma2)

    expect_equal(pars$sigma2_mle, expected_sigma2_mle, tolerance = 1e-12)
    expect_equal(pars$sigma2, expected_sigma2, tolerance = 1e-12)
    expect_true(pars$floor_active)
    expect_equal(
        treesurgeon:::logLik_lognormal_profile(state, sigma_min = sigma_min),
        expected_loglik,
        tolerance = 1e-12
    )
})


test_that("a fitted model respects and reports both constraints", {
    sigma_min <- 0.05
    # Keep the bound well above the expected optimum in this integration test
    # so it verifies the constrained fit without relying on a boundary warning.
    root_max_age <- 1000
    model <- make_constraint_test_model(
        sigma_min = sigma_min,
        root_max_age = root_max_age
    )

    fit <- treesurgeon::fit_timescale_ml(
        model,
        method = "Nelder-Mead",
        control = list(maxit = 3000, reltol = 1e-8)
    )

    root_age <- fit$node_ages[model$root]
    root_range <- root_max_age - model$minimum_ages[model$root]

    expect_s3_class(fit, "timescale_ml_fit")
    expect_true(all(is.finite(fit$node_ages)))
    expect_true(is.finite(fit$logLik))
    expect_gt(root_age, model$minimum_ages[model$root])
    expect_true(root_age <= root_max_age + 1e-10)

    expect_equal(fit$sigma_min, sigma_min)
    expect_equal(fit$root_max_age, root_max_age)
    expect_equal(
        fit$sigma2,
        max(fit$sigma2_mle, sigma_min^2),
        tolerance = 1e-10
    )
    expect_identical(
        fit$sigma_floor_active,
        fit$sigma2_mle < sigma_min^2
    )
    expect_identical(
        fit$root_max_active,
        (root_max_age - root_age) <= 1e-4 * root_range
    )
})


test_that("upper-bound diagnostics distinguish boundary and interior ages", {
    lower <- 1
    upper <- 10

    expect_true(treesurgeon:::is_near_upper_bound(9.9995, lower, upper))
    expect_true(treesurgeon:::is_near_upper_bound(10, lower, upper))
    expect_false(treesurgeon:::is_near_upper_bound(9.9, lower, upper))
    expect_false(treesurgeon:::is_near_upper_bound(10.1, lower, upper))
})


test_that("a fit with a tight root upper bound stays feasible and reports its status", {
    model <- make_constraint_test_model(
        sigma_min = 0.05,
        root_max_age = 1.001
    )

    warning_messages <- character()
    fit <- withCallingHandlers(
        treesurgeon::fit_timescale_ml(
            model,
            method = "Nelder-Mead",
            control = list(maxit = 5000, reltol = 1e-10)
        ),
        warning = function(w) {
            warning_messages <<- c(warning_messages, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    root_age <- fit$node_ages[model$root]
    root_min <- model$minimum_ages[model$root]

    expect_true(is.finite(fit$logLik))
    expect_gt(root_age, root_min)
    expect_lte(root_age, model$root_max_age + 1e-10)
    expect_identical(
        fit$root_upper_active,
        treesurgeon:::is_near_upper_bound(
            root_age, root_min, model$root_max_age
        )
    )
    expect_identical(
        any(grepl("hard upper bound", warning_messages, fixed = TRUE)),
        fit$root_upper_active
    )
})

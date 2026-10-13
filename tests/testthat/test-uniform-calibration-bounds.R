test_that("uniform calibration maxima constrain transformed node ages", {
    model <- make_constraint_test_model(root_max_age = 1000)

    calibrated_node <- which(vapply(
        model$calibrations,
        function(calibration) {
            !is.null(calibration) && identical(calibration$distribution, "uniform")
        },
        logical(1)
    ))
    expect_length(calibrated_node, 1L)
    calibrated_node <- calibrated_node[[1]]
    calibration <- model$calibrations[[calibrated_node]]
    parameter_index <- match(calibrated_node, model$internal_nodes_preorder)

    for (root_theta in c(-5, 0, 5)) {
        for (node_theta in c(-5, 0, 5)) {
            theta <- rep(0, model$n_node)
            theta[1] <- root_theta
            theta[parameter_index] <- node_theta

            ages <- treesurgeon:::theta_to_node_ages(theta, model)
            expect_gt(ages[calibrated_node], calibration$age_min)
            expect_lt(ages[calibrated_node], calibration$age_max)
            expect_true(all(ages[model$parent] > ages[model$child]))
            expect_true(is.finite(treesurgeon:::calibration_log_density(
                ages[calibrated_node], calibration
            )))
        }
    }
})


test_that("initial ages respect uniform maxima even with a very old root bound", {
    model <- make_constraint_test_model(root_max_age = 1000)
    ages <- treesurgeon:::initial_node_ages(model)

    calibrated_node <- which(vapply(
        model$calibrations,
        function(calibration) {
            !is.null(calibration) && identical(calibration$distribution, "uniform")
        },
        logical(1)
    ))[[1]]
    calibration <- model$calibrations[[calibrated_node]]

    expect_lt(ages[calibrated_node], calibration$age_max)
    expect_gt(ages[calibrated_node], calibration$age_min)

    theta <- treesurgeon:::node_ages_to_theta(ages, model)
    expect_true(is.finite(treesurgeon:::objective(theta, model)))

    recovered <- treesurgeon:::theta_to_node_ages(theta, model)
    expect_equal(recovered, ages, tolerance = 1e-10)
})


test_that("a uniform calibration maximum also bounds the root when no root_max_age is set", {
    tree <- make_constraint_test_tree()
    root <- setdiff(tree$edge[, 1], tree$edge[, 2])[[1]]
    calibrations <- vector("list", ape::Ntip(tree) + ape::Nnode(tree))
    calibrations[[root]] <- list(
        distribution = "uniform",
        age_min = 1,
        age_max = 10
    )

    model <- treesurgeon::timescale_ml(
        tree = tree,
        tip_ages = rep(0, ape::Ntip(tree)),
        calibrations = calibrations,
        root_max_age = NULL
    )

    theta <- rep(0, model$n_node)
    theta[1] <- 5
    ages <- treesurgeon:::theta_to_node_ages(theta, model)

    expect_lt(ages[root], 10)
    expect_gt(ages[root], model$minimum_ages[root])
})


test_that("gamma and exponential age_max values remain soft maxima", {
    tree <- make_constraint_test_tree()
    root <- setdiff(tree$edge[, 1], tree$edge[, 2])[[1]]
    non_root_internal <- setdiff(unique(tree$edge[, 1]), root)[[1]]

    soft_calibrations <- list(
        gamma = list(
            distribution = "gamma", age_min = 1, age_max = 10,
            shape = 2, rate = 0.5
        ),
        exponential = list(
            distribution = "exponential", age_min = 1, age_max = 10,
            rate = 0.2
        )
    )

    for (calibration in soft_calibrations) {
        calibrations <- vector("list", ape::Ntip(tree) + ape::Nnode(tree))
        calibrations[[non_root_internal]] <- calibration
        model <- treesurgeon::timescale_ml(
            tree = tree,
            tip_ages = rep(0, ape::Ntip(tree)),
            calibrations = calibrations,
            root_max_age = 20
        )

        theta <- rep(0, model$n_node)
        theta[1] <- 0
        node_parameter_index <- match(
            non_root_internal,
            model$internal_nodes_preorder
        )
        theta[node_parameter_index] <- 5
        ages <- treesurgeon:::theta_to_node_ages(theta, model)

        expect_gt(ages[non_root_internal], calibration$age_max)
        expect_true(is.finite(treesurgeon:::calibration_log_density(
            ages[non_root_internal], calibration
        )))
    }
})


test_that("infeasible uniform upper bounds are rejected during construction", {
    tree <- make_constraint_test_tree()
    calibrations <- make_constraint_test_calibrations(tree)
    calibrated_node <- which(vapply(
        calibrations,
        function(calibration) {
            !is.null(calibration) && identical(calibration$distribution, "uniform")
        },
        logical(1)
    ))[[1]]
    calibrations[[calibrated_node]]$age_max <- 1

    expect_error(
        treesurgeon::timescale_ml(
            tree = tree,
            tip_ages = rep(0, ape::Ntip(tree)),
            calibrations = calibrations
        ),
        "Infeasible uniform calibration"
    )
})


test_that("root_max_age and a uniform root maximum use the tighter bound", {
    tree <- make_constraint_test_tree()
    root <- setdiff(tree$edge[, 1], tree$edge[, 2])[[1]]
    calibrations <- vector("list", ape::Ntip(tree) + ape::Nnode(tree))
    calibrations[[root]] <- list(
        distribution = "uniform",
        age_min = 1,
        age_max = 10
    )

    model <- treesurgeon::timescale_ml(
        tree = tree,
        tip_ages = rep(0, ape::Ntip(tree)),
        calibrations = calibrations,
        root_max_age = 20
    )

    theta <- rep(0, model$n_node)
    theta[1] <- 20
    ages <- treesurgeon:::theta_to_node_ages(theta, model)
    expect_lt(ages[root], 10)
    expect_equal(treesurgeon:::node_upper_age(model, root), 10)
})

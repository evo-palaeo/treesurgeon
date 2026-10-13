# Small deterministic fixture used by the constraint tests.

make_constraint_test_tree <- function() {
    ape::read.tree(text = "((A:1,B:1):1,C:2);")
}

make_constraint_test_calibrations <- function(tree = make_constraint_test_tree()) {
    n_nodes_total <- ape::Ntip(tree) + ape::Nnode(tree)
    calibrations <- vector("list", n_nodes_total)

    internal_nodes <- unique(tree$edge[, 1])
    root <- setdiff(tree$edge[, 1], tree$edge[, 2])
    non_root_internal <- setdiff(internal_nodes, root)

    stopifnot(length(root) == 1L, length(non_root_internal) >= 1L)

    # A uniform calibration gives the non-root internal node a hard minimum
    # age of 1, making the minimum permissible root age easy to reason about.
    calibrated_node <- non_root_internal[[1]]
    calibrations[[calibrated_node]] <- list(
        distribution = "uniform",
        age_min = 1,
        age_max = 10
    )

    calibrations
}

make_constraint_test_model <- function(sigma_min = 0.01,
                                       root_max_age = NULL,
                                       tree = make_constraint_test_tree()) {
    treesurgeon::timescale_ml(
        tree = tree,
        tip_ages = rep(0, ape::Ntip(tree)),
        calibrations = make_constraint_test_calibrations(tree),
        sigma_min = sigma_min,
        root_max_age = root_max_age
    )
}

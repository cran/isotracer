### * Setup

library(magrittr)

ITERS <- 20
DRAWS <- 10
FORCE_STAN_COVR <- TRUE

n_cores <- min(2, parallel::detectCores())
n_chains <- max(n_cores, 2)

run_mcmc <- function(...) {
  isotracer:::run_mcmc(..., cores = n_cores, chains = n_chains)
}

new_networkModel <- function() {
    isotracer::new_networkModel(quiet = TRUE)
}

### * Two compartments

test_that("Model with two compartments works", {
    m <- new_networkModel() %>%
        set_topo("A -> B") %>%
        set_init(tibble::tibble(comp = c("A", "B"),
                                size = c(100, 100),
                                prop = c(0.5, 0)),
                 comp = "comp", size = "size", prop = "prop")
    expect_setequal(params(m), c("eta", "lambda_A", "lambda_B",
                                 "upsilon_A_to_B", "zeta"))
    # Run model with set parameters
    m %<>% set_params(c("upsilon_A_to_B" = 0.1,
                        "eta" = 0.05, "zeta" = 0.05,
                        "lambda_A" = 0, "lambda_B" = 0)) %>%
        project(end = 50)
    traj <- m$trajectory[[1]]
    last <- length(traj$timepoints[[1]])
    expect_equal(traj$sizes[[1]][1, ], c("A" = 100, "B" = 100))
    x <- traj$sizes[[1]][last, "A"]
    expect_true(0 < x & x < 1)
    x <- traj$sizes[[1]][last, "B"]
    expect_true(199 < x & x < 200)
    expect_equal(traj$proportions[[1]][1, ], c("A" = 0.5, "B" = 0))
    x <- traj$proportions[[1]][last, "A"]
    expect_true(x == 0.5)
    x <- traj$proportions[[1]][last, "B"]
    expect_true(0.249 < x & x < 0.251)
    expect_error(plot(m), NA)
    # Run model fit
    obs <- m %>% sample_from(at = seq(1, 50, length.out = 5))
    m %<>% set_obs(obs, comp = "comp", size = "size", prop = "prop",
                   time = "time")
    plot(m)
    if (!covr::in_covr() | FORCE_STAN_COVR) {
        expect_error(
            capture_warnings(capture_output({f <- run_mcmc(m, iter = ITERS)})),
            NA)
        expect_error(plot(f), NA)
        # Posterior predictions
        expect_error({p <- predict(m, f, draws = DRAWS, cores = n_cores)}, NA)
        expect_error(plot(p), NA)
    }
})

### * Two compartments, one split

test_that("Model with two compartments, one split works", {
    m <- new_networkModel() %>%
        set_topo("A -> B") %>%
        set_split("B") %>%
        set_init(tibble::tibble(comp = c("A", "B"),
                                size = c(100, 100),
                                prop = c(0.5, 0)),
                 comp = "comp", size = "size", prop = "prop")
    expect_setequal(params(m), c("eta", "lambda_A", "lambda_B", "portion.act_B",
                                 "upsilon_A_to_B", "zeta"))
    # Run model with set parameters
    m %<>% set_params(c("upsilon_A_to_B" = 0.1,
                        "eta" = 0.01, "zeta" = 0.01,
                        "lambda_A" = 0, "lambda_B" = 0,
                        "portion.act_B" = 0.25)) %>%
        project(end = 50)
    traj <- m$trajectory[[1]]
    last <- length(traj$timepoints[[1]])
    expect_equal(traj$sizes[[1]][1, ], c("A" = 100, "B" = 100))
    x <- traj$sizes[[1]][last, "A"]
    expect_true(0 < x & x < 1)
    x <- traj$sizes[[1]][last, "B"]
    expect_true(199 < x & x < 200)
    expect_equal(traj$proportions[[1]][1, ], c("A" = 0.5, "B" = 0))
    x <- traj$proportions[[1]][last, "A"]
    expect_true(x == 0.5)
    x <- traj$proportions[[1]][last, "B"]
    expect_true(0.249 < x & x < 0.251)
    expect_error(plot(m), NA)
    # Run model fit
    obs <- m %>% sample_from(at = seq(1, 50, length.out = 5))
    m %<>% set_obs(obs, comp = "comp", size = "size", prop = "prop",
                   time = "time")
    plot(m)
    if (!covr::in_covr() | FORCE_STAN_COVR) {
        expect_error(
            capture_warnings(capture_output({f <- run_mcmc(m, iter = ITERS)})),
            NA)
        expect_error(plot(f), NA)
        # Posterior predictions
        expect_error({p <- predict(m, f, draws = DRAWS, cores = n_cores)}, NA)
        expect_error(plot(p), NA)
    }
})

### * Three compartments

test_that("Model with three compartments work", {
    m <- new_networkModel() %>%
        set_topo("A -> B -> C") %>%
        set_init(tibble::tibble(comp = c("A", "B", "C"),
                                size = c(100, 100, 100),
                                prop = c(0.5, 0, 0)),
                 comp = "comp", size = "size", prop = "prop")
    expect_setequal(params(m), c("eta", "lambda_A", "lambda_B", "lambda_C",
                                 "upsilon_A_to_B", "upsilon_B_to_C", "zeta"))
    # Run model with set parameters
    m %<>% set_params(c("upsilon_A_to_B" = 0.1,
                        "upsilon_B_to_C" = 0.1,
                        "eta" = 0.01, "zeta" = 0.01,
                        "lambda_A" = 0, "lambda_B" = 0,
                        "lambda_C" = 0)) %>%
        project(end = 50)
    traj <- m$trajectory[[1]]
    last <- length(traj$timepoints[[1]])
    expect_equal(traj$sizes[[1]][1, ], c("A" = 100, "B" = 100, "C" = 100))
    x <- traj$sizes[[1]][last, "A"]
    expect_true(0 < x & x < 1)
    x <- traj$sizes[[1]][last, "B"]
    expect_true(3 < x & x < 4)
    x <- traj$sizes[[1]][last, "C"]
    expect_true(295 < x & x < 296)
    expect_equal(traj$proportions[[1]][1, ], c("A" = 0.5, "B" = 0, "C" = 0))
    x <- traj$proportions[[1]][last, "A"]
    expect_true(x == 0.5)
    x <- traj$proportions[[1]][last, "B"]
    expect_true(0.4 < x & x < 0.42)
    x <- traj$proportions[[1]][last, "C"]
    expect_true(0.16 < x & x < 0.17)
    expect_error(plot(m), NA)
    # Run model fit
    obs <- m %>% sample_from(at = seq(1, 50, length.out = 5))
    m %<>% set_obs(obs, comp = "comp", size = "size", prop = "prop",
                   time = "time")
    plot(m)
    if (!covr::in_covr() | FORCE_STAN_COVR) {
        expect_error(
            capture_warnings(capture_output({f <- run_mcmc(m, iter = ITERS)})),
            NA)
        expect_error(plot(f), NA)
        # Posterior predictions
        expect_error({p <- predict(m, f, draws = DRAWS, cores = n_cores)}, NA)
        expect_error(plot(p), NA)
    }
})

### * Three compartments, one split

test_that("Model with three compartments, one split works", {
    m <- new_networkModel() %>%
        set_topo("A -> B -> C") %>%
        set_split("B") %>%
        set_init(tibble::tibble(comp = c("A", "B", "C"),
                                size = c(100, 100, 100),
                                prop = c(0.5, 0, 0)),
                 comp = "comp", size = "size", prop = "prop")
    expect_setequal(params(m), c("eta", "lambda_A", "lambda_B", "lambda_C",
                                 "portion.act_B", "upsilon_A_to_B", "upsilon_B_to_C",
                                 "zeta"))
    # Run model with set parameters
    m %<>% set_params(c("upsilon_A_to_B" = 0.1,
                        "upsilon_B_to_C" = 0.1,
                        "eta" = 0.01, "zeta" = 0.01,
                        "lambda_A" = 0, "lambda_B" = 0,
                        "lambda_C" = 0,
                        "portion.act_B" = 0.25)) %>%
        project(end = 50)
    traj <- m$trajectory[[1]]
    last <- length(traj$timepoints[[1]])
    expect_equal(traj$sizes[[1]][1, ], c("A" = 100, "B" = 100, "C" = 100))
    x <- traj$sizes[[1]][last, "A"]
    expect_true(0 < x & x < 1)
    x <- traj$sizes[[1]][last, "B"]
    expect_true(78 < x & x < 79)
    x <- traj$sizes[[1]][last, "C"]
    expect_true(220 < x & x < 221)
    expect_equal(traj$proportions[[1]][1, ], c("A" = 0.5, "B" = 0, "C" = 0))
    x <- traj$proportions[[1]][last, "A"]
    expect_true(x == 0.5)
    x <- traj$proportions[[1]][last, "B"]
    expect_true(0.02 < x & x < 0.03)
    x <- traj$proportions[[1]][last, "C"]
    expect_true(0.21 < x & x < 0.22)
    expect_error(plot(m), NA)
    # Run model fit
    obs <- m %>% sample_from(at = seq(1, 50, length.out = 5))
    m %<>% set_obs(obs, comp = "comp", size = "size", prop = "prop",
                   time = "time")
    plot(m)
    if (!covr::in_covr() | FORCE_STAN_COVR) {
        expect_error(
            capture_warnings(capture_output({f <- run_mcmc(m, iter = ITERS)})),
            NA)
        expect_error(plot(f), NA)
        # Posterior predictions
        expect_error({p <- predict(m, f, draws = DRAWS, cores = n_cores)}, NA)
        expect_error(plot(p), NA)
    }
})

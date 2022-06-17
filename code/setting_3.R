library(tidyverse)
library(tdnn)
library(kknn)
library(argparser)
library(kableExtra)

ROUNDING_DIGITS <- 4
OUTPUT_DIR <- "../output"

# Simulation Parameters ---------------------------------------------------

# Parse arguments provided to Rscript command
parser <- arg_parser("Setting 3 params")
parser <-
  add_argument(
    parser,
    "--dimension",
    type = "integer",
    default = 3,
    help = "Dimension for the setting [default 3]"
  )


opt <- parse_args(parser)

p <- opt$dimension
data_type <- "unif"
draw_random_data <- "fixed"
est_variance <- TRUE
verbose <- TRUE

# Sets the number of threads to be used for the Rcpp code in the tdnn package
num_threads <- 1


seed_val <- 1234
num_reps <- 1000
n_test <- 100
n <- 1000

c_seq <- c(2, 4, 6, 8, 10, 15, 20, 25, 30)
c_seq_fixed <- 2
s_1_seq <- seq(1, 100, 1)
dnn_s_seq <- seq(10, 100, 1)
alpha <- 0.001




grid_start_end <- c(0, 1)

fixed_test_vector <- c(c(0.2, 0.4, 0.6), rep(0.5, p - 3))

dgp_function <- function(x) {
  prod(sapply(c(1, 2, 3), function(i) {
    (1 + 1 / (1 + exp(-20 * (x[i] - 1 / 3))))
  }))
}

baseline <- function(x) {
  0.125 * (x[1] - 1)
}

# Simulation helper functions ---------------------------------------------
make_results_df <- function(truth, predictions, variance, method) {
  result_df <- data.frame(cbind(truth, predictions, variance)) %>%
    mutate(method = method)
  colnames(result_df) <-
    c("truth", "predictions", "variance", "method")
  result_df
}

knn_reg <- function(X,
                    Y,
                    X_test,
                    k_val = 20,
                    est_variance = T) {
  knn_mu <- FNN::knn.reg(X, X_test, Y, k = k_val)$pred
  if (est_variance) {
    knn_mu2 <- FNN::knn.reg(X, X_test, Y^2, k = k_val)$pred
    knn_var <- (knn_mu2 - knn_mu^2) / (k_val - 1)
  } else {
    knn_var <- NA
  }
  list(
    pred_knn = knn_mu,
    variance_knn = knn_var
  )
}


make_tuned_knn_results <-
  function(X, Y, X_test, W_0, kmax, estimate_variance = T) {
    X <- X[, as.logical(W_0)]
    X_test <- X_test[, as.logical(W_0)]
    training_data <- cbind(data.frame(X), data.frame(Y = Y))
    testing_data <- data.frame(X_test)

    tuned_knn <- train.kknn(Y ~ ., data = training_data, kmax = 100)
    tuned_k <- tuned_knn$best.parameters$k
    tuned_results <-
      knn_reg(X, as.numeric(Y), X_test, k_val = tuned_k, estimate_variance)
    tuned_results[["best_k"]] <- tuned_k
    return(tuned_results)
  }

tune_knn_trt_effect <-
  function(X_trt,
           X_ctl,
           Y_trt,
           Y_ctl,
           X_test,
           W_0,
           kmax,
           estimate_variance = T) {
    tuned_trt_knn <-
      make_tuned_knn_results(X_trt, Y_trt, X_test, W_0, kmax, estimate_variance)
    tuned_ctl_knn <-
      make_tuned_knn_results(X_ctl, Y_ctl, X_test, W_0, kmax, estimate_variance)

    trt_effect <- tuned_trt_knn$pred_knn - tuned_ctl_knn$pred_knn
    variance <- tuned_trt_knn$variance_knn + tuned_ctl_knn$variance_knn
    list(
      trt_effect = trt_effect,
      variance = variance,
      best_k_trt = tuned_trt_knn$best_k,
      best_k_ctl = tuned_ctl_knn$best_k
    )
  }

tune_dnn_trt_effect <-
  function(X_trt,
           X_ctl,
           Y_trt,
           Y_ctl,
           X_test,
           W_0,
           s_seq,
           num_threads,
           estimate_variance = T) {
    trt_est <- tdnn:::tune_dnn_no_dist_thread(
      X_trt,
      Y_trt,
      X_test,
      s_seq = s_seq,
      W0_ = W_0,
      estimate_variance = estimate_variance,
      num_threads = num_threads
    )

    ctl_est <- tdnn:::tune_dnn_no_dist_thread(
      X_ctl,
      Y_ctl,
      X_test,
      s_seq = s_seq,
      W0_ = W_0,
      estimate_variance = estimate_variance,
      num_threads = num_threads
    )
    trt_effect <- trt_est$estimate_loo - ctl_est$estimate_loo
    variance <- as.numeric(trt_est$variance) + as.numeric(ctl_est$variance)
    list(
      trt_effect = trt_effect,
      variance = variance,
      s_1_trt = trt_est$s_1_B_NN,
      s_1_ctl = ctl_est$s_1_B_NN
    )
  }

tune_tdnn_trt_effect <-
  function(X_trt,
           X_ctl,
           Y_trt,
           Y_ctl,
           X_test,
           W_0,
           c_seq,
           num_threads,
           estimate_variance = T) {
    trt_est <-
      tdnn:::tune_de_dnn_no_dist_vary_c_cpp_thread(
        X_trt,
        Y_trt,
        X_test,
        c = c_seq,
        W0_ = W_0,
        estimate_variance = estimate_variance,
        num_threads = num_threads,
        bootstrap_iter = 500
      )

    ctl_est <-
      tdnn:::tune_de_dnn_no_dist_vary_c_cpp_thread(
        X_ctl,
        Y_ctl,
        X_test,
        c = c_seq,
        W0_ = W_0,
        estimate_variance = estimate_variance,
        num_threads = num_threads,
        bootstrap_iter = 500
      )
    trt_effect <- trt_est$estimate_loo - ctl_est$estimate_loo
    variance <- as.numeric(trt_est$variance) + as.numeric(ctl_est$variance)
    list(
      trt_effect = trt_effect,
      variance = variance,
      s_1_trt = trt_est$s_1_B_NN,
      c_trt = trt_est$c_B_NN,
      s_1_ctl = ctl_est$s_1_B_NN,
      c_ctl = ctl_est$c_B_NN
    )
  }

# Draw random test data ---------------------------------------------------
set.seed(seed_val)
X_test_random <- matrix(runif(n_test * p, 0, 1), n_test, p)

mu_rand <-
  apply(X_test_random, MARGIN = 1, baseline) + apply(X_test_random, MARGIN = 1, dgp_function)
mu_fixed <-
  baseline(fixed_test_vector) + dgp_function(fixed_test_vector)

# Main simulation loop ----------------------------------------------------
run_sim <-
  function(i,
           n,
           c_seq,
           X_test_random,
           mu_rand,
           mu_fixed,
           draw_random_data,
           est_variance = F,
           verbose = F) {
    if (verbose) print(glue::glue("rep {i}/{num_reps}, data type: {data_type}"))

    X <- matrix(runif(n * p, 0, 1), n, p)
    X_test_fixed <- matrix(fixed_test_vector, 1, p)
    W <- rbinom(n, 1, 0.5)
    epsi <- matrix(rnorm(n), n, 1)

    Y <- apply(X, 1, baseline) + (W - 0.5) * apply(X, MARGIN = 1, dgp_function) + epsi

    X_trt <- X[W == 1, ]
    Y_trt <- Y[W == 1]
    X_ctl <- X[W == 0, ]
    Y_ctl <- Y[W == 0]

    if (p > 3) {
      W_0 <- tdnn:::feature_screen(X, Y, alpha = alpha)
    } else {
      W_0 <- rep(1, p)
    }
    if (verbose) print(W_0)

    # TDNN Estimates ----------------------------------------------------
    tdnn_rand_estimate <-
      tdnn:::tune_treatment_effect_thread(
        X,
        Y,
        W,
        X_test_random,
        W0_ = W_0,
        c = c_seq,
        estimate_variance = est_variance,
        bootstrap_iter = 500,
        num_threads = num_threads
      )
    tdnn_fixed_estimate_R <-
      tune_tdnn_trt_effect(
        X_trt,
        X_ctl,
        Y_trt,
        Y_ctl,
        X_test_fixed,
        W_0,
        c_seq_fixed,
        num_threads = num_threads,
        estimate_variance = est_variance
      )
    # Collect TDNN estimates ----------------------------------------------------------
    tdnn_fixed_results_r <-
      data.frame(
        estimate = tdnn_fixed_estimate_R$trt_effect,
        loss = (tdnn_fixed_estimate_R$trt_effect - mu_fixed)^
          2,
        bias = (tdnn_fixed_estimate_R$trt_effect - mu_fixed),
        s_1_trt = tdnn_fixed_estimate_R$s_1_trt,
        s_1_ctl = tdnn_fixed_estimate_R$s_1_ctl,
        variance = if (est_variance) {
          as.numeric(tdnn_fixed_estimate_R$variance)
        } else {
          NA
        },
        method = "tdnn_fixed"
      ) %>% mutate(
        covered = case_when(
          str_ends(method, "fixed") ~ abs(estimate - mu_fixed) <= qnorm(1 - 0.05 /
            2) * sqrt(variance)
        ),
        width = 2 * qnorm(1 - 0.05 / 2) * sqrt(variance)
      )

    tdnn_rand_results <-
      data.frame(
        estimate = tdnn_rand_estimate$treatment_effect,
        loss = (tdnn_rand_estimate$treatment_effect - mu_rand)^
          2,
        bias = (tdnn_rand_estimate$treatment_effect - mu_rand),
        s_1_ctl = tdnn_rand_estimate$s_1_B_NN_ctl,
        c_ctl = tdnn_rand_estimate$c_B_NN_ctl,
        s_1_trt = tdnn_rand_estimate$s_1_B_NN_trt,
        c_trt = tdnn_rand_estimate$c_B_NN_trt,
        variance = if (est_variance) {
          as.numeric(tdnn_rand_estimate$variance)
        } else {
          NA
        },
        method = "tdnn_rand"
      ) %>% mutate(
        covered = case_when(
          str_ends(method, "rand") ~ abs(estimate - mu_rand) <= qnorm(1 - 0.05 / 2) * sqrt(variance)
        ),
        width = 2 * qnorm(1 - 0.05 / 2) * sqrt(variance)
      )

    # DNN Estimates ----------------------------------------------------
    dnn_fixed_estimate <-
      tune_dnn_trt_effect(
        X_trt,
        X_ctl,
        Y_trt,
        Y_ctl,
        X_test_fixed,
        W_0,
        dnn_s_seq,
        num_threads = num_threads,
        estimate_variance = est_variance
      )
    dnn_rand_estimate <-
      tune_dnn_trt_effect(
        X_trt,
        X_ctl,
        Y_trt,
        Y_ctl,
        X_test_random,
        W_0,
        dnn_s_seq,
        num_threads = num_threads,
        estimate_variance = est_variance
      )
    # Collect DNN estimates ----------------------------------------------------------
    dnn_rand_results <-
      data.frame(
        estimate = dnn_rand_estimate$trt_effect,
        loss = (dnn_rand_estimate$trt_effect - mu_rand)^
          2,
        bias = (dnn_rand_estimate$trt_effect - mu_rand),
        s_1_trt = dnn_rand_estimate$s_1_trt,
        s_1_ctl = dnn_rand_estimate$s_1_ctl,
        variance = if (est_variance) {
          as.numeric(dnn_rand_estimate$variance)
        } else {
          NA
        },
        method = "dnn_rand"
      ) %>% mutate(
        covered = case_when(
          str_ends(method, "rand") ~ abs(estimate - mu_rand) <= qnorm(1 - 0.05 / 2) * sqrt(variance)
        ),
        width = 2 * qnorm(1 - 0.05 / 2) * sqrt(variance)
      )
    dnn_fixed_results <-
      data.frame(
        estimate = dnn_fixed_estimate$trt_effect,
        loss = (dnn_fixed_estimate$trt_effect - mu_fixed)^
          2,
        bias = (dnn_fixed_estimate$trt_effect - mu_fixed),
        s_1_trt = dnn_fixed_estimate$s_1_trt,
        s_1_ctl = dnn_fixed_estimate$s_1_ctl,
        variance = if (est_variance) {
          as.numeric(dnn_fixed_estimate$variance)
        } else {
          NA
        },
        method = "dnn_fixed"
      ) %>% mutate(
        covered = case_when(
          str_ends(method, "fixed") ~ abs(estimate - mu_fixed) <= qnorm(1 - 0.05 /
            2) * sqrt(variance)
        ),
        width = 2 * qnorm(1 - 0.05 / 2) * sqrt(variance)
      )
    # KNN Estimates ----------------------------------------------------
    knn_fixed_estimate <-
      tune_knn_trt_effect(
        X_trt,
        X_ctl,
        Y_trt,
        Y_ctl,
        X_test_fixed,
        W_0,
        kmax = 200,
        estimate_variance = est_variance
      )
    knn_rand_estimate <-
      tune_knn_trt_effect(
        X_trt,
        X_ctl,
        Y_trt,
        Y_ctl,
        X_test_random,
        W_0,
        kmax = 200,
        estimate_variance = est_variance
      )
    # Collect KNN estimates ----------------------------------------------------------
    knn_rand_results <-
      data.frame(
        estimate = knn_rand_estimate$trt_effect,
        loss = (knn_rand_estimate$trt_effect - mu_rand)^
          2,
        bias = (knn_rand_estimate$trt_effect - mu_rand),
        k_trt = knn_rand_estimate$best_k_trt,
        k_ctl = knn_rand_estimate$best_k_ctl,
        variance = if (est_variance) {
          as.numeric(knn_rand_estimate$variance)
        } else {
          NA
        },
        method = "knn_rand"
      ) %>% mutate(
        covered = case_when(
          str_ends(method, "rand") ~ abs(estimate - mu_rand) <= qnorm(1 - 0.05 / 2) * sqrt(variance)
        ),
        width = 2 * qnorm(1 - 0.05 / 2) * sqrt(variance)
      )
    knn_fixed_results <-
      data.frame(
        estimate = knn_fixed_estimate$trt_effect,
        loss = (knn_fixed_estimate$trt_effect - mu_fixed)^
          2,
        bias = (knn_fixed_estimate$trt_effect - mu_fixed),
        k_trt = knn_fixed_estimate$best_k_trt,
        k_ctl = knn_fixed_estimate$best_k_ctl,
        variance = if (est_variance) {
          as.numeric(knn_fixed_estimate$variance)
        } else {
          NA
        },
        method = "knn_fixed"
      ) %>% mutate(
        covered = case_when(
          str_ends(method, "fixed") ~ abs(estimate - mu_fixed) <= qnorm(1 - 0.05 /
            2) * sqrt(variance)
        ),
        width = 2 * qnorm(1 - 0.05 / 2) * sqrt(variance)
      )

    # Combine results ----------------------------------------------------
    results <- bind_rows(
      tdnn_rand_results,
      tdnn_fixed_results_r,
      dnn_rand_results,
      dnn_fixed_results,
      knn_rand_results,
      knn_fixed_results
    )

    return(results)
  }

# Run simulation loop ----------------------------------------------------------
set.seed(seed_val)
results <- map_dfr(1:num_reps, function(i) {
  run_sim(
    i,
    n,
    c_seq,
    X_test_random,
    mu_rand,
    mu_fixed,
    draw_random_data,
    est_variance = est_variance,
    verbose = verbose
  )
})



# View and Save Results ---------------------------------------------------
results %>%
  group_by(method) %>%
  summarize(
    MSE = mean(loss),
    Bias = mean(bias),
    Variance = mean(variance),
    Coverage = mean(covered),
    Width = mean(width)
  ) %>%
  arrange(MSE) %>%
  print()


results_data <- list(
  results_df = results,
  c_seq = c_seq,
  mu_rand = mu_rand,
  mu_fixed = mu_fixed,
  p = p,
  n = n,
  num_reps = num_reps,
  n_test = n_test,
  data_type = data_type,
  fixed_test_vector = fixed_test_vector,
  draw_random_data = draw_random_data
)
file_path <-
  file.path(OUTPUT_DIR, glue::glue(
    "tdnn_setting_3_p_{p}_{data_type}_data_{draw_random_data}_{str_replace_all(strftime(Sys.time(), '%Y-%m-%d_%H_%M'),'-', '_')}.rds"
  ))
results_data %>%
  write_rds(., file_path)

# Make LaTeX results table ------------------------------------------------
hlinesep <- function(x, y = character()) {
  if (!length(x)) {
    return(y)
  }
  hlinesep(x[-length(x)], c(rep("", x[length(x)] - 1), "\\hline", y))
}


setting_3_results_df <- results %>%
  mutate(
    method = case_when(
      method == "tdnn_fixed" ~ "TDNN fixed vector",
      method == "tdnn_rand" ~ "TDNN random vector",
      method == "dnn_fixed" ~ "DNN fixed vector",
      method == "dnn_rand" ~ "DNN random vector",
      method == "knn_fixed" ~ "KNN fixed vector",
      method == "knn_rand" ~ "KNN random vector"
    )
  ) %>%
  rename(Method = method) %>%
  mutate(p = p, point = rep(1:(n_test + 1), times = 3 * num_reps))

fixed_results <- setting_3_results_df %>%
  filter(str_ends(Method, "fixed vector")) %>%
  # mutate(point = rep(1:1000, each = 3)) %>%
  group_by(Method, p) %>%
  summarize(
    # MSE = mean((estimate - mu_fixed)^2),
    MSE = mean(loss),
    mean_est = mean(estimate),
    Variance = mean(variance),
    Coverage = mean(covered),
    Width = mean(width),
    .groups = "keep"
  ) %>%
  summarize(
    Method = Method, MSE = MSE,
    `Squared Bias` = (mean_est - mu_fixed)^2,
    Variance = Variance,
    Coverage = Coverage,
    Width = Width
  ) %>%
  mutate(Method = str_remove(Method, " fixed vector"))

random_results <- setting_3_res_df %>%
  ungroup() %>%
  filter(str_ends(Method, "random vector")) %>%
  group_by(Method, p, point) %>%
  summarize(
    MSE = mean(loss),
    mean_est = mean(estimate),
    Variance = mean(variance),
    Coverage = mean(covered),
    Width = mean(width)
  ) %>%
  summarize(
    MSE = MSE,
    `Squared Bias` = (mean_est - mu_rand)^2,
    Variance = Variance,
    Coverage = Coverage,
    Width = Width
  ) %>%
  summarize(across(everything(), mean), .groups = "drop") %>%
  select(-one_of(c("Method", "p")))

combined_results <- bind_cols(fixed_results, random_results)

combined_results %>%
  knitr::kable(
    digits = c(0, 0, rep(ROUNDING_DIGITS, 10)),
    "latex", booktabs = T, align = "lccccccccccc", linesep = hlinesep(c(3, 3, 3, 3))
  ) %>%
  kableExtra::add_header_above(c(" " = 2, "Fixed Test Point" = 5, "Random Test Points" = 5)) %>%
  gsub("\\.\\.\\.\\d{1,2}", "", .) %>%
  kableExtra::save_kable(., file.path(OUTPUT_DIR, glue::glue("tdnn_setting_3_combined_results_table.tex")))

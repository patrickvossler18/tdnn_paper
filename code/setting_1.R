library(tidyverse)
library(tdnn)
library(kknn)
library(FNN)
library(glue)
library(kableExtra)
library(argparser)

ROUNDING_DIGITS <- 4
OUTPUT_DIR <- "output"

# Simulation Parameters ---------------------------------------------------

# Parse arguments provided to Rscript command
parser <- arg_parser("Setting 1 params")
parser <-
    add_argument(
        parser,
        "--data_type",
        type = "character",
        default = "normal",
        help = "Choose between normal and unif for the covariate distribution [default normal]"
    )

opt <- parse_args(parser)

data_type = opt$data_type

draw_random_data <- "fixed"
est_variance <- TRUE
verbose <- TRUE

# Sets the number of threads to be used for the Rcpp code in the tdnn package
num_threads <- 1


seed_val <- 1234
num_reps <- 1000
n_test <- 100

c_seq = c(2, 4, 6, 8, 10, 15, 20, 25, 30)
c_seq_fixed <- 2
s_1_seq <- seq(1, 100, 1)
dnn_s_seq <- seq(10, 100, 1)

n_fixed <- 1000
n <- 1000
p <- 3
grid_start_end <- c(0, 1)

if (data_type == "unif") {
    fixed_test_vector <- c(0.5, 0.5, 0.5)
} else if (data_type == "normal") {
    fixed_test_vector <- c(0.5, -0.5, 0.5)
}

dgp_function <- function(x) {
    # This function will take a row of the matrix as input and return the
    # transformed value. To use this on a matrix, we can use the apply function
    # with margin=1
    (x[1] - 1) ^ 2 + (x[2] + 1) ^ 3 - 3 * x[3]
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
    knn_mu <- knn.reg(X, X_test, Y, k = k_val)$pred
    if (est_variance) {
        knn_mu2 <- knn.reg(X, X_test, Y ^ 2, k = k_val)$pred
        knn_var <- (knn_mu2 - knn_mu ^ 2) / (k_val - 1)
    } else {
        knn_var <- NA
    }
    list(pred_knn = knn_mu,
         variance_knn = knn_var)
}

make_tuned_knn_results <-
    function(X, Y, X_test, kmax, estimate_variance = T) {
        training_data <- cbind(data.frame(X), data.frame(Y = Y))
        testing_data <- data.frame(X_test)
        
        tuned_knn <- train.kknn(Y ~ ., data = training_data, kmax = 100)
        tuned_k <- tuned_knn$best.parameters$k
        tuned_results <-
            knn_reg(X, as.numeric(Y), X_test, k_val = tuned_k, estimate_variance)
        tuned_results[["best_k"]] <- tuned_k
        return(tuned_results)
    }



# Draw random test data ---------------------------------------------------
set.seed(seed_val)
# fix X_test_random
if (data_type == "unif") {
    X_test_random <- matrix(runif(n_test * p, 0, 1), n_test, p)
} else if (data_type == "normal") {
    X_test_random <- matrix(rnorm(n_test * p), n_test, p)
}

mu_rand <- apply(X_test_random, MARGIN = 1, dgp_function)
mu_fixed <- dgp_function(fixed_test_vector)

# Define main simulation loop ----------------------------------------------------
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
        if(verbose) print(glue::glue("rep {i}/{num_reps}, data type: {data_type}"))
        if (data_type == "unif") {
            X <- matrix(runif(n * p, 0, 1), n, p)
        } else if (data_type == "normal") {
            X <- matrix(rnorm(n * p), n, p)
        }
        X_test_fixed <- matrix(fixed_test_vector, 1, p)
        
        epsi <- matrix(rnorm(n), n, 1)
        Y <- apply(X, MARGIN = 1, dgp_function) + epsi
        
        if (p > 3) {
            W_0 <- tdnn:::feature_screen(X, Y)
        } else {
            W_0 <- rep(1, p)
        }
        if (verbose) print(W_0)
        
        # TDNN estimates ----------------------------------------------------------
        tdnn_rand_estimate <-
            tdnn:::tune_de_dnn_no_dist_vary_c_cpp_thread(
                X,
                Y,
                X_test_random,
                W0_ = W_0,
                c = c_seq,
                estimate_variance = est_variance,
                bootstrap_iter = 500,
                B_NN = 20,
                num_threads = num_threads
            )
        tdnn_fixed_estimate <-
            tdnn:::tune_de_dnn_no_dist_vary_c_cpp_thread(
                X,
                Y,
                X_test_fixed,
                W0_ = W_0,
                c = c_seq_fixed,
                estimate_variance = est_variance,
                bootstrap_iter = 500,
                B_NN = 20
            )
        # Collect TDNN results ----------------------------------------------------------
        tdnn_rand_results <-
            data.frame(
                estimate = tdnn_rand_estimate$estimate_loo,
                loss = (tdnn_rand_estimate$estimate_loo - mu_rand) ^ 2,
                bias = (tdnn_rand_estimate$estimate_loo - mu_rand),
                s_1 = tdnn_rand_estimate$s_1_B_NN,
                c = tdnn_rand_estimate$c_B_NN,
                variance = if (est_variance) {
                    as.numeric(tdnn_rand_estimate$variance)
                } else {
                    NA
                },
                method = "tdnn_rand"
            )
        tdnn_fixed_results <-
            data.frame(
                estimate = tdnn_fixed_estimate$estimate_loo,
                loss = (tdnn_fixed_estimate$estimate_loo - mu_fixed) ^ 2,
                bias = (tdnn_fixed_estimate$estimate_loo - mu_fixed),
                s_1 = tdnn_fixed_estimate$s_1_B_NN,
                c = tdnn_fixed_estimate$c_B_NN,
                variance = if (est_variance) {
                    as.numeric(tdnn_fixed_estimate$variance)
                } else {
                    NA
                },
                method = "tdnn_fixed"
            )
        # DNN estimates ----------------------------------------------------------
        dnn_fixed_estimate <-
            tdnn:::tune_dnn_no_dist_thread(
                X,
                Y,
                X_test_fixed,
                s_seq = dnn_s_seq,
                W0_ = W_0,
                estimate_variance = est_variance,
                num_threads = num_threads
            )
        dnn_rand_estimate <-
            tdnn:::tune_dnn_no_dist_thread(
                X,
                Y,
                X_test_random,
                s_seq = dnn_s_seq,
                W0_ = W_0,
                estimate_variance =
                    est_variance,
                num_threads = num_threads
            )
        
        # Collect DNN estimates ----------------------------------------------------------
        dnn_rand_results <-
            data.frame(
                estimate = dnn_rand_estimate$estimate_loo,
                loss = (dnn_rand_estimate$estimate_loo - mu_rand) ^ 2,
                bias = (dnn_rand_estimate$estimate_loo - mu_rand),
                s_1 = dnn_rand_estimate$s_1_B_NN,
                variance = if (est_variance) {
                    as.numeric(dnn_rand_estimate$variance)
                } else {
                    NA
                },
                method = "dnn_rand"
            )
        dnn_fixed_results <-
            data.frame(
                estimate = dnn_fixed_estimate$estimate_loo,
                loss = (dnn_fixed_estimate$estimate_loo - mu_fixed) ^ 2,
                bias = (dnn_fixed_estimate$estimate_loo - mu_fixed),
                s_1 = dnn_fixed_estimate$s_1_B_NN,
                variance = if (est_variance) {
                    as.numeric(dnn_fixed_estimate$variance)
                } else {
                    NA
                },
                method = "dnn_fixed"
            )
        # KNN estimates ----------------------------------------------------------
        # Use screened set of covariates
        X_filtered <- X[, as.logical(W_0)]
        X_test_random_filtered <- X_test_random[, as.logical(W_0)]
        X_test_fixed_filtered <- X_test_fixed[, as.logical(W_0)]
        knn_rand_estimate <-
            make_tuned_knn_results(
                X_filtered,
                Y,
                X_test_random_filtered,
                kmax = 200,
                estimate_variance = est_variance
            )
        knn_fixed_estimate <-
            make_tuned_knn_results(
                X_filtered,
                Y,
                X_test_fixed_filtered,
                kmax = 200,
                estimate_variance = est_variance
            )
        # Collect KNN estimates ----------------------------------------------------------
        knn_rand_results <-
            data.frame(
                estimate = knn_rand_estimate$pred_knn,
                loss = (knn_rand_estimate$pred_knn - mu_rand) ^ 2,
                bias = (knn_rand_estimate$pred_knn - mu_rand),
                k = knn_rand_estimate$best_k,
                variance = if (est_variance) {
                    as.numeric(knn_rand_estimate$variance)
                } else {
                    NA
                },
                method = "knn_rand"
            )
        knn_fixed_results <-
            data.frame(
                estimate = knn_fixed_estimate$pred_knn,
                loss = (knn_fixed_estimate$pred_knn - mu_fixed) ^
                    2,
                bias = (knn_fixed_estimate$pred_knn - mu_fixed),
                k = knn_fixed_estimate$best_k,
                variance = if (est_variance) {
                    as.numeric(knn_fixed_estimate$variance)
                } else {
                    NA
                },
                method = "knn_fixed"
            )
        # Combine results ----------------------------------------------------
        results <-
            bind_rows(
                tdnn_rand_results,
                tdnn_fixed_results,
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
        n_fixed,
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
        Estimate = mean(estimate),
        Estimate_SE = sd(estimate),
        Variance = mean(variance)
    ) %>%
    print()


results_data <- list(
    results_df = results,
    c_seq_fixed = c_seq_fixed,
    c_seq = c_seq,
    mu_rand = mu_rand,
    mu_fixed = mu_fixed,
    p = p,
    n = n_fixed,
    num_reps = num_reps,
    n_test = n_test,
    data_type = data_type,
    fixed_test_vector = fixed_test_vector,
    draw_random_data = draw_random_data
)
file_path <-
    file.path(OUTPUT_DIR, glue::glue(
        "tdnn_setting_1_{data_type}_data_{draw_random_data}_{str_replace_all(strftime(Sys.time(), '%Y-%m-%d_%H_%M'),'-', '_')}.rds")
    )
results_data %>%
    write_rds(., file_path)


# Make LaTeX tables -------------------------------------------------------
setting_1_results_df <- results %>%
    mutate(
        data_type = data_type,
        method = case_when(
            method == "tdnn_fixed" ~ "TDNN fixed vector",
            method == "tdnn_rand" ~ "TDNN random vector",
            method == "dnn_fixed" ~ "DNN fixed vector",
            method == "dnn_rand" ~ "DNN random vector",
            method == "knn_fixed" ~ "KNN fixed vector",
            method == "knn_rand" ~ "KNN random vector"
        )
    ) %>%
    rename(Method = method)

results_fixed <- setting_1_results_df %>%
    filter(str_ends(Method, "fixed vector")) %>%
    group_by(Method) %>%
    summarize(
        MSE = mean((estimate - mu_fixed) ^ 2),
        mean_est = mean(estimate),
        Variance = mean(variance),
        .groups = "keep"
    ) %>%
    summarize(
        Method = Method,
        MSE = MSE,
        `Squared Bias` = (mean_est - mu_fixed_normal) ^ 2,
        Variance = Variance,
    ) %>%
    mutate(Method = str_remove(Method, "fixed vector"))

results_random <- setting_1_results_df %>%
    ungroup() %>%
    filter(str_ends(Method, " random vector")) %>%
    mutate(point = rep(1:n_test, times = 3 * num_reps)) %>%
    group_by(Method, point) %>%
    summarize(
        MSE = mean(loss),
        mean_est = mean(estimate),
        Variance = mean(variance)
    ) %>%
    summarize(
        MSE = MSE,
        `Squared Bias` = (mean_est - mu_rand_normal) ^ 2,
        Variance = Variance
    ) %>%
    summarize(across(everything(), mean)) %>%
    select(-one_of(c("Method")))
combined_results <- bind_cols(results_fixed, results_random)

combined_results

# Save table to file ------------------------------------------------------
combined_results %>%
    knitr::kable(
        digits = c(0, rep(ROUNDING_DIGITS, 6)),
        "latex",
        booktabs = T,
        align = "lccccccccc"
    ) %>%
    kableExtra::add_header_above(c(
        " " = 1,
        "Fixed Test Point" = 3,
        "Random Test Points" = 3
    )) %>%
    gsub("\\.\\.\\.\\d", "", .) %>%
    kableExtra::save_kable(.,
                           file.path(OUTPUT_DIR,
                                     glue::glue("tdnn_setting_1_{data_type}_combined_results_table.tex")))

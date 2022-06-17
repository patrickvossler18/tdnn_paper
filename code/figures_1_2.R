library(tidyverse)
library(tdnn)
library(ggpubr)
library(cowplot)

# Plotting helper functions and constants ---------------------------------
OUTPUT_DIR <- "../output"
AXIS_TEXT_SIZE <- 22
HLINE_TEXT_SIZE <- 8

save_plot_pdf <-
    function(plot_obj,
             filename,
             width = 7,
             height = 7) {
        pdf(filename, width = width, height = height)
        print(plot_obj)
        dev.off()
    }

save_group_plot <-
    function(plot_obj_list,
             filename,
             width = 14,
             height = 14) {
        pdf(filename, width = width, height = height)
        print(plot_grid(plotlist = plot_obj_list))
        dev.off()
    }

save_group_plot_inset <-
    function(plot_obj_list,
             inset_plot_obj,
             filename,
             width = 14,
             height = 14) {
        grid_plots <- plot_grid(plotlist = plot_obj_list)
        pdf(filename, width = width, height = height)
        print()
        dev.off()
    }


# Algorithm functions -----------------------------------------------------

dgp_function <- function(x) {
    # This function will take a row of the matrix as input and return the
    # transformed value. To use this on a matrix, we can use the apply function
    # with margin=1
    (x[1] - 1) ^ 2 + (x[2] + 1) ^ 3 - 3 * x[3]
}

dnn <-
    function(X, Y,
             X_test,
             s.size) {
        n <- nrow(X)
        p <- ncol(X)
        # Weight
        ord <- matrix(1:n, n, 1)
        # (n-k s-1) over (n s)
        weight <- choose(n - ord, s.size - 1) / choose(n, s.size)
        # Distance
        X.dis <- X - kronecker(matrix(1, n, 1), X_test)
        EuDis <- (X.dis ^ 2) %*% matrix(1, p, 1)
        # Ascending small->large
        noise <- matrix(rnorm(1), n, 1)
        TempD <- data.frame(EuDis, Y, noise)[order(EuDis, noise),]
        # Estimator
        U <- sum(TempD$Y * weight)
        return(U)
    }

de.dnn <- function(X,
                   Y,
                   X.test,
                   s.size = 2,
                   bc.p = 2) {
    n <- nrow(X)
    p <- ncol(X)
    # Weight
    ord <- matrix(1:n, n, 1)
    # (n-k s-1) over (n s)
    weight1 <- choose(n - ord, s.size - 1) / choose(n, s.size)
    weight2 <-
        choose(n - ord, bc.p * s.size - 1) / choose(n, bc.p * s.size)
    # Distance
    X.dis <- X - kronecker(matrix(1, n, 1), X.test)
    EuDis <- (X.dis ^ 2) %*% matrix(1, p, 1)
    # Ascending small->large
    noise <- matrix(rnorm(1), n, 1)
    TempD <- data.frame(EuDis, Y, noise)[order(EuDis, noise),]
    # Estimator
    U1 <- sum(TempD$Y * weight1)
    U2 <- sum(TempD$Y * weight2)
    Magic <-
        solve(matrix(c(1, 1, 1, (1 / bc.p) ^ (2 / min(
            p, 3
        ))), 2, 2)) %*% matrix(c(1, 0), 2, 1)
    U <- Magic[1, 1] * U1 + Magic[2, 1] * U2
    return(U)
}

knn <- function(X, Y,
                X_test,
                k.size = 2) {
    n <- nrow(X)
    p <- ncol(X)
    # Weight
    ord <- matrix(1:n, n, 1)
    weight <- as.numeric(ord <= k.size) * 1 / k.size
    # Distance
    X.dis <- X - kronecker(matrix(1, n, 1), X_test)
    EuDis <- (X.dis ^ 2) %*% matrix(1, p, 1)
    # Ascending small->large
    noise <- matrix(rnorm(n), n, 1)
    TempD <- data.frame(EuDis, Y, noise)[order(EuDis, noise),]
    # Estimator
    U <- sum(TempD$Y * weight)
    return(U)
}

# Simulation Parameters ---------------------------------------------------
n <- 1000
s_1_seq <- seq(1, 200, 1)
p <- 3

seed_val <- 1234

fixed_test_vector <- c(0.5, -0.5, 0.5)
X_test_fixed <- matrix(fixed_test_vector, 1, p)
mu <- dgp_function(fixed_test_vector)
fixed_c <- 2
param_df <- expand_grid(c = fixed_c, s_1 = s_1_seq)

k_grid <- seq(1, 200, 1)

n_reps <- 500

# Run simulation loop ----------------------------------------------------------
set.seed(seed_val)

rep_results <- map_df(seq(n_reps), function(rep) {
    print(rep)
    X <- matrix(rnorm(n * p), n, p)
    epsi <- matrix(rnorm(n), n, 1)
    Y <- apply(X, MARGIN = 1, dgp_function) + epsi
    
    vary_c_m_results <- pmap_df(param_df, function(c, s_1) {
        bind_rows(list(
            data.frame(
                estimate = dnn(X, Y, X_test_fixed, s_1),
                s_1 = s_1,
                type = "dnn"
            ),
            data.frame(
                estimate = de.dnn(X, Y, X_test_fixed, s_1, bc.p = c),
                s_1 = s_1,
                type = "de_dnn"
            )
        ))
    })
    
    # Tuned de-DNN results
    tuned_results <-
        tdnn:::tune_de_dnn(X,
                           Y,
                           X_test_fixed,
                           c = 2,
                           B_NN = 20,
                           scale_p = 1)
    tuned_df <-
        data.frame(
            c = fixed_c,
            s_1 = tuned_results$s_1_B_NN,
            estimate = tuned_results$estimate_loo,
            type = "tuned"
        )
    
    # k-NN over grid of k values
    knn_results <- map_dfr(k_grid, function(k) {
        knn_est <- knn(X, Y, X_test_fixed, k.size = k)
        data.frame(
            c = NA,
            k = k,
            s_1 = NA,
            estimate = knn_est,
            type = "knn"
        )
    })
    bind_rows(vary_c_m_results, tuned_df, knn_results) %>% mutate(id = rep)
})



# Calculate aggregated statistics -----------------------------------------

# TUNED RESULTS
tuned_rep_results <- rep_results %>% filter(type == "tuned") %>%
    summarize(MSE = mean((estimate-mu)^2),
              Bias = mean((estimate - mu)))

# KNN RESULTS
knn_rep_results <- rep_results %>% filter(type == "knn") %>% group_by(k, type) %>% 
    summarize(MSE = mean((estimate-mu)^2),
              Bias = mean((estimate - mu)))

# DNN and DE-DNN RESULTS
dnn_rep_df <- rep_results %>% filter(type == "dnn") %>% group_by(s_1,type) %>% summarize(
    MSE = mean((estimate-mu)^2),
    Bias = mean(estimate-mu)
)

de_dnn_rep_df <- rep_results %>% filter(type == "de_dnn") %>% group_by(s_1, type) %>% summarize(
    MSE = mean((estimate-mu)^2),
    Bias = mean(estimate-mu)
)

min_dnn_MSE <- dnn_rep_df %>% ungroup() %>% filter(MSE == min(MSE)) 
min_de_dnn_MSE <- de_dnn_rep_df %>% ungroup() %>% filter(MSE == min(MSE))
min_knn_MSE <- knn_rep_results %>% ungroup() %>% filter(MSE == min(MSE))

dnn_de_dnn_data <- bind_rows(list(dnn_rep_df, de_dnn_rep_df))

min_dnn_s <- min_dnn_MSE$s_1

# FIGURE 1 PLOTS -------------------------------------------------------------------------

# DNN MSE plot -------------------------------------------------------------------------
dnn_plot_full <- dnn_de_dnn_data %>%
    filter(type %in% c("dnn")) %>% ungroup() %>% mutate(type = "DNN") %>% 
    ggplot(aes(x = s_1, y = MSE, color = type)) + geom_line(size=1) +
    geom_texthline(aes(yintercept = min_dnn_MSE$MSE, linetype = "dnn", label = glue::glue("DNN minimum = {round(min_dnn_MSE$MSE,4)}")),
                   hjust = 1, vjust = -0.2, size = HLINE_TEXT_SIZE) +
    scale_linetype_manual(
        name = "Minimum value",
        values = c(3)
    ) +
    scale_color_manual(
        values = c("DNN" = "black")
    )+
    guides(color = "none", linetype = "none") + 
    xlab("Subsampling scale s") +
    ylab("Mean squared error (MSE)") +
    theme_pubr(base_size = AXIS_TEXT_SIZE)

dnn_plot_full

# DNN bias plot -------------------------------------------------------------------------
dnn_bias_plot <- dnn_de_dnn_data %>%
    filter(type %in% c("dnn")) %>% ungroup() %>% mutate(type = "DNN") %>% 
    ggplot(aes(x = s_1, y = Bias, color = type)) + geom_line(size=1) +
    scale_color_manual(
        values = c("DNN" = "black")
    )+
    guides(color = "none")+
    xlab("Subsampling scale s") +
    ylab("Bias of DNN") +
    theme_pubr(base_size = AXIS_TEXT_SIZE,legend = "none")
dnn_bias_plot

## TDNN MSE plot -------------------------------------------------------------------------
min_de_dnn_s <- min_de_dnn_MSE$s_1
# dnn, de-dnn, tuned de-dnn zoomed in

tdnn_full_zoom_plot <- dnn_de_dnn_data %>%
    filter(between(s_1, min_de_dnn_s - 5, min_de_dnn_s + 40)) %>%
    filter(type %in% c("de_dnn")) %>% ungroup() %>% mutate(type = "TDNN") %>% 
    ggplot(aes(x = s_1, y = MSE, color = type)) + geom_line(size=1) +
    geom_texthline(aes(yintercept = min_de_dnn_MSE$MSE, linetype = "de-dnn", label = glue::glue("TDNN minimum = {round(min_de_dnn_MSE$MSE,4)}")),
                   hjust = 1, vjust = 0.5, size = HLINE_TEXT_SIZE)+
    geom_texthline(aes(yintercept = min_dnn_MSE$MSE, linetype = "dnn", label = glue::glue("DNN\n minimum = {round(min_dnn_MSE$MSE,4)}")),
                   hjust = 1, vjust = 0, size = HLINE_TEXT_SIZE) +
    geom_texthline(aes(yintercept = tuned_rep_results$MSE, linetype = "tuned de-dnn", label = glue::glue("tuned TDNN minimum = {round(tuned_rep_results$MSE,4)}")),
                   hjust = 1, vjust = 0, size = HLINE_TEXT_SIZE)+
    scale_linetype_manual(
        name = "Minimum value",
        values = c(2, 3, 4),
    ) +
    scale_color_manual(
        values = c("TDNN" = "black")
    )+
    guides(color = "none", linetype = "none") + 
    xlab("Subsampling scale s") +
    ylab("Mean squared error (MSE)") +
    scale_x_continuous(breaks = seq(0,45,5))+
    theme_pubr(base_size = AXIS_TEXT_SIZE)

tdnn_full_zoom_plot

## TDNN bias plot -------------------------------------------------------------------------
tdnn_bias_plot <- dnn_de_dnn_data %>%
    filter(between(s_1, 1, min_de_dnn_s + 40)) %>%
    filter(type %in% c("de_dnn")) %>% ungroup() %>% mutate(type = "TDNN") %>% 
    ggplot(aes(x = s_1, y = Bias, color = type)) + geom_line(size=1) +
    scale_color_manual(
        values = c("TDNN" = "black")
    )+
    guides(color = "none", linetype = "none") +
    xlab("Subsampling scale s") +
    ylab("Bias of TDNN") +
    scale_x_continuous(breaks = seq(0,45,5))+
    theme_pubr(base_size = AXIS_TEXT_SIZE)
tdnn_bias_plot


# dnn tdnn grid plot with inset -------------------------------------------

dnn_plot_zoomed_in_mod <- dnn_de_dnn_data %>%
    filter(between(s_1, min_dnn_s - 40, min_dnn_s + 50)) %>%
    filter(type %in% c("dnn")) %>% ungroup() %>% mutate(type = "DNN") %>% 
    ggplot(aes(x = s_1, y = MSE, color = type)) + geom_line(size=1) +
    geom_hline(aes(yintercept = min_dnn_MSE$MSE, linetype = "dnn")) +
    scale_linetype_manual(
        name = "Minimum value",
        values = c(3)
    ) +
    scale_color_manual(
        values = c("DNN" = "black")
    )+
    guides(color = "none", linetype = "none") + 
    xlab("Subsampling scale s") +
    ylab("Mean squared error (MSE)") +
    theme_pubr(base_size = AXIS_TEXT_SIZE) +
    theme(axis.title=element_blank())

dnn_tdnn_grid_plot_inset <- ggdraw(plot_grid(plotlist = list(dnn_bias_plot, dnn_plot_full, tdnn_bias_plot, tdnn_full_zoom_plot)))+ 
    draw_plot(dnn_plot_zoomed_in_mod , x = .3, y = .32, scale = 0.37)

dnn_tdnn_grid_plot_inset

# FIGURE 2 PLOTS -------------------------------------------------------------------------
# KNN MSE plot -------------------------------------------------------------------------
knn_mse_plot_full <- knn_rep_results %>%
    ungroup() %>% mutate(type = "KNN") %>% 
    ggplot(aes(x = k, y = MSE, color = type)) + geom_line(size=1) +
    geom_texthline(aes(yintercept = min_knn_MSE$MSE, linetype = "KNN", label = glue::glue("KNN minimum = {round(min_knn_MSE$MSE,4)}")),
                   hjust = 1, vjust = 0.3, size = HLINE_TEXT_SIZE)+
    scale_linetype_manual(
        name = "Minimum value",
        values = c(5),
    ) +
    scale_color_manual(
        values = c("KNN" = "black")
    )+
    guides(color = "none", linetype = "none") + 
    xlab("Neighborhood size k") +
    ylab("Mean squared error (MSE)") +
    # scale_x_continuous(breaks = seq(0,60,5))+
    theme_pubr(base_size = AXIS_TEXT_SIZE)

knn_mse_plot_full

# KNN Bias Plot -------------------------------------------------------------------------
knn_bias_plot <- knn_rep_results %>%
    filter(between(k, min_k - 10, min_k + 50)) %>% ungroup() %>% mutate(type = "KNN") %>% 
    ggplot(aes(x = k, y = Bias, color = type)) + geom_line(size=1) +
    scale_color_manual(
        values = c("KNN" = "black")
    )+
    guides(color = "none", linetype = "none") +
    xlab("Neighborhood size k") +
    ylab("Bias of KNN") +
    scale_x_continuous(breaks = seq(0,60,5))+
    theme_pubr(base_size = AXIS_TEXT_SIZE)

knn_bias_plot


# Save plots to file -------------------------------------------------------------------------
save_plot_pdf(dnn_tdnn_grid_plot_inset,
              file.path(OUTPUT_DIR,"figure_1.pdf"),width=16,height=14)
save_group_plot(list(knn_bias_plot, knn_mse_plot_full),
                file.path(OUTPUT_DIR, "figure_2.pdf"), width=16, height = 7)

library(tidyverse)
library(parallel)
library(energy)


################
# Abalone data #
################
# replace I with 0, M with 1, F with 2
# The Abalone data can be found at https://archive.ics.uci.edu/ml/datasets/abalone

aba.data <- read.table("abalone.txt", sep = ",") # X in R^8 and Y in R
X <- as.matrix(aba.data[, 1:8])
Y <- as.matrix(aba.data[, 9], ncol = 1)


# DNN estimator with the input of ordered Y according to the distance of corresponding feature vector to a  given point#
dnn <- function(ordered_Y, n, p, s.size = 2) {
  # Weight
  ord1 <- matrix(1:(n - s.size + 1), n - s.size + 1, 1)
  # (n-k s-1) over (n s)
  weight1 <- rbind(s.size * exp(lgamma(n - ord1 + 1) + lgamma(n - s.size + 1) - lgamma(n + 1) - lgamma(n - ord1 - s.size + 2)), matrix(0, nrow = s.size - 1, ncol = 1)) # choose(n - ord, s.size - 1) / choose(n, s.size)
  U1 <- sum(ordered_Y * weight1)
  return(U1)
}

# TDNN estimator with the input of ordered Y according to the distance of corresponding feature vector to a  given point#
de.dnn <- function(ordered_Y,
                   n,
                   p,
                   s.size = 2,
                   bc.p = 2) {
  # Weight
  ord1 <- matrix(1:(n - s.size + 1), n - s.size + 1, 1)
  ord2 <- matrix(1:(n - ceiling(bc.p * s.size) + 1), n - ceiling(bc.p * s.size) + 1, 1)
  # (n-k s-1) over (n s)
  weight1 <- rbind(s.size * exp(lgamma(n - ord1 + 1) + lgamma(n - s.size + 1) - lgamma(n + 1) - lgamma(n - ord1 - s.size + 2)), matrix(0, nrow = s.size - 1, ncol = 1)) # choose(n - ord, s.size - 1) / choose(n, s.size)
  weight2 <- rbind(ceiling(bc.p * s.size) * exp(lgamma(n - ord2 + 1) + lgamma(n - ceiling(bc.p * s.size) + 1) - lgamma(n + 1) - lgamma(n - ord2 - ceiling(bc.p * s.size) + 2)), matrix(0, nrow = ceiling(bc.p * s.size) - 1, ncol = 1)) # choose(n - ord, bc.p * s.size - 1) / choose(n, bc.p * s.size)

  # Estimator
  U1 <- sum(ordered_Y * weight1)
  U2 <- sum(ordered_Y * weight2)
  Magic <- solve(matrix(c(1, 1, 1, (1 / bc.p)^(2 / p)), 2, 2)) %*% matrix(c(1, 0), 2, 1)
  U <- Magic[1, 1] * U1 + Magic[2, 1] * U2
  return(U)
}

# DNN estimator with the input of X, Y and the given point#
de.dnn0 <- function(X,
                    Y,
                    X.test,
                    s.size = 2,
                    bc.p = 2) {
  n <- nrow(X)
  p <- ncol(X)
  # Weight
  ord1 <- matrix(1:(n - s.size + 1), n - s.size + 1, 1)
  ord2 <- matrix(1:(n - ceiling(bc.p * s.size) + 1), n - ceiling(bc.p * s.size) + 1, 1)
  # (n-k s-1) over (n s)
  weight1 <- rbind(s.size * exp(lgamma(n - ord1 + 1) + lgamma(n - s.size + 1) - lgamma(n + 1) - lgamma(n - ord1 - s.size + 2)), matrix(0, nrow = s.size - 1, ncol = 1)) # choose(n - ord, s.size - 1) / choose(n, s.size)
  weight2 <- rbind(ceiling(bc.p * s.size) * exp(lgamma(n - ord2 + 1) + lgamma(n - ceiling(bc.p * s.size) + 1) - lgamma(n + 1) - lgamma(n - ord2 - ceiling(bc.p * s.size) + 2)), matrix(0, nrow = ceiling(bc.p * s.size) - 1, ncol = 1)) # choose(n - ord, bc.p * s.size - 1) / choose(n, bc.p * s.size)
  # Distance
  X.dis <- X - kronecker(matrix(1, n, 1), X.test)
  EuDis <- (X.dis^2) %*% matrix(1, p, 1)
  # Ascending small->large
  noise <- matrix(rnorm(1), n, 1)
  TempD <- data.frame(EuDis, Y, noise)[order(EuDis, noise), ]
  # Estimator
  U1 <- sum(TempD$Y * weight1)
  U2 <- sum(TempD$Y * weight2)
  Magic <- solve(matrix(c(1, 1, 1, (1 / bc.p)^(2 / p)), 2, 2)) %*% matrix(c(1, 0), 2, 1)
  U <- Magic[1, 1] * U1 + Magic[2, 1] * U2
  return(U)
}

# TDNN estimator with the input of X, Y and the given point#
dnn0 <- function(X,
                 Y,
                 X.test,
                 s.size = 2) {
  n <- nrow(X)
  p <- ncol(X)
  # Weight
  ord1 <- matrix(1:(n - s.size + 1), n - s.size + 1, 1)
  # (n-k s-1) over (n s)
  weight1 <- rbind(s.size * exp(lgamma(n - ord1 + 1) + lgamma(n - s.size + 1) - lgamma(n + 1) - lgamma(n - ord1 - s.size + 2)), matrix(0, nrow = s.size - 1, ncol = 1)) # choose(n - ord, s.size - 1) / choose(n, s.size)
  # Distance
  X.dis <- X - kronecker(matrix(1, n, 1), X.test)
  EuDis <- (X.dis^2) %*% matrix(1, p, 1)
  # Ascending small->large
  noise <- matrix(rnorm(1), n, 1)
  TempD <- data.frame(EuDis, Y, noise)[order(EuDis, noise), ]
  # Estimator
  U1 <- sum(TempD$Y * weight1)
  return(U1)
}




## PCA: dimension reduction
library("stats")
pca.f <- function(X, m) {
  pca.X <- prcomp(X)
  pc <- X %*% (pca.X$rotation)
  pc <- as.matrix(pc[, 1:m], ncol = m)
  var.ratio <- sum(((pca.X$sdev)^2)[1:m]) / sum(pca.X$sdev^2)
  return(list(pc = pc, var.ratio = var.ratio))
}


#############
# rescale X #
#############
scaling <- function(X) {
  n <- nrow(X)
  scaled_X <- sqrt(n) * X %*% diag(apply(X^2, 2, sum)^(-1 / 2))
  return(scaled_X)
}



########################################################
########## approximate point-wise MSE ##################
########################################################

tune_s <- function(X, Y, X_test, s_seq, c) {
  t <- length(s_seq)
  tuning <- matrix(0, length(s_seq), 1)
  tuning <- sapply(s_seq, function(s) {
    de.dnn0(X, Y, X_test, s.size = s + 1, bc.p = c)
  })
  s.choice <- which(diff(abs(diff(tuning) / tuning[1:t - 1])) > -0.01)[1] + 3
  if (is.na(s.choice)) {
    s.choice <- 3
  }
  return(s.choice)
}

tune.f <- function(train.set.X1, train.set.Y1, test.set.X1, test.set.Y1, c_seq, s_seq, B, scale_p = 1) {
  train.number1 <- nrow(train.set.X1)
  p <- ncol(train.set.X1)

  tdnn_dnn_results <- bind_rows(lapply(1:nrow(test.set.X1), function(l) {
    if (l %% 20 == 0) {
      print(l)
    }
    rrr <- bind_rows(lapply(c_seq, function(fixed.c) {
      max_s_1 <- trunc((train.number1 - 3) / fixed.c) - 1
      s_1_seq0_1 <- seq(1, min(round(sqrt(train.number1)), max_s_1), 1)
      s.re1 <- tune_s(train.set.X1, train.set.Y1, matrix(test.set.X1[l, ], 1, p), s_1_seq0_1, fixed.c)
      estimate_re1 <- de.dnn0(train.set.X1, train.set.Y1, matrix(test.set.X1[l, ], 1, p), s.size = s.re1, bc.p = fixed.c)
      s_1_seq1 <- seq(min(ceiling(s.re1), max_s_1), min(ceiling(s.re1 * 2), max_s_1), 1)
      MSE.re1 <- (estimate_re1 - test.set.Y1[l, ])^2

      param_df1 <- tidyr::expand_grid(c = fixed.c, s_1 = s_1_seq1)
      X.dis1 <- train.set.X1 - kronecker(matrix(1, train.number1, 1), matrix(test.set.X1[l, ], 1, p))
      EuDis1 <- (X.dis1^2) %*% matrix(1, p, 1)
      B.index1 <- (sort(EuDis1, index.return = T)$ix)[1:B]

      boot_reps1 <- lapply(1:B, function(b) {
        # get the bth LOO sample indices
        # Make train data and validation data
        X_train <- train.set.X1[-B.index1[b], ]
        Y_train <- as.matrix(train.set.Y1[-B.index1[b], ], ncol = 1)

        X_val <- matrix(train.set.X1[B.index1[b], ], 1, ncol(train.set.X1))
        Y_val <- as.matrix(train.set.Y1[B.index1[b], ], ncol = 1)

        n_train <- nrow(X_train)
        p_train <- ncol(X_train)

        X_train_dis <- X_train - kronecker(matrix(1, n_train, 1), X_val)
        EuDis_train <- (X_train_dis^2) %*% matrix(1, p_train, 1)
        index.order <- sort(EuDis_train, index = T)$ix
        ordered_Y_train <- Y_train[index.order]

        # for each LOO sample, calculate estimates for each s_1 value
        # this loops through the parameter combinations and returns a data frame with the results
        tdnn.res <- pmap_df(param_df1, function(c, s_1) {
          param_estimate <- de.dnn(ordered_Y_train, n_train, p_train, s_1, c)
          bind_rows(list(
            data.frame(
              estimate = param_estimate,
              s_1 = s_1,
              c = c,
              y_val = Y_val,
              MSE.tdnn = (param_estimate - Y_val)^2 * exp(-sum((X_val - matrix(test.set.X1[l, ], 1, p))^2) / scale_p)
            )
          ))
        })

        dnn.res <- pmap_df(tidyr::expand_grid(c = 1, s = s_seq), function(c, s) {
          est.dnn <- dnn(ordered_Y_train, n_train, p_train, s)
          bind_rows(list(data.frame(
            MSE.dnn = (est.dnn - Y_val)^2 * exp(-sum((X_val - matrix(test.set.X1[l, ], 1, p))^2) / scale_p),
            s = s,
            y_val = Y_val
          )))
        })
        return(list(tdnn.res = tdnn.res, dnn.res = dnn.res))
      })

      boot_rep_results_tdnn <- matrix(0, nrow = length(c_seq) * length(s_1_seq1), ncol = 5)
      boot_rep_results_dnn <- matrix(0, nrow = length(s_seq), ncol = 3)
      for (i in 1:B) {
        u <- boot_reps1[[i]]
        boot_rep_results_tdnn <- boot_rep_results_tdnn + u$tdnn.res
        boot_rep_results_dnn <- boot_rep_results_dnn + u$dnn.res
      }
      boot_rep_results_tdnn <- boot_rep_results_tdnn / B
      boot_rep_results_dnn <- boot_rep_results_dnn / B
      tdnn.min <- boot_rep_results_tdnn[which.min(boot_rep_results_tdnn[, 5]), ]
      dnn.min <- boot_rep_results_dnn[which.min(boot_rep_results_dnn[, 1]), ]

      tdnn.min %>% mutate(
        dnn.MSE = dnn.min$MSE.dnn,
        dnn.s = dnn.min$s,
        fixed.c = fixed.c
      )
    }))
    tdnn.min <- rrr[rrr$MSE.tdnn == min(rrr$MSE.tdnn), ]
    tdnn.tuned.MSE1 <- (de.dnn0(train.set.X1, train.set.Y1, matrix(test.set.X1[l, ], 1, p), tdnn.min$s_1, tdnn.min$c) - test.set.Y1[l, ])^2
    dnn.tuned.MSE1 <- (dnn0(train.set.X1, train.set.Y1, matrix(test.set.X1[l, ], 1, p), tdnn.min$dnn.s) - test.set.Y1[l, ])^2

    data.frame(
      tdnn.tuned.MSE = tdnn.tuned.MSE1,
      dnn.tuned.MSE = dnn.tuned.MSE1,
      fixed.c = tdnn.min$fixed.c,
      c = tdnn.min$c,
      s_1 = tdnn.min$s_1,
      s = tdnn.min$dnn.s
    )
  }))
  return(tdnn_dnn_results)
}

real.mse.pointwise.pca <- function(X, Y, ratio = 0.75, B = 50, c_seq, s_seq, n.rep = 10, scale_p = 1, m) {
  label <- X[, 1]
  pca.result <- pca.f(X[, -1], m)
  X <- pca.result$pc
  var.ratio <- pca.result$var.ratio

  # three categories
  X <- cbind(label, X)
  X_1 <- X[X[, 1] == 1, -1]
  Y_1 <- as.matrix(Y[X[, 1] == 1, ], ncol = 1)
  n1 <- nrow(X_1)
  X_2 <- X[X[, 1] == 0, -1]
  Y_2 <- as.matrix(Y[X[, 1] == 0, ], ncol = 1)
  n2 <- nrow(X_2)
  X_3 <- X[X[, 1] == 2, -1]
  Y_3 <- as.matrix(Y[X[, 1] == 2, ], ncol = 1)
  n3 <- nrow(X_3)
  p <- ncol(X_1)
  # scaling
  X_1 <- scaling(X_1)
  X_2 <- scaling(X_2)
  X_3 <- scaling(X_3)

  train.number1 <- ceiling(ratio * n1)
  train.number2 <- ceiling(ratio * n2)
  train.number3 <- ceiling(ratio * n3)

  real.mse.res <- lapply(1:n.rep, function(j) {
    train.index1 <- sample(n1, train.number1)
    train.set.X1 <- X_1[train.index1, ]
    train.set.Y1 <- as.matrix(Y_1[train.index1, ], ncol = 1)
    test.set.X1 <- X_1[-train.index1, ]
    test.set.Y1 <- as.matrix(Y_1[-train.index1, ])

    train.index2 <- sample(n2, train.number2)
    train.set.X2 <- X_2[train.index2, ]
    train.set.Y2 <- as.matrix(Y_2[train.index2, ], ncol = 1)
    test.set.X2 <- X_2[-train.index2, ]
    test.set.Y2 <- as.matrix(Y_2[-train.index2, ])

    train.index3 <- sample(n3, train.number3)
    train.set.X3 <- X_3[train.index3, ]
    train.set.Y3 <- as.matrix(Y_3[train.index3, ], ncol = 1)
    test.set.X3 <- X_3[-train.index3, ]
    test.set.Y3 <- as.matrix(Y_3[-train.index3, ], ncol = 1)

    tdnn_dnn_results1 <- apply(tune.f(train.set.X1, train.set.Y1, test.set.X1, test.set.Y1, c_seq, s_seq, B, scale_p), 2, mean)
    tdnn_dnn_results2 <- apply(tune.f(train.set.X2, train.set.Y2, test.set.X2, test.set.Y2, c_seq, s_seq, B, scale_p), 2, mean)
    tdnn_dnn_results3 <- apply(tune.f(train.set.X3, train.set.Y3, test.set.X3, test.set.Y3, c_seq, s_seq, B, scale_p), 2, mean)

    return(list(tdnn_dnn_results1 = tdnn_dnn_results1, tdnn_dnn_results2 = tdnn_dnn_results2, tdnn_dnn_results3 = tdnn_dnn_results3))
  })

  tdnn_dnn1 <- rep(0, 6)
  tdnn_dnn2 <- rep(0, 6)
  tdnn_dnn3 <- rep(0, 6)
  for (i in 1:n.rep) {
    tdnn_dnn1 <- tdnn_dnn1 + real.mse.res[[i]]$tdnn_dnn_results1
    tdnn_dnn2 <- tdnn_dnn2 + real.mse.res[[i]]$tdnn_dnn_results2
    tdnn_dnn3 <- tdnn_dnn3 + real.mse.res[[i]]$tdnn_dnn_results3
  }
  tdnn_dnn1 <- tdnn_dnn1 / n.rep
  tdnn_dnn2 <- tdnn_dnn2 / n.rep
  tdnn_dnn3 <- tdnn_dnn3 / n.rep

  return(rbind(tdnn_dnn1, tdnn_dnn2, tdnn_dnn3))
}



## setting of parameters
ratio <- 0.75
scale_p <- 1
B <- 50
c_seq <- c(1.2, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20)
s_seq <- seq(from = 50, to = 250, by = 5)

# We repeat the data splitting procedure for 50 times and obtain averaged prediction errors over those 50 repetitions
### To increase computing efficiency, we separate the 50 repetitions into 50 parallel computation, and the random seeds are chosen from seq(from=1234, to = 1724, by = 10)
## In our computation, the following procedure was run in parallel in 50 different R sessions and used a different seed from seq(from=1234, to = 1724, by = 10) for different R sessions, then we take average of the results from 50 parallel R sessions
set.seed(1234)
n.rep <- 1
res.point_pca_m_3 <- real.mse.pointwise.pca(X, Y, ratio, B, c_seq, s_seq, n.rep, scale_p, m = 3)

X_1 <- X[X[, 1] == 1, -1]
X_2 <- X[X[, 1] == 0, -1]
X_3 <- X[X[, 1] == 2, -1]
n1 <- nrow(X_1)
n2 <- nrow(X_2)
n3 <- nrow(X_3)
tdnn_dnn_results <- (res.point_pca_m_3[1, 1:2] * n1 + res.point_pca_m_3[2, 1:2] * n2 + res.point_pca_m_3[3, 1:2] * n3) / (n1 + n2 + n3)

## print the results of prediction errors for TDNN and DNN
tdnn_dnn_results

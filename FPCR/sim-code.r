## Simulations for Fast NCV
## Shuyang (Henry) Cao
## June 2026
## File 03a: Code for simulation of FPCR for Fast NCV 

if (!exists("FAST_DEMO")) FAST_DEMO <- FALSE  
# Define parameters dynamically based on mode
if (FAST_DEMO) {
  message("--> RUNNING IN FAST_DEMO MODE (Toy Example)")
  r <- 3 # Number of replications for resampling in NCV
} else {
  message("--> RUNNING IN PRODUCTION MODE (Full Manuscript Runs)")
  r <- 200  # Number of replications for resampling in NCV
}

simname <- "sims-FPCR"
cat("R is chosen as", r, "for nested cross-validation algorithm. \n")

## import packages
pkgs <- c(
    'parallel', 'Matrix', 'refund', 'fda', 'mgcv', 
    'mvtnorm'
)
suppressPackageStartupMessages({
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste0("Could not find package ", pkg, ", installing from CRAN.\n"))
      
      # Set a secure CRAN mirror automatically so the terminal doesn't pause to ask you to choose one
      install.packages(pkg, repos = "https://cloud.r-project.org/")
      
      # Try loading it again, quietly
      require(pkg, character.only = TRUE, quietly = TRUE)
    }
  }
})

## Set paths 
basepath <- getwd()
resultspath <- file.path(basepath, "results")
if (!dir.exists(resultspath)) dir.create(resultspath)
figurespath <- file.path(basepath, "figures")
if (!dir.exists(figurespath)) dir.create(figurespath)
simresultsname <- paste0(simname, ".RData")

## Define all the functions needed
# Loss function -- square loss
sqloss <- function(y1, y2) {
    return((y1-y2)^2)
}

# NCV function for one run
fastncv <- function(X, y, model_mat, penalty_mat, K, loss=sqloss, lambda, store_mat, Res) {
    # Randomly assign K folds
    n <- length(y)
    nk <- n / K
    idx <- sample.int(n)
    kidx <- rep(1:K,each=nk)
    foldidx <- split(idx,kidx) 
    # construct identity matrix for usage later
    I <- diag(nk)
    II <- diag(2 * nk)
    # compute all residues for one-fold and store them in an array
    res_array <- array(0, dim = c(K, K, nk)) 
    # start computing residuals
    for (i in 1:K) {
        for (j in i:K) {
            if (j == i) {
                temp_idx <- c(foldidx[[i]])
                res_array[i, i, ] <- solve(I - tcrossprod(model_mat[temp_idx, ] %*% store_mat, model_mat[temp_idx, ]), Res[temp_idx])
            } else {
                temp_idx <- c(foldidx[[i]], foldidx[[j]])
                temp_res <- solve(II - tcrossprod(model_mat[temp_idx, ] %*% store_mat, model_mat[temp_idx, ]), Res[temp_idx])
                res_array[i, j, ] <- temp_res[(nk+1):(2*nk)]
                res_array[j, i, ] <- temp_res[1:nk]
            }
        }
    }
    # computev NCV prediction interval
    ## First part of the formula
    part_a <- lapply(1:K, function(idx) {
        (mean(loss(res_array[idx, idx, ], 0)) - mean(loss(res_array[idx, -idx, ], 0)))^2
    })
    ## Second part of the formula
    part_b <- lapply(1:K, function(idx) {
        var(loss(res_array[idx, idx, ], 0)) / nk
    })
    ## compute prediction interval
    pred_int_vec <- unlist(part_a) - unlist(part_b)
    pred_int <- mean(pred_int_vec)
    ## compute point estimate
    pe_vec <- unlist(lapply(1:K, function(idx) {
        loss(res_array[idx, -idx, ], 0)
    }))
    pe <- mean(pe_vec)
    return(list(ncv_predint_vec = pred_int_vec, ncv_predint = pred_int, 
    ncv_pe_vec = pe_vec, ncv_pe = pe))
}

# Main NCV function with parallelization
fastncv_pal <- function(X, y, model_mat, penalty_mat, K, loss=sqloss, lambda, R) {
    # This is where we fit the model ONCE on the full data.
    # However, we do not want to store any n by n hat matrix.
    # Construct the penalty matrix and store it for later usage
    store_mat <- solve(crossprod(model_mat) + lambda * penalty_mat)
    # compute residuals of the full model
    Res <- y - model_mat %*% (store_mat %*% crossprod(model_mat, y))

    all <- mclapply(1:R, function(i){fastncv(X, y, model_mat, penalty_mat, K, loss, lambda, store_mat, Res)}, 
    mc.cores = detectCores() - 1)
    ncv_preds <- c()
    ncv_pes <- c()
    for (i in 1:R) {
        ncv_preds <- c(ncv_preds, all[[i]]$ncv_predint_vec)
        ncv_pes <- c(ncv_pes, all[[i]]$ncv_pe_vec)
     }

     ## Bias Correction 
     cv_means <- c()
     repos <- R / 5 ## less repeats than NCV;
     ## reference to Bates et al. (2024) expression (15)
     for (r in 1:repos) {
          n <- length(y)
          nk <- n / K
          idx <- sample.int(n)
          kidx <- rep(1:K,each=nk)
          foldidx <- split(idx,kidx) 
          # construct identity matrix for usage later
          I <- diag(nk)
          for (i in 1:K) {
                temp_idx <- c(foldidx[[i]])
                one_fold_res <- solve(I - tcrossprod(model_mat[temp_idx, ] %*% store_mat, model_mat[temp_idx, ]), Res[temp_idx])
                cv_means <- c(cv_means, loss(one_fold_res, 0))
          }
     }
return(list(ncv_predint = max(0, mean(ncv_preds)), ncv_pe = mean(ncv_pes), 
     bias_correct = (1 + (K - 2) / K) * (mean(ncv_pes) - mean(cv_means))))
}

# Function for fitting FPCR model
fit_fpcr <- function(X, y, pci, spline_num = 40, lam) {
  N <- nrow(X)
  P <- ncol(X)
  
  # 1. Create the evaluation grid across the SIGNAL FEATURES (1 to P)
  signal_grid <- seq(0, 1, length.out = P)
  helper_dat <- data.frame(x = signal_grid)
  
  # 2. Generate the basis matrix B (P x spline_num) and Penalty S (spline_num x spline_num)
  # We use absorb.cons = FALSE to keep the raw basis functions intact
  helper_obj <- smoothCon(s(x, bs = "bs", k = spline_num, fx = FALSE), 
                          data = helper_dat, absorb.cons = FALSE)
  
  B <- helper_obj[[1]]$X        # Dimensions: P x spline_num
  S <- helper_obj[[1]]$S[[1]]   # Dimensions: spline_num x spline_num
  
  # 3. Map X into the Spline Space
  # Z matrix dimensions: N x spline_num
  Z_spline <- X %*% B 
  
  # 4. Perform SVD on the spline-mapped space for Dimension Reduction (FPCR)
  Z_svd <- svd(Z_spline)
  eigen_values <- (Z_svd$d)^2
  V <- Z_svd$v                  # Dimensions: spline_num x spline_num
  VA <- V[, 1:pci, drop = FALSE] # Take first 'pci' components (spline_num x pci)
  variance_proportion <- sum(eigen_values[1:pci]) / sum(eigen_values)
  # 5. Final Reduced Predictor Matrix (N x pci)
  Z_reduced <- Z_spline %*% VA
  
  # 6. Penalized Least Squares for the Reduced Coefficients
  # Penalty is projected into the reduced PCA space: VA^T %*% S %*% VA
  S_reduced <- crossprod(VA, S %*% VA)
  
  # Solve: (Z_red^T %*% Z_red + lam * S_red)^-1 %*% Z_red^T %*% y
  pci_coef <- solve(crossprod(Z_reduced) + lam * S_reduced) %*% crossprod(Z_reduced, y)
  
  # 7. Calculate the Hat Matrix (H) for leverage/predictions
  H <- Z_reduced %*% solve(crossprod(Z_reduced) + lam * S_reduced) %*% t(Z_reduced)
  
  return(list(H = H, 
              B = B, 
              S = S, 
              W = S_reduced, 
              Z = Z_reduced, 
              VA = VA, 
              var_prop = variance_proportion,
              fit_coef = pci_coef, 
              lambda = lam, 
              eigen_values = eigen_values))
}

# Function for generating data
gendat2 <- function(N, P, noise_sd, rho = 2) {
  # Time grid and quadrature step size
  grid <- seq(0.01, 1, length.out = P)
  dt <- 1 / P
  
  # True parameter function (true_w)
  true_w <- sin(2 * pi * grid)
  
  # Covariance matrix for X(t)
  dist_mat <- as.matrix(dist(grid))
  Cov_mat <- exp(-dist_mat / rho)
  
  # Generate predictor matrix X
  X <- rmvnorm(n = N, mean = rep(0, P), sigma = Cov_mat)
  
  # Compute continuous signal with step size dt
  signal <- as.vector((X %*% true_w) * dt)
  
  # Generate response y
  y <- signal + rnorm(N, mean = 0, sd = noise_sd)
  
  return(list(
    y = y,
    X = X * dt,
    grid = grid,
    true_w = true_w,
    signal = signal
  ))
}

## REML
REML <- function(mod, y) {
  lam <- mod$lambda
  Z <- mod$Z
  W <- mod$W
  n <- dim(Z)[1]
  V_lam <- diag(n) + Z %*% solve(W) %*% t(Z) / lam
  log_score <- -1/2 * (log(det(V_lam)) + log(sum(solve(V_lam))) + 
                         (length(y) - 1) * log(crossprod(y, solve(V_lam) %*% y)))
  return(log_score)
}

## Leave-one-out cross-validation (LOO)
loo_withse <- function(X, Y, H) {
  Y_hat <- H %*% Y
  pred_diff <- Y - Y_hat
  loo_score <- 0
  n <- length(Y)
  for (i in 1:n) {
    loo_score <- append(loo_score, (pred_diff[i] / (1 - H[i, i]))^2)
  }
  return(list(score = mean(loo_score), 
              se = sd(loo_score)))
}

## Generalized cross-validation (GCV)
GCV <- function(X, Y, H) {
  Y_hat <- H %*% Y
  pred_diff <- Y - Y_hat
  loo_score <- 0
  n <- length(Y)
  hii <- mean(length(H[1, ]) -  sum(diag(H)))
  for (i in 1:n) {
    loo_score <- append(loo_score, (pred_diff[i] / hii)^2)
  }
  return(mean(loo_score))
}

# Function for computing true prediction error of FPCR model
comp_generr2 <- function(B_matrix, fitted_coef, VA) {
  newdat <- gendat2(N = 100, P = 100, noise_sd = 0.1)
  newX <- newdat$X
  newY <- newdat$y
    
  pred_Y <- newX %*% B_matrix %*% VA %*% fitted_coef
  gen_err <- mean((pred_Y - newY)^2)
  return(gen_err)
}

# Main function for simulation
all_in_one <- function(dat, pcis, lambdas) {
  X <- dat$X
  y <- dat$y
  lam_length <- length(lambdas)
  ## Do NCV pe + mse, GCV, CV and REML for all combinations of pci and lambda
  cat("Start computing scores and standard errors...\n")
  ncv_pes <- data.frame(matrix(NA, ncol = length(pcis), nrow = lam_length))
  ncv_pes_debias <- data.frame(matrix(NA, ncol = length(pcis), nrow = lam_length))
  ncv_mses <- data.frame(matrix(NA, ncol = length(pcis), nrow = lam_length))
  gcv_pes <- data.frame(matrix(NA, ncol = length(pcis), nrow = lam_length))
  cv_pes <- data.frame(matrix(NA, ncol = length(pcis), nrow = lam_length))
  loo_se <- data.frame(matrix(NA, ncol = length(pcis), nrow = lam_length))
  reml_scores <- data.frame(matrix(NA, ncol = length(pcis), nrow = lam_length))
  var_prop <- c()
  for (num in pcis) {
    cat("Running for component", num, "\n")
    one_pe <- c()
    one_pe_debias <- c()
    one_mse <- c()
    one_gcv <- c()
    one_cv <- c()
    one_loo_se <- c()
    one_reml <- c()
    cat("Start computing quantities for all tuning parameter lambdas...\n")
    for (i in seq_along(lambdas)) {
      lambda <- lambdas[i]
      fpcr_mod <- fit_fpcr(X, y, pci = num, lam = lambda)
      hat_mat <- fpcr_mod$H
      one_gcv <- c(one_gcv, GCV(X, y, hat_mat))
      loo_output <- loo_withse(X, y, hat_mat)
      one_cv <- c(one_cv, loo_output$score)
      one_loo_se <- c(one_loo_se, loo_output$se)
      
      one_reml <- c(one_reml, REML(mod = fpcr_mod, y = y))
      ncv_result <- fastncv_pal(X, y, fpcr_mod$Z, fpcr_mod$W, K = length(y) / 5, loss = sqloss, 
                                lambda = lambda, R = r)
      one_pe <- c(one_pe, ncv_result$ncv_pe)
      one_pe_debias <- c(one_pe_debias, ncv_result$ncv_pe - ncv_result$bias_correct)
      one_mse <- c(one_mse, ncv_result$ncv_predint)
      
      ## report progress
      pct <- round((i / lam_length) * 100)
      cat("\rProgress: ", pct, "% finished, ", (100 - pct), "% left to do.", sep="")
      flush.console()
    }
    cat("\n\n")
    
    var_prop <- c(var_prop, fpcr_mod$var_prop)
    reml_scores[, which(pcis == num)] <- one_reml
    cv_pes[, which(pcis == num)] <- one_cv
    loo_se[, which(pcis == num)] <- one_loo_se
    gcv_pes[, which(pcis == num)] <- one_gcv
    ncv_pes[, which(pcis == num)] <- one_pe
    ncv_pes_debias[, which(pcis == num)] <- one_pe_debias
    ncv_mses[, which(pcis == num)] <- one_mse
  }
  
  ## compute the prediction errors
  cat("Start computing prediction errors...\n")
  pred_list <- list()
  for (num in pcis) {
      cat("Running for component", num, "\n")
      reps <- 1000
      gen_pes <- data.frame(matrix(NA, ncol = length(lambdas), nrow = reps))
      cat("Start computing prediction errors for each tuning parameter...\n")
      for (j in seq_along(lambdas)) {
        lambda <- lambdas[j]
        pes <- c()

        fpcr_mod <- fit_fpcr(X, y, pci=num, lam = lambda)
        fit_coef <- fpcr_mod$fit_coef
        B <- fpcr_mod$B
        VA <- fpcr_mod$VA

        pes <- mclapply(1:reps, function(i) {comp_generr2(B, fit_coef, VA)})
        gen_pes[, j] <- unlist(pes)

        ## report progress
        pct <- round((j / lam_length) * 100)
        cat("\rProgress: ", pct, "% finished, ", (100 - pct), "% left to do.", sep="")
        flush.console()
      }
      cat("\n\n")

      pred_list[[which(num == pcis)]] <- gen_pes
  }
  cat("Simulation is done. \n")
  return(list(reml = reml_scores, 
              cv = cv_pes, 
              loo_se = loo_se,
              ncv_pes = ncv_pes, 
              ncv_pes_debias = ncv_pes_debias,
              ncv_mses = ncv_mses, 
              gcv = gcv_pes, 
              preds = pred_list, 
              var_prop = var_prop))
}


## Simple test cases for the code
# dat <- gendat2(100, 1)
# X <- dat$X
# y <- dat$y
# mod <- fit_fpcr(X, y, pci = 5, lam = 1)
# output <- fastncv_pal(X, y, mod$Z, mod$W, K = 5, loss = sqloss, lambda = 1, R = 100)

# Generate Data
set.seed(20894177)
alldats2 <- list()
for (j in 1:10) {
   newdat <- gendat2(100, 100, 0.1)
   alldats2[[j]] <- newdat
}

# Set up parameters for the main simulation
pcis <- c(5, 10)
lambdas <- 10^seq(-1, 2.2, length.out = 100)
data_index <- c(2) # datasets to run the main simulation on

## Do simulations
tm <- Sys.time()
all_results <- vector("list", length = length(data_index))
for (index in 1:length(data_index)) {
  cat("Running simulation for dataset", data_index[index], "\n")
  all_results[[index]] <- all_in_one(alldats2[[data_index[index]]], pcis, lambdas)
}
simtime <- as.numeric(difftime(Sys.time(), tm, units='secs'))
cat("Finished simulations, they took", simtime, "seconds.\n")
cat("Saving simulations...\n")
save(all_results, pcis, lambdas, data_index, file = file.path(resultspath, simresultsname)) # Save the simulation results to a file
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")

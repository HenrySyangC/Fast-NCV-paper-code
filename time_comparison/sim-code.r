## Simulations for Fast NCV
## Shuyang (Henry) Cao
## June 2026
## File 01a: Code for simulation of computation comparison

if (!exists("FAST_DEMO")) FAST_DEMO <- FALSE  
# Define parameters dynamically based on mode
if (FAST_DEMO) {
  message("--> RUNNING IN FAST_DEMO MODE (Toy Example)")
  repos <- 5
} else {
  message("--> RUNNING IN PRODUCTION MODE (Full Manuscript Runs)")
  repos <- 100
}

simname <- "sims-time-comparison"
cat("Number of replications for benchmarking computation time is chosen as", repos, ".\n")

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

## Define helper functions
#### loss function -- square loss
sqloss <- function(y1, y2) {
    return((y1-y2)^2)
}

#### Main Functions -- Two ways of running NCV
#### fastncv uses the proposed formula; ncv does not
fastncv <- function(X, y, model_mat, penalty_mat, K, loss=sqloss, lambda = 1) {
    # Randomly assign K folds
    n <- length(y)
    nk <- n / K
    idx <- sample.int(n)
    kidx <- rep(1:K,each=nk)
    foldidx <- split(idx,kidx) 
    # Construct the penalty matrix and store it for later usage
    store_mat <- solve(crossprod(model_mat) + lambda * penalty_mat)
    # compute residuals of the full model
    Res <- y - model_mat %*% (store_mat %*% crossprod(model_mat, y))
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

## Original NCV function without the fast formula
ncv <- function(X, y, model_mat, penalty_mat, K, loss=sqloss, lambda = 1) {
    # Randomly assign K folds
    n <- length(y)
    nk <- n / K
    idx <- sample.int(n)
    kidx <- rep(1:K,each=nk)
    foldidx <- split(idx,kidx) 
    # Fit the full model
    H <- tcrossprod(model_mat %*% solve(crossprod(model_mat) + lambda * penalty_mat), model_mat)
    # compute residuals of the full model
    Res <- y - H %*% y
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

                sub_model_mat <- model_mat[-temp_idx, ]
                sub_y <- y[-temp_idx]

                beta <- solve(crossprod(sub_model_mat) + lambda * penalty_mat, crossprod(sub_model_mat, sub_y))
                res_array[i, i, ] <- y[temp_idx] - model_mat[temp_idx, ] %*% beta
            } else {
                temp_idx <- c(foldidx[[i]], foldidx[[j]])
                sub_model_mat <- model_mat[-temp_idx, ]
                sub_y <- y[-temp_idx]

                beta <- solve(crossprod(sub_model_mat) + lambda * penalty_mat, crossprod(sub_model_mat, sub_y))
                temp_res <- y[temp_idx] - model_mat[temp_idx, ] %*% beta
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

#### fit_models is the function that conducts the time comparsion simulation
fit_models <- function(n, p, K, tms) {
    # n: sample size
    # p: number of co-variates
    # K: number of folds for NCV
    # tms: number of times to run each computation, for timing
    # Simulate data
    x <- sort(runif(n, 0, 1))
    y <- sin(2 * pi * x) + rnorm(n, sd = 0.2)
    y <- y - mean(y)

    # Use smoothCon to generate the basis and penalty
    # 's(x)' specifies a thin plate regression spline by default
    sm <- smoothCon(s(x, k = p), data = data.frame(x = x))[[1]]

    X <- sm$X     # The Basis Matrix
    S <- sm$S[[1]]   # The Penalty Matrix
    ## Method 0: Slow Nested CV with Penalized Regression
    method0time <- microbenchmark::microbenchmark({
    ncv(x, y, X, S, K)
    }, times = tms)
    
    ## Method 1: Fast Nested CV with Linear Regression
    method1time <- microbenchmark::microbenchmark({
    fastncv(x, y, X, S, K)
    }, times = tms)
    
    # Clean the results
    method0time <- data.frame(n = n, p = p, scaling = n / K, method = "slow", time = method0time$time)
    method1time <- data.frame(n = n, p = p, scaling = n / K, method = "fast", time = method1time$time)

    rbind(method0time, method1time) 
}

## Set up simulation details
ntodo <- c(100, 200, 400, 800, 1600)
s <- 5 # scaling factor for K, K = n / s
ptodo <- c(5, 10, 20, 40)
resultslist <- list()

## Do simulations
tm <- Sys.time()
#### Regime 1: p is fixed, n varies
idx <- 1
cat("Running simulations for regime 1: p is fixed. \n")
for (n in ntodo){
    for (p in ptodo) {
    cat("Running simulations for n =", n, "and p =", p, "...\n")
    resultslist[[idx]] <- fit_models(n = n, p = p, K = n / s, tms = repos)
    idx <- idx + 1
    }
}

#### Regime 2: p scales with n
ntodo <- c(100, 200, 400, 800, 1600, 3200) # update ntodo for regime 2
indices <- 1:length(ptodo)
n_indices <- indices
groups <- c(1, 2, 3) # group indices for p
cat("Running simulations for regime 2: p scales with n. \n")
for (index in groups) {
    for (i in indices) {
        cat("Running simulations for n =", ntodo[n_indices[i]], "and p =", ptodo[i], "...\n")
        idx <- idx + 1
        resultslist[[idx]] <- fit_models(n = ntodo[n_indices[i]], p = ptodo[i], K = ntodo[n_indices[i]] / s, tms = repos)
        }
    n_indices <- n_indices + 1
}
simtime <- as.numeric(difftime(Sys.time(), tm, units='secs'))
cat("Finished simulations, they took", simtime, "seconds.\n")
cat("Saving simulations...\n")
save(resultslist, file = file.path(resultspath, simresultsname)) # Save the simulation results to a file
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
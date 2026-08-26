## Simulations for Fast NCV
## Shuyang (Henry) Cao
## June 2026
## File 02a: Code for simulation of computation cutoff point for Fast NCV 

if (!exists("FAST_DEMO")) FAST_DEMO <- FALSE  
# Define parameters dynamically based on mode
if (FAST_DEMO) {
  message("--> RUNNING IN FAST_DEMO MODE (Toy Example)")
  repos <- 3 # Number of replications for computation of each (n, p) pair
} else {
  message("--> RUNNING IN PRODUCTION MODE (Full Manuscript Runs)")
  repos <- 10 # Number of replications for computation of each (n, p) pair
}

simname <- "sims-comp-cutoff"
cat("Running", repos, "replications for benchmarking computation time of each (n, p) pair.\n")

## import packages ##
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

## Set paths ##
basepath <- getwd()
resultspath <- file.path(basepath, "results")
if (!dir.exists(resultspath)) dir.create(resultspath)
figurespath <- file.path(basepath, "figures")
if (!dir.exists(figurespath)) dir.create(figurespath)
simresultsname <- paste0(simname, ".RData")

## define helper functions
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
    # compute NCV prediction interval
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

#### original ncv implementation that does not use the proposed formula
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

find_cutoff <- function(n, p) {
    A <- 6*p^2*n + p^3 + 15*p^2 + 6*p*n
    B <- -18*p^2*n + 6*n + p^3 + 15*p^2 + 6*p*n
    C <-  -24*p*n^2 - 60*n^2
    D <- -8*n^3 + 12*p*n^2 + 30*n^2
    E <- 6*n^3
    f <- function(K) {
        (A*K^4) + (B*K^3) + (C*K^2) + (D*K) + E
    }
    result <- uniroot(f, interval = c(2, n))
    return(result$root)
}


## Function that computes the computation cutoff point
## and record the computation time for comparison
comp_cutoff <- function(n, p) {
    x <- sort(runif(n, 0, 1))
    y <- sin(2 * pi * x) + rnorm(n, sd = 1)
    y <- y - mean(y)

    sm <- smoothCon(s(x, k = p), data = data.frame(x = x))[[1]]

    X <- sm$X     # The Basis Matrix
    S <- sm$S[[1]]   # The Penalty Matrix (List because there can be multiple)

    cutoff_K <- find_cutoff(n = n, p = p)

    method0time_a <- microbenchmark::microbenchmark({
        ncv(x, y, X, S, K = 5, lambda = 2)
    }, times = repos)
    method0time_b <- microbenchmark::microbenchmark({
        ncv(x, y, X, S, K = 10, lambda = 2)
    }, times = repos)
    method0time_c <- microbenchmark::microbenchmark({
        ncv(x, y, X, S, K = 20, lambda = 2)
    }, times = repos)
    method0time_d <- microbenchmark::microbenchmark({
        ncv(x, y, X, S, K = 50, lambda = 2)
    }, times = repos)
    method0time_e <- microbenchmark::microbenchmark({
        ncv(x, y, X, S, K = 100, lambda = 2)
    }, times = repos)

    method1time_a <- microbenchmark::microbenchmark({
        fastncv(x, y, X, S, K = 5, lambda = 2)
    }, times = repos)
    method1time_b <- microbenchmark::microbenchmark({
        fastncv(x, y, X, S, K = 10, lambda = 2)
    }, times = repos)
    method1time_c <- microbenchmark::microbenchmark({
        fastncv(x, y, X, S, K = 20, lambda = 2)
    }, times = repos)
    method1time_d <- microbenchmark::microbenchmark({
        fastncv(x, y, X, S, K = 50, lambda = 2)
    }, times = repos)
    method1time_e <- microbenchmark::microbenchmark({
        fastncv(x, y, X, S, K = 100, lambda = 2)
    }, times = repos)

    result <- c(mean(round(method0time_a$time / method1time_a$time, 3)), 
                mean(round(method0time_b$time / method1time_b$time, 3)), 
                mean(round(method0time_c$time / method1time_c$time, 3)), 
                mean(round(method0time_d$time / method1time_d$time, 3)), 
                mean(round(method0time_e$time / method1time_e$time, 3)))
    return(list(result = result, cutoff = cutoff_K))
}

# Set simulation details
## (n, p) pairs for simulation
## 9 groups of (n, p) with their predicted cutoff point from f(K)=0
K_vals <- c(5, 10, 20, 50,  100)
ntodo <- c(1000, 1000, 1000, 2000, 2000, 2000, 4000, 4000, 4000)
ptodo <- c(10, 20, 50, 10, 20, 50, 10, 20, 50)

# Do simulation
tm <- Sys.time()
resultslist <- list(); length(resultslist) <- length(ntodo)
for (i in 1:length(ntodo)) {
    cat("Running simulation for (n, p) = (", ntodo[i], ", ", ptodo[i], ")...\n")
    resultslist[[i]] <- comp_cutoff(n = ntodo[i], p = ptodo[i])
}
simtime <- as.numeric(difftime(Sys.time(), tm, units='secs'))
cat("Finished simulations, they took", simtime, "seconds.\n")
cat("Saving simulations...\n")
save(resultslist, file = file.path(resultspath, simresultsname)) # Save the simulation results to a file
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
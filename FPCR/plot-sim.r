## Simulations for Fast NCV
## Shuyang (Henry) Cao
## June 2026
## File 03b: Code for plot simulation results of FPCR for Fast NCV 

simname <- "sims-FPCR"

pkgs <- c(
    'tidyverse', 'scales', 'ggplot2' # For summarizing results
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

## Load the simulation results
load(file.path(resultspath, paste0(simname, ".RData")))

# Data Processing
# Define vectors for datasets and components
num_datasets <- length(data_index)
num_components <- length(pcis)
full_range <- 1:length(lambdas)      # Full data index range for finding true optimums


# Data Visualization CONFIGURATION
## Choose the range of lambda values to visualize
start_visual_idx <- 1  # Change this to shift the start window (e.g., 30, 50)
end_visual_idx   <- 100 # Change this to shift the end window (e.g., 80, 90)

pred_data_list   <- list()
method_data_list <- list()
opt_lines_list   <- list()
reml_data_list   <- list()

# Track unique combinations to order factor levels later
unique_variance_labels <- c()

# Nested loop: Loop through each dataset, then each component
for (k in 1:num_datasets) {
  one_data <- all_results[[k]]
  
  reml       <- one_data$reml
  loo        <- one_data$cv
  ncv_debias <- one_data$ncv_pes_debias
  gcv        <- one_data$gcv
  preds      <- one_data$preds
  ncvmses    <- one_data$ncv_mses
  var_prop   <- one_data$var_prop
  
  for (pind in 1:num_components) {
    ds_label  <- paste("Dataset", k)
    
    # FIX 1: Include component index to guarantee unique labels even if percentages match
    var_label <- paste0("Comp ", pind, ": ", round(var_prop[pind] * 100, 2), "% Variance Explained")
    
    if (k == 1) {
      unique_variance_labels <- c(unique_variance_labels, var_label)
    }
    
    # 1. Prepare prediction error data (Full Range)
    pred_data_list[[length(pred_data_list) + 1]] <- data.frame(
      dataset    = ds_label,
      component  = var_label,  
      lambda_idx = rep(full_range, each = nrow(preds[[pind]])),
      pred_error = c(as.matrix(preds[[pind]][, full_range]))
    )
    
    # 2. Prepare method estimates and intervals (Full Range)
    method_data_list[[length(method_data_list) + 1]] <- data.frame(
      dataset    = ds_label,
      component  = var_label,  
      lambda_idx = full_range,
      ncv_pe     = ncv_debias[full_range, pind],
      ncv_lower  = ncv_debias[full_range, pind] - 1.96 * sqrt(ncvmses[full_range, pind]),
      ncv_upper  = ncv_debias[full_range, pind] + 1.96 * sqrt(ncvmses[full_range, pind])
    )
    
    # 3. Find global optimal lambda indices across the entire range
    opt_lines_list[[length(opt_lines_list) + 1]] <- data.frame(
      dataset   = ds_label,
      component = var_label,   
      reml_opt  = which.max(reml[, pind]),
      loo_opt   = which.min(loo[, pind]),
      gcv_opt   = which.min(gcv[, pind])
    )

    reml_data_list[[length(reml_data_list) + 1]] <- data.frame(
      dataset    = ds_label,
      component  = var_label,
      lambda_idx = full_range,
      reml_val   = -reml[full_range, pind] # Flipped curve
    )
  }
}

# Bind everything into master data frames
pred_data      <- do.call(rbind, pred_data_list)
method_data    <- do.call(rbind, method_data_list)
opt_lines_data <- do.call(rbind, opt_lines_list)
reml_data      <- do.call(rbind, reml_data_list)

# Filter out indices outside global visualization threshold window
pred_data   <- pred_data[pred_data$lambda_idx >= start_visual_idx & pred_data$lambda_idx <= end_visual_idx, ]
method_data <- method_data[method_data$lambda_idx >= start_visual_idx & method_data$lambda_idx <= end_visual_idx, ]
reml_data   <- reml_data[reml_data$lambda_idx >= start_visual_idx & reml_data$lambda_idx <= end_visual_idx, ]

# FIX 2: Force unique() on levels to prevent factor duplicate error
unique_levels <- unique(unique_variance_labels)
pred_data$component      <- factor(pred_data$component, levels = unique_levels)
method_data$component    <- factor(method_data$component, levels = unique_levels)
opt_lines_data$component <- factor(opt_lines_data$component, levels = unique_levels)
reml_data$component      <- factor(reml_data$component, levels = unique_levels)

# Helper function for dynamic X-axis
get_lambda_axis_spec <- function(lambdas, n_breaks = 8) {
  log_range <- log10(range(lambdas))
  
  exp_min <- ceiling(log_range[1])
  exp_max <- floor(log_range[2])
  exponents <- seq(exp_min, exp_max, by = 1)
  
  if (length(exponents) > n_breaks) {
    step <- ceiling(length(exponents) / n_breaks)
    exponents <- exponents[seq(1, length(exponents), by = step)]
  }
  
  target_lambdas <- 10^exponents
  idx_breaks <- sapply(target_lambdas, function(val) which.min(abs(lambdas - val)))
  idx_labels <- lapply(exponents, function(e) bquote(10^.(e)))
  
  list(breaks = idx_breaks, labels = idx_labels)
}

idx <- get_lambda_axis_spec(lambdas)
idx_breaks <- idx$breaks
idx_labels <- idx$labels

p1 <- ggplot(pred_data, aes(x = lambda_idx, y = pred_error)) +
  
  # Boxplot layer
  geom_boxplot(
    aes(group = lambda_idx, fill = "Generalization Error"), 
    alpha = 0.7, width = 0.5, outlier.shape = NA, linewidth = 0.3
  ) +
  
  # NCV Upper and Lower bounds
  geom_line(
    data = method_data, 
    aes(x = lambda_idx, y = ncv_upper, linetype = "NCV 95% Upper/Lower", group = 1),
    color = "gray40", linewidth = 0.4, inherit.aes = FALSE
  ) +
  geom_line(
    data = method_data, 
    aes(x = lambda_idx, y = ncv_lower, linetype = "NCV 95% Upper/Lower", group = 1),
    color = "gray40", linewidth = 0.4, inherit.aes = FALSE
  ) +
  
  # Point Estimate & Optimal Lines
  geom_line(
    data = method_data, 
    aes(x = lambda_idx, y = ncv_pe, linetype = "NCV Point Est.", group = 1),
    color = "black", linewidth = 0.8, inherit.aes = FALSE
  ) +
  geom_vline(
    data = opt_lines_data, 
    aes(xintercept = reml_opt, linetype = "REML Opt."), 
    color = "black", linewidth = 0.5
  ) +
  geom_vline(
    data = opt_lines_data, 
    aes(xintercept = loo_opt, linetype = "LOO Opt."), 
    color = "black", linewidth = 0.5
  ) +
  geom_vline(
    data = opt_lines_data, 
    aes(xintercept = gcv_opt, linetype = "GCV Opt."), 
    color = "black", linewidth = 0.5
  ) +
  
  facet_wrap(~ component, ncol = 1) +  

  # Manual scales
  scale_fill_manual(
    name = NULL, 
    values = c("Generalization Error" = "gray93")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "NCV 95% Upper/Lower" = "solid",
      "NCV Point Est."      = "solid",      
      "REML Opt."           = "dashed",
      "LOO Opt."            = "longdash",
      "GCV Opt."            = "twodash"
    )
  ) +
  
  # Clean mathematical power-of-10 x-axis mapping
  scale_x_continuous(
    breaks = idx_breaks,
    labels = idx_labels
  ) + 

  # FIX 3: Tight Y-axis scaling using standard R functions without over-expanding
  scale_y_continuous(
    expand = expansion(mult = c(0.01, 0.02)),
    breaks = function(limits) pretty(limits, n = 5)
  ) +

  coord_cartesian(xlim = c(1, 100)) +
  
  labs(
    title = "Generalization Error Analysis",
    x = expression(Tuning ~ Parameter ~ (lambda)),
    y = "Generalization Error"
  ) +
  
  # Journal Base Theme
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 11, color = "black", hjust = 0.5, margin = margin(b = 6)),
    axis.title.x     = element_text(face = "bold", size = 9.5, color = "black", margin = margin(t = 4)),
    axis.title.y     = element_text(face = "bold", size = 9.5, color = "black", margin = margin(r = 4)),
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 8, color = "black"),
    axis.text.y      = element_text(size = 8, color = "black"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 8.5, color = "black"), 
    
    # Legend settings
    legend.position   = "bottom",
    legend.box        = "vertical",
    legend.box.just   = "center",
    legend.text       = element_text(size = 6.8),
    legend.margin     = margin(t = -2, b = 0),
    legend.box.margin = margin(t = -2, b = 0),
    legend.key.width  = unit(0.55, "lines"),
    legend.key.height = unit(0.75, "lines"),
    legend.spacing.x  = unit(0.12, "lines"),
    
    panel.spacing.y   = unit(0.4, "lines")
  ) +
  
  guides(
    fill = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2, nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0.6))
  )

p2 <- ggplot(reml_data, aes(x = lambda_idx, y = reml_val)) +
  
  # REML Curve layer
  geom_line(
    aes(group = dataset, color = "Log-REML Curve", linetype = "Log-REML Curve"), 
    linewidth = 0.6
  ) +
  
  # Vertical Optimal Lines
  geom_vline(
    data = opt_lines_data, 
    aes(xintercept = reml_opt, color = "REML Opt.", linetype = "REML Opt."), 
    linewidth = 0.5
  ) +
  geom_vline(
    data = opt_lines_data, 
    aes(xintercept = loo_opt, color = "LOO Opt.", linetype = "LOO Opt."), 
    linewidth = 0.5
  ) +
  geom_vline(
    data = opt_lines_data, 
    aes(xintercept = gcv_opt, color = "GCV Opt.", linetype = "GCV Opt."), 
    linewidth = 0.5
  ) +
  
  # Scale Mappings
  scale_color_manual(
    name = NULL,
    values = c(
      "Log-REML Curve" = "black",
      "REML Opt."      = "black",
      "LOO Opt."       = "black",
      "GCV Opt."       = "black"
    )
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Log-REML Curve" = "solid",
      "REML Opt."      = "dashed",
      "LOO Opt."       = "longdash",
      "GCV Opt."       = "twodash"
    )
  ) +
  
  facet_wrap(~ component, ncol = 1) +
  
  # Clean mathematical power-of-10 x-axis mapping
  scale_x_continuous(
    breaks = idx_breaks,
    labels = idx_labels
  ) + 
  coord_cartesian(xlim = c(1, 100)) +
  
  labs(
    title = "Negative Log-REML Curves",
    x = expression(Tuning ~ Parameter ~ (lambda)),
    y = "Negative Log-REML Score"
  ) +
  
  
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 11, color = "black", hjust = 0.5, margin = margin(b = 6)),
    axis.title.x     = element_text(face = "bold", size = 9.5, color = "black", margin = margin(t = 4)),
    axis.title.y     = element_text(face = "bold", size = 9.5, color = "black", margin = margin(r = 4)),
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 8, color = "black"), # Horizontal alignment for power labels
    axis.text.y      = element_text(size = 8, color = "black"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 8.5, color = "black"), 
    
    # Legend settings
    legend.position   = "bottom",
    legend.box        = "vertical",
    legend.box.just   = "center",
    legend.text       = element_text(size = 6.8),
    legend.margin     = margin(t = -2, b = 0),
    legend.box.margin = margin(t = -2, b = 0),
    legend.key.width  = unit(0.55, "lines"),
    legend.key.height = unit(0.75, "lines"),
    legend.spacing.x  = unit(0.12, "lines"),
    
    panel.spacing.y   = unit(0.4, "lines")
  ) +
  
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0.6)),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

# Save the plots as a PDF
ggsave(
  filename = file.path(figurespath, paste0(simname, "_pred_error.pdf")),
  plot = p1, 
  width = 3.3, 
  height = 5.0, 
  units = "in", 
  useDingbats = FALSE
)
cat("Saved plot to file:",file.path(figurespath, paste0(simname, "_pred_error.pdf")), "\n")

ggsave(
  filename = file.path(figurespath, paste0(simname, "_reml.pdf")),
  plot = p2, 
  width = 3.3, 
  height = 5.0, 
  units = "in", 
  useDingbats = FALSE
)
cat("Saved plot to file:",file.path(figurespath, paste0(simname, "_reml.pdf")), "\n")
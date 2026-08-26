## Simulations for Fast NCV
## Shuyang (Henry) Cao
## June 2026
## File 01b: Code for simulation of computation comparison

simname <- "sims-time-comparison"

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

## Data processing
results <- bind_rows(resultslist)
results$group <- paste0("(n = ", results$n, ", p = ", results$p, ")")

# Prepare ratio data with n and p preserved
ratio_data_full <- results %>%
    group_by(n, p) %>%
    reframe(
        ratio = time[method == "slow"] / time[method == "fast"],
        .groups = "drop" # Drop the groups here
    ) %>%
    # Now this filters outliers relative to the WHOLE dataset
    filter(
        ratio >= quantile(ratio, 0.025, na.rm = TRUE) & 
        ratio <= quantile(ratio, 0.975, na.rm = TRUE)
    ) %>%
    as.data.frame()

# Create p scaling groups for regime 2
p_scale_groups <- data.frame(
    n = c(100, 200, 400, 800,  200, 400, 800, 1600, 400, 800, 1600, 3200),
    p = c(5, 10, 20, 40, 5, 10, 20, 40, 5, 10, 20, 40),
    p_scale_group = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
)

log2_even_breaks <- function(x) {
  rng <- range(x, na.rm = TRUE)
  min_pow <- floor(log2(rng[1]))
  max_pow <- ceiling(log2(rng[2]))
  2^(min_pow:max_pow)
}
# ============================================================
# Regime 1: n increases, p fixed
# ============================================================
ratio_regime1 <- ratio_data_full %>%
    mutate(p = as.factor(p))

p1 <- ggplot(ratio_regime1, aes(x = factor(n), y = ratio)) +
  
  geom_boxplot(
    linewidth = 0.4, 
    outlier.size = 1, 
    outlier.alpha = 0.6,
    fill = "gray97"
  ) +
  
  # Forced even 2^0, 2^1, 2^2, 2^3... breaks
  scale_y_continuous(
    trans = "log2",
    breaks = log2_even_breaks,
    labels = scales::label_log(base = 2)
  ) + 
  
  facet_wrap(
    ~ p, 
    nrow = 1, 
    labeller = labeller(p = c("5" = "p = 5", "10" = "p = 10", "20" = "p = 20", "40" = "p = 40"))
  ) +
  
  labs(
    title = "NCV Time Comparison: n Increases, p Fixed",
    x = "Sample Size (n)",
    y = expression(Log[2] ~ "Relative Computation Time (Slow/Fast)")
  ) +
  
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 11, color = "black", hjust = 0.5, margin = margin(b = 6)),
    axis.title.x     = element_text(face = "bold", size = 9.5, color = "black", margin = margin(t = 4)),
    axis.title.y     = element_text(face = "bold", size = 9.5, color = "black", margin = margin(r = 4)),
    axis.text.x      = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 8, color = "black"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 8.5, color = "black"),
    
    panel.spacing.x  = unit(0.5, "lines")
  )

# ============================================================
# Regime 2: p scales with n
# ============================================================
ratio_regime2 <- ratio_data_full %>%
    inner_join(p_scale_groups, by = c("n", "p")) %>%
    mutate(p_col = as.character(p_scale_group))

p2 <- ggplot(ratio_regime2, aes(x = factor(n), y = ratio)) +
  
  geom_boxplot(
    linewidth = 0.4, 
    outlier.size = 1, 
    outlier.alpha = 0.6,
    fill = "gray97"
  ) +
  
  # Forced even 2^0, 2^1, 2^2, 2^3... breaks
  scale_y_continuous(
    trans = "log2",
    breaks = log2_even_breaks,
    labels = scales::label_log(base = 2)
  ) +
  
  facet_wrap(
    ~ p_col, 
    nrow = 1, 
    scales = "free_x", 
    labeller = labeller(p_col = c(
      "1" = "n/p = 20",
      "2" = "n/p = 40",
      "3" = "n/p = 80"
    ))
  ) +
  
  labs(
    title = "NCV Time Comparison: p Scales with n",
    x = "Sample Size (n)",
    y = expression(Log[2] ~ "Relative Computation Time (Slow/Fast)")
  ) +
  
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 11, color = "black", hjust = 0.5, margin = margin(b = 6)),
    axis.title.x     = element_text(face = "bold", size = 9.5, color = "black", margin = margin(t = 4)),
    axis.title.y     = element_text(face = "bold", size = 9.5, color = "black", margin = margin(r = 4)),
    axis.text.x      = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 8, color = "black"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 8.5, color = "black"),
    
    panel.spacing.x  = unit(0.5, "lines")
  )

## Save the plot as a PDF
ggsave(
  filename = file.path(figurespath, paste0(simname, "_pfixed.pdf")),
  plot = p1, 
  width = 6.85, # Full page width (174 mm)
  height = 3.5, # Compact height appropriate for a single row of panels
  units = "in",
  useDingbats = FALSE
)
cat("Saved plot to file:",file.path(figurespath, paste0(simname, "_pfixed.pdf")), "\n")

ggsave(
  filename = file.path(figurespath, paste0(simname, "_pscaled.pdf")),
  plot = p2, 
  width = 6.85, 
  height = 3.5, 
  units = "in",
  useDingbats = FALSE
)
cat("Saved plot to file:",file.path(figurespath, paste0(simname, "_pscaled.pdf")), "\n")

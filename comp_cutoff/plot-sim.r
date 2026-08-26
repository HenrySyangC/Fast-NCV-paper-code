## Simulations for Fast NCV
## Shuyang (Henry) Cao
## June 2026
## File 02b: Code for plotting computation cutoff results for Fast NCV

simname <- "sims-comp-cutoff"

pkgs <- c(
    'tidyverse', 'purrr', 'dplyr', 'scales', 'ggplot2', 
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

## Load the simulation results
load(file.path(resultspath, paste0(simname, ".RData")))

## Paramters from simulation to help plotting
K_vals <- c(5, 10, 20, 50, 100)
ntodo <- c(1000, 1000, 1000, 2000, 2000, 2000, 4000, 4000, 4000)
ptodo <- c(10, 20, 50, 10, 20, 50, 10, 20, 50)

## Data Processing for Plotting
num_elements <- length(resultslist)

df_long <- map_dfr(seq_len(num_elements), function(i) {
  data.frame(
    K_vals = K_vals,
    Dataset = paste0("Result", i), 
    Relative_Time = resultslist[[i]]$result,
    Cutoff_K = resultslist[[i]]$cutoff
  )
})

unique_datasets <- paste0("Result", seq_along(ntodo))
display_names   <- paste0("(", ntodo, ",", ptodo, ")")
facet_labels    <- setNames(display_names, unique_datasets)

y_max <- max(df_long$Relative_Time, na.rm = TRUE)
y_breaks <- sort(unique(c(1, pretty(c(0, y_max)))))

## Create the plot
p <- ggplot(df_long, aes(x = K_vals, y = Relative_Time, group = Dataset)) +
  
  # LAYER 1: Clean Minimal Background Theme
  theme_minimal() +
  
  # LAYER 2: Separate data into a 3x3 grid based on Dataset (Changed ncol = 3)
  facet_wrap(~ Dataset, ncol = 3, labeller = as_labeller(facet_labels)) +
  
  # LAYER 3: Highly Distinguished ALL-BLACK Baseline ("44" = thick block-dash)
  geom_hline(yintercept = 1, color = "black", linetype = "44", linewidth = 1.2, alpha = 0.8) +
  
  # LAYER 4: Vertical Cutoff Lines
  geom_vline(
    aes(xintercept = Cutoff_K, linetype = "Predicted Comp. Cutoff K"), 
    color = "black", 
    alpha = 0.6, 
    linewidth = 0.8
  ) +
  
  # LAYER 5: Main Data Trajectories 
  geom_line(
    aes(color = "Relative Comp. Time"), 
    linetype = "solid", 
    linewidth = 0.9
  ) +
  geom_point(
    aes(color = "Relative Comp. Time"), 
    shape = 21, 
    fill = "white", 
    size = 2.5
  ) +
  
  # Explicitly injecting 1 into the official Y-axis ticks
  scale_y_continuous(breaks = y_breaks) +
  
  scale_color_manual(
    name = NULL,
    values = c("Relative Comp. Time" = "black")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Predicted Comp. Cutoff K" = "dashed") 
  ) +
  
  # Labels and Formatting
  labs(
    title = "Change of Relative Computation Time by (n, p) Configurations",
    x = "K",
    y = "Relative Computation Time"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, margin = margin(b = 12)),
    legend.position = "bottom", # Moved to bottom to give the 3-column grid more horizontal room
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray94"),
    axis.line.y = element_line(color = "black", linewidth = 0.4),
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),
    
    # Clean facet header strip formatting
    strip.background = element_rect(fill = "gray97", color = "gray88"),
    strip.text = element_text(face = "bold", size = 12, color = "black"), # Adjusted slightly for 3 columns
    
    # Generates explicit spacing padding between the 9 plots
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(1.2, "lines")
  ) +
  
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(shape = 21, fill = "white", linetype = "solid")
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(color = "black") 
    )
  )

p <- ggplot(df_long, aes(x = K_vals, y = Relative_Time, group = Dataset)) +
  
  # Base theme set to 9pt for tight double-column rendering
  theme_minimal(base_size = 9, base_family = "sans") +
  
  facet_wrap(~ Dataset, ncol = 3, labeller = as_labeller(facet_labels)) +
  
  # Reference Lines
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = 1, color = "black", linetype = "44", linewidth = 0.8, alpha = 0.8) +
  
  geom_vline(
    aes(xintercept = Cutoff_K, linetype = "Predicted Comp. Cutoff K"), 
    color = "black", 
    alpha = 0.6, 
    linewidth = 0.6
  ) +
  
  # Data Trajectories
  geom_line(
    aes(color = "Relative Comp. Time"), 
    linetype = "solid", 
    linewidth = 0.6
  ) +
  geom_point(
    aes(color = "Relative Comp. Time"), 
    shape = 21, 
    fill = "white", 
    size = 1.8
  ) +
  
  scale_y_continuous(breaks = y_breaks) +
  
  scale_color_manual(
    name = NULL,
    values = c("Relative Comp. Time" = "black")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Predicted Comp. Cutoff K" = "dashed") 
  ) +
  
  labs(
    title = "Change of Relative Computation Time by (n, p) Configurations",
    x = "K",
    y = "Relative Computation Time"
  ) +
  
  # Springer Journal Typography Guidelines (8pt - 11pt target range)
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 11, margin = margin(b = 6)),
    axis.title.x     = element_text(size = 9.5, margin = margin(t = 4)),
    axis.title.y     = element_text(size = 9.5, margin = margin(r = 4)),
    axis.text        = element_text(size = 8, color = "black"),
    
    legend.position  = "bottom",
    legend.text      = element_text(size = 8.5),
    legend.margin    = margin(t = -2),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    axis.line.y      = element_line(color = "black", linewidth = 0.3),
    
    # Compact facet headers for multi-panel layout
    strip.background = element_rect(fill = "gray95", color = "gray85", linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 8.5, color = "black"), 
    
    panel.spacing.x  = unit(0.6, "lines"),
    panel.spacing.y  = unit(0.6, "lines")
  ) +
  
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(shape = 21, fill = "white", linetype = "solid", size = 1.5)
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(color = "black", linewidth = 0.6) 
    )
  )

## Save the plot as a PDF
ggsave(
  filename = file.path(figurespath, paste0(simname, "_comp_cutoff_plot.pdf")),
  plot = p, 
  width = 6.85,
  height = 6.5,
  units = "in",
  useDingbats = FALSE # Ensures standard font embed compatibility for Springer PDF check
)
cat("Saved plot to file:",file.path(figurespath, paste0(simname, "_comp_cutoff_plot.pdf")), "\n")
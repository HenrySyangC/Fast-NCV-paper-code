# Parse command line arguments (Default to FAST_DEMO = TRUE for safety)
args <- commandArgs(trailingOnly = TRUE)
FAST_DEMO <- if ("--test" %in% args) TRUE else FALSE

if (FAST_DEMO) {
  message("\n=======================================================")
  message("  Rnning TOY EXAMPLE simulations for testing")
  message("  Simulations will finish quickly.")
  message("  To run full paper results, use: Rscript run_all.R --full")
  message("=======================================================\n")
} else {
  message("\n=======================================================")
  message("  Rnning FULL simulations for reproduction")
  message("  Warning: This may take several days to complete.")
  message("=======================================================\n")
}
message("=== Starting Simulations Reproduction ===")

simulations <- c("time_comparison", "comp_cutoff", "FPCR")

for (sim in simulations) {
  message(paste0("\n Processing: ", sim, " ..."))
  source(file.path(sim, "sim-code.R"), chdir = TRUE)
  source(file.path(sim, "plot-sim.R"), chdir = TRUE)
}

message("\n All simulations and plots are completed!")
# Reproduction Code for Fast Computation of Nested Cross-Validation for Penalized Regression

This repository contains the R source code required to reproduce the simulation results and figures presented in:

> **Fast Computation of Nested Cross-Validation for Penalized Regression**  
> Shuyang Cao & Alex Stringer  
> *Preprint (2026)*  
> Link: 

---

## 📂 Repository Structure

```
├── .Rprofile 
├── .gitignore            
├── README.md
├── renv.lock  # Package version blueprint
├── renv/      # System configuration folder for package reproducibility 
│   ├── activate.R
│   └── settings.json
├── run_all.R  # run all simulations with produced plots
├── time_comparison/
│   ├── sim-code.R   # Runs time comparison simulation and saves output to a folder /results
│   └── plot-sim.R   # Generate Figure 1 and Figure 2
├── comp_cutoff/
│   ├── sim-code.R   # Runs simulation for predicting computation cutoff point
│   └── plot-sim.R     # Generate Figure 3
├── FPCR/
│   ├── sim-code.R   # Runs FPCR simulation for tuning parameters
│   └── plot-sim.R     # Generate Figure 4(a) and Figure 4(b)
```

---

## 1. System Requirements & Dependencies

### Software Prerequisites
* **R Version:** `>= 4.2.0`
* **Operating System:** Tested on macOS Sonoma / Ubuntu 22.04 / Windows 11

### Required R Packages
Ensure the following packages are installed before running the scripts:

```R
install.packages(c("ggplot2", "dplyr", "purrr", "scales", "tidyverse", "parallel", 'Matrix', 'refund', 'fda', 'mgcv', 'mvtnorm'))
```

* There is **no need** to install these packages beforehand, since the scripts will check whether they are already installed and install them if not. 
* Optional: If you want version control for packages I used, you can use `renv` with the renv folder I provided. Navigate to the simulation folder (ex: 'sims' in Step 3) and run `renv::restore()` before reproducing any simulation results.

---

## 2. Paper-to-Code Matching
### time_comparison (Section 4--Figure 1 and Figure 2)
+ Compare computation time between the proposed formula and refitting the model
+ Replicate the simulation: `sim-code.R`
+ Visualize the results: `plot-sim.R`

### comp_cutoff (Section 4--Figure 3)
+ Predict the computation cutoff point of K (Proposition 2)
+ Replicate the simulation: `sim-code.R`
+ Visualize the results: `plot-sim.R`

### FPCR (Section 5--Figure 4(a) and Figure 4(b))
+ Tuning parameter selection with plot of generalization error and REML
+ Replicate the simulation: `sim-code.R`
+ Visualize the results: `plot-sim.R`

---

## 3. How to Reproduce the Results

Clone or download this repository into a folder named `sims`:

   ```bash
   git clone https://github.com/HenrySyangC/Fast-NCV.git sims
   ```

Navigate into the folder:

   ```bash
   cd sims
   ```

### Option 1: Running All Simulations Together

Run the script run_all.R when navigating into ~/sims

```bash
Rscript run_all.R
```
```bash
# Simulations results and figures can be found in each individual simulation folder,
# which are created automatically after running the script. 
# Ex: ~/sims/time_comparison/figures
```

### Option 2: Running Individual Simulations

To reproduce a specific simulation and figure, navigate to the respective directory and run the scripts in order:

```bash
# Example for Simulation time_comparison
# For the other two simulations, replace time_comparison with comp_cutoff or FPCR
cd time_comparison

# Step 1: Run simulation and save results to ~/sims/time_comparison/results
# A folder named 'results' will be created automatically
Rscript sim-code.R

# Step 2: Generate plots and save it to ~/sims/time_comparison/figures
# A folder named 'figures' will be created automatically
Rscript plot-sim.R
```

### Estimated Computational Time

| Simulation | Target Figure | Est. Runtime |
| :--- | :--- | :--- |
| **time_comparison** | Figure 1, Figure 2 | ~1.5 days |
| **comp_cutoff** | Figure 3 | ~5 hours | 
| **FPCR** | Figure 4(a), Figure 4(b) | ~10 mins |

* ⚠️ Benchmarked on a **MacBook (M2)**. 
* 💡 *Simulation **FPCR** supports parallel processing—running with more cores will reduce runtime.*

---

## 4. Citation & License

If you use this code in your research, please cite our preprint:

```bibtex
@article{yourname2026title,
  title={Your Paper Title},
  author={Your Name and Co-authors},
  journal={arXiv / BioRxiv / SSRN},
  year={2026}
}
```

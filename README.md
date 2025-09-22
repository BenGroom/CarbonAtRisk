# Carbon Removal Portfolios and Storage Security – Code and Figures

This repository contains the R and Python code used to generate simulations, figures, and analysis for the paper. It covers two strands of analysis:

1. **Geological storage security** (R): Monte Carlo and sensitivity analysis of long-term leakage risk (based on Alcalde et al. 2018, extended here).
2. **Carbon removal portfolio analysis** (Python): Simulation and optimisation of DACCS, Biochar, and Forest portfolios, including Carbon-at-Risk (CaR) illustrations and Markowitz efficient frontier analysis.

All code lives in the `code/` folder. When run, outputs are saved to the `Outputs/` folder.

---

## Repository structure

project-root/
│
├── code/
│ ├── SSC_storage_security.R # Geological storage security model (R)
│ ├── figure1b_truncNormal.py # Main text Fig. 1b (CaR definition, truncated normal)
│ ├── figure1b_SI_logNormal.py # SI Fig. 1b (CaR definition, log-normal)
│ ├── portfolio_fig4ab.py # Figs. 4a and 4b (portfolio simulation, scatterplots, likelihoods, equalisation helper)
│ └── markowitz_portfolios_SI.py # SI Markowitz mean–variance optimisation
│
├── Outputs/ # empty folder; figures/CSVs appear here when code runs
│ └── README.md # explains what files are generated
│
└── README.md # this file

markdown
Copy code

---

## Figures generated

**Main text**
- `CaRDefinition1b.png` – Fig. 1b (CaR definition, truncated normal).
- `Fig_4a_portfolios_p5_CAR.png` – Fig. 4a (portfolio likelihood distributions with CaR arrows).
- `Fig_4b_Forest_Share.png` – Fig. 4b (scatterplot of portfolio risk-return, coloured by forest share).
- `Fig_3.png` – Fig. 3 (Carbon-at-Risk distribution at 200 years from SSC Monte Carlo).

**Supplementary Information**
- `CaRDefinition1bSI.png` – SI Fig. 1b (CaR definition, log-normal).
- `Fig_4b_SI_Cost_per_ton.png` – SI Fig. 4b (scatterplot of portfolio risk-return, coloured by cost per tonne).
- `Fig_S8.png` – SSC results at 1000 years (% CO₂ loss distribution).
- `Fig_S9.png` – SSC results at 10,000 years (% CO₂ loss distribution).
- `Fig_S11.png` – Markowitz scatter + efficient frontier (positive correlation).
- `Fig_S12.png` – Markowitz scatter + efficient frontier (zero correlation).
- `Fig_S13.png` – Markowitz scatter + efficient frontier (negative correlation).
- `Fig_S14.png` – Markowitz overlay of all efficient frontiers.
- `Fig_S15.png` – Markowitz composition stackplots along frontiers.
- `sample_portfolios.csv` – representative portfolios (low, mid, high risk per scenario).

---

## Code descriptions

### R code
- **`SSC_storage_security.R`**  
  Implements the Storage Security Calculator (SSC).  
  Functions:  
  - `SSCBase` – base-case leakage/trapping simulation.  
  - `SSCMC` – Monte Carlo uncertainty analysis.  
  - `Basic`, `BasicTable`, `FigLoss`, `MinMaxSA` – interrogation and plotting tools.  

  **Figures generated:**  
  - **Main text Fig. 3** – Carbon-at-Risk distribution at 200 years (Monte Carlo).  
  - **Fig. S8** – Monte Carlo % CO₂ loss distribution at 1000 years.  
  - **Fig. S9** – Monte Carlo % CO₂ loss distribution at 10,000 years.

### Python code
- **`figure1b_truncNormal.py`**  
  Defines CaR and “X% sure carbon removed” using a truncated normal distribution (main text Fig. 1b).
- **`figure1b_SI_logNormal.py`**  
  Alternative CaR illustration using a log-normal distribution (SI Fig. 1b).
- **`portfolio_fig4ab.py`**  
  Portfolio simulations of DACCS, Biochar, Forest projects. Produces:  
  - Fig. 4a (likelihood distributions with CaR).  
  - Fig. 4b (forest share scatterplot).  
  - SI Fig. 4b (cost per tonne scatterplot).  
  Also includes a helper function for “forest equalisation” (cost/survival adjustments) and console output of portfolio compositions and CaR.
- **`markowitz_portfolios_SI.py`**  
  Markowitz efficient frontier analysis with within- and between-type variance/correlation structure.  
  Produces scatterplots, frontier overlays, composition plots, and sample portfolio tables.

---

## Dependencies

### R
- R ≥ 3.6  
- Packages: `ggplot2`, `dplyr`

Install with:
```r
install.packages(c("ggplot2", "dplyr"))
Python
Python ≥ 3.8

numpy, pandas, matplotlib, scipy

seaborn (optional, for SI figures)

Install with:

bash
Copy code
pip install numpy pandas matplotlib scipy seaborn
Usage
Clone the repo and run scripts locally in RStudio or Python.

Example (Python):

bash
Copy code
cd code
python figure1b_truncNormal.py
python portfolio_fig4ab.py
python markowitz_portfolios_SI.py
Figures and CSVs will appear in ../Outputs/.

Example (R, in RStudio):

r
Copy code
source("code/SSC_storage_security.R")
BasicTable <- SSCBase()
plot(BasicTable[,"time"], BasicTable[,"injCO2"], type="l")
Notes
Outputs are not stored in GitHub to keep the repo light. Running the scripts will populate the Outputs/ folder.

All figures in the main text and SI can be reproduced with these scripts.

Large-scale simulations (e.g. Monte Carlo with 10,000 runs) may take several minutes.

References
Alcalde, J., et al. (2018). Assessing the capability of CO₂ storage security to deliver on climate mitigation. IJGGC, 91, 102834.

Lee et al. (2025) 'The Carbon at Risk measure can unlock financialmarkets for gigaton-scale carbon removal'. [This paper]

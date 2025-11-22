# IMD Vaccination Cost-Effectiveness Model

**Developed by the Modelling & Health Economics Unit, PHAC**

This repository contains an R-based decision-analytic model for evaluating the cost-effectiveness of Invasive Meningococcal Disease (IMD) vaccination strategies. It is designed as a base model for future economic evaluations by the Health Econ Team.

## 1. Project Structure

The project is organized into the following directories:

*   **`R/`**: Contains all the core R functions and modules.
    *   `01_setup.R`: Package loading and path configuration.
    *   `02_globals.R`: Global constants (cohort size, discount rates, states, strategies).
    *   `03_utils.R`: Utility functions (VE decay, clamping).
    *   `04_inputs.R`: Excel data loading and parameter processing.
    *   `05_model_core.R`: Markov model engine, state transitions, and cost/QALY calculations.
    *   `06_psa.R`: Probabilistic Sensitivity Analysis (PSA) logic.
    *   `06b_psa_visualization.R`: Functions for plotting PSA distributions.
    *   `07_owsa.R`: One-Way Sensitivity Analysis (OWSA) logic.
*   **`scripts/`**: Execution scripts.
    *   `run_analysis.R`: **Main entry point**. Runs Base Case, PSA, and OWSA for both Healthcare and Societal perspectives.
    *   `01_validate_model.R`: Comprehensive validation suite to check parameters, structure, and execution.
    *   `02_analyze_psa_results.R`: In-depth statistical analysis of PSA outputs (correlations, CVs).
    *   `03_model_diagnostics.R`: Deep dive into model mechanics (trace analysis, flow tracking).
    *   `04_quick_results_check.R`: Fast post-run sanity checks and one-page summary.
    *   `05_render_report.R`: Generates a polished Word report from analysis results.
*   **`reports/`**: RMarkdown templates and generated reports.
    *   `imd_analysis_report.Rmd`: Template for the cost-effectiveness report.
*   **`data/`**: Input data (Excel files).
    *   `raw/`: Place `IMD Data.xls` and `vaccine_costs.xlsx` here.
*   **`outputs/`**: Generated results.
    *   `tables/`: CSV files for costs, QALYs, ICERs.
    *   `figures/`: Plots for ICER frontiers, CEAC, EVPI, and Tornado diagrams.

## 2. Usage

To run the complete analysis:

1.  Ensure your input data is in `data/raw/`.
2.  Open `scripts/run_analysis.R`.
3.  Adjust configuration parameters if needed:
    *   `N_SIMULATIONS`: Number of PSA iterations (default: 2500).
    *   `COHORT_SIZE`: Size of the simulated cohort.
    *   `ANALYSIS_PERSPECTIVE`: "healthcare", "societal", or "both".
4.  Run the script:
    ```r
    source("scripts/run_analysis.R")
    ```

The script will:
*   Load all modules.
*   Run the deterministic Base Case.
*   Run the Probabilistic Sensitivity Analysis (PSA) using parallel processing.
*   Run the One-Way Sensitivity Analysis (OWSA).
*   Save all tables and figures to the `outputs/` directory.

### Validation & Diagnostics

To verify model integrity and explore results:

*   **Validation**: `source("scripts/01_validate_model.R")`
*   **Quick Check**: `source("scripts/04_quick_results_check.R")`
*   **PSA Analysis**: `source("scripts/02_analyze_psa_results.R")`
*   **Mechanics**: `source("scripts/03_model_diagnostics.R")`

### Reporting

To generate a formatted Word document report:

*   `source("scripts/05_render_report.R")`

## 3. Key Features

*   **Dual Perspective**: Automatically handles both Healthcare and Societal perspectives with correct cost logic.
*   **Parallel PSA**: Uses parallel processing for efficient simulation of large iterations.
*   **Reproducibility**: Uses fixed seeds for PSA sampling, ensuring common random numbers across perspectives.
*   **Robust Validation**: Includes dedicated scripts for validation, diagnostics, and result analysis.
*   **Automated Reporting**: RMarkdown script to generate submission-ready reports.

## 4. Data Availability

This model relies on two external Excel files for input parameters:

*   **`IMD Data.xls`**: Contains epidemiological data, utility weights, and other model parameters.
*   **`vaccine_costs.xlsx`**: Contains vaccine prices and administration costs.

**Note:** These files are **not included** in this repository as they contain confidential information (e.g., negotiated vaccine prices). To run the model, users must provide their own versions of these files with the expected structure and column names.

## 5. License & Citation

### License
This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.
This license allows for free use, modification, and distribution of the code, provided that the original copyright notice and authorship are preserved.

### Citation
If you use this code or model structure in your work, please cite the authors and the Public Health Agency of Canada:

> Ximenes, R., & Gebretekle, G. (2025). *IMD Vaccination Cost-Effectiveness Model*. Modelling & Health Economics Unit, Public Health Agency of Canada.

## 6. Authorship

**Raphael Ximenes**
*raphael.ximenes@phac-aspc.gc.ca*

**Gebremedhin Gebretekle**
*Gebremedhin.Gebretekle@phac-aspc.gc.ca*

**Modelling & Health Economics Unit**
Centre for Immunization Surveillance and Programs (CISP) | Centre de surveillance et de programmes d’immunisation (CISP)
Public Health Agency of Canada (PHAC) | Agence de la santé publique du Canada (ASPC)

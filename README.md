# PM10 and Precipitation Analysis: Piece-wise Polynomial & Bootstrap Approach

This repository contains a statistical analysis project developed for the **Nonparametric Statistics** course (Academic Year 2024-2025). The goal is to model the relationship between **PM10** pollution levels and **Precipitation** using flexible regression techniques and non-parametric validation.

---

## Project Overview

The analysis explores how precipitation affects the maximum weekly concentration of PM10. Instead of a simple linear regression, the implemented model accounts for potential structural changes in the relationship when precipitation exceeds a specific threshold (the median).

### Key Phases:
1. **Preprocessing & Aggregation**: Data cleaning and calculation of weekly maxima to filter daily variability and focus on pollution peaks.
2. **Piece-wise Modeling**: Fitting a second-degree polynomial model (`poly(degree = 2)`) with a **cutoff** point at the median of the precipitation values.
3. **Uncertainty & Bootstrap**: Estimation of coefficient distributions using **Non-parametric Bootstrap** (Case Resampling) with $10,000$ iterations.
4. **Parallel Computing**: Leveraging the `parallel` and `pbapply` packages to optimize computation time by distributing simulations across CPU cores.

---

## The Model

The regression model is defined as follows:

$$y = \beta_0 + \text{poly}(x, 2) + \beta_3 \mathbb{I}(x > c) + \beta_4 (x - c)\mathbb{I}(x > c) + \epsilon$$

Where:
* $x$: Weekly maximum precipitation (`max_Precip`).
* $c$: The **cutoff** value (median of $x$).
* $\mathbb{I}(x > c)$: An indicator function (dummy variable) for values above the threshold.
* $(x - c)\mathbb{I}(x > c)$: The interaction term allowing for a slope change.

---

## Requirements & Setup

### Libraries
To run the analysis, ensure you have the following R packages installed:
```R
install.packages(c("data.table", "dplyr", "progress", "pbapply", "parallel", "broom"))

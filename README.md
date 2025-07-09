# SimOmics

**SimOmics** is an R package that provides tools for simulating realistic multi-omics datasets. It enables researchers to generate synthetic omics blocks (e.g., transcriptomics, proteomics) with defined latent factors, noise, and block-wise correlation structures. This package is useful for benchmarking data integration methods such as `mixOmics`, teaching, and reproducibility testing.

---

## 🚀 Features

- Simulate multiple omics blocks with customizable dimensions
- Inject shared or independent latent structures
- Control inter-block correlation with block covariance structure
- Add Gaussian noise to simulate real-world signal-to-noise ratios
- Plot PCA, correlation heatmaps, and latent components
- Designed to interface with integration packages like `mixOmics`
- Export or benchmark using reproducible datasets

---

## 📦 Installation

You can install the development version from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install from GitHub (replace with your GitHub username/repo)
devtools::install_github("biosciences/SimOmics")
```

---

## 🧬 Example Usage

### Simulate and Visualize Data

```r
library(SimOmics)
library(mixOmics)

# Simulate data with 2 omics blocks
sim_data <- simulate_multiomics(
  n = 200,
  block_dims = list(transcriptome = 1000, proteome = 200),
  n_factors = 3,
  block_corr = 0.4,
  noise_sd = 0.5,
  seed = 123
)

# PCA of transcriptome block
plot_simulated_data(sim_data, type = "pca", block = "transcriptome")

# Correlation heatmap across all blocks
plot_simulated_data(sim_data, type = "correlation")
```

---

### Supervised Integration with mixOmics

```r
# Create synthetic class labels
Y <- factor(rep(c("A", "B"), each = 100))

# Perform supervised block integration
res <- block.plsda(X = sim_data$X_blocks, Y = Y, ncomp = 2)

# Visualize
plotIndiv(res, legend = TRUE)
```

---

## 📁 Directory Structure

```txt
SimOmics/
├── R/                     # Core R functions
├── man/                   # Function documentation
├── vignettes/             # Long-form tutorial (Rmd)
├── tests/testthat/        # Unit tests
├── paper.md               # JOSS paper
├── paper.bib              # Bibliography for JOSS paper
└── README.md              # This file
```

---

## 📚 Why Simulated Data Matters

In benchmarking and method development, real datasets often lack ground truth. Simulated data allows:
- Full control over latent effects and confounding
- Stress-testing of statistical integration methods
- Reproducible method comparison
- Fast prototyping and teaching

> 🔬 SimOmics was designed with biological realism in mind, generating omics datasets that reflect sparsity, structure, and integration complexity found in real studies.

---

## 🧪 Testing

Run all tests with:

```r
devtools::test()
```

---

## ✨ Citation

If you use SimOmics in your research, please cite the accompanying JOSS paper (in submission).

---

## 📄 License

MIT © 2025 Kaitao Lai

[![crp7](https://img.shields.io/badge/CRP-%237-%23ED9145?labelColor=%23ED9145&color=%2321609D)](https://cosmostatistics-initiative.org/residence-programs/crp7/)
[![arXiv](https://img.shields.io/badge/arXiv-astro--ph%2F2404.18165-%23ED9145?labelColor=%23ED9145&color=%2321609D)](https://arxiv.org/abs/2410.21962)
[![GitHub](https://img.shields.io/github/license/RafaelSdeSouza/capivara)](https://github.com/RafaelSdeSouza/capivara/blob/main/LICENSE) 
# [<img align="left" src="images/capivara.jpeg" width="45">](https://cosmostatistics-initiative.org/) Capivara: A Spectral-based Segmentation Method for IFU Data Cubes

Capivara implements a spectral-based segmentation method for Integral Field Unit (IFU) data cubes. Designed with astronomers in mind, it facilitates the decomposition of spectral data into regions of similar physical properties, leveraging advanced matrix operations via **torch** for GPU acceleration.

## Installation

Install Capivara from GitHub using the following commands:

```R
install.packages('remotes')
remotes::install_github("RafaelSdeSouza/capivara")
library(capivara)
```

## Usage

### Basic Usage

Here is a simple example of using **capivara** to process IFU data:

```R
library(capivara)

# Example data cube
set.seed(42)
n_row <- 50
n_col <- 50
n_wave <- 100
cube <- array(rnorm(n_row * n_col * n_wave), dim = c(n_row, n_col, n_wave))

# Segment the cube
result <- capivara_segment(cube)

# Plotting a decomposed region
library(ggplot2)
region <- result$regions[[1]]
region_df <- as.data.frame(as.table(region))

ggplot(region_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()
```

### Astronomical Example: MaNGA Data Cube

This example demonstrates how to use **capivara** to process an IFU datacube from the MaNGA survey, with dimensions [74, 74, 4563]. The data is flattened for segmentation:

```R
require(capivara)
# Read the MaNGA datacube
cube <- "manga-7443-12703-LOGCUBE.fits"

# Apply Capivara segmentation
result <- capivara_segment(data_2D)

# Extract and visualize a specific region
region_map <- matrix(result$regions[[1]], nrow = n_row, ncol = n_col)
region_df <- melt(region_map)

# Plot the segmented region
ggplot(region_df, aes(Var1, Var2, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "C") +
  theme_minimal()
```

## Dependencies

- **torch**: Efficient tensor computations with GPU support.
- **ggplot2**: Visualization.
- **FITSio**: Reading and handling FITS files.
- **reshape2**: Data manipulation.

## References

1. **MaNGA Survey**: Bundy, Kevin, et al. "Overview of the SDSS-IV MaNGA Survey: Mapping Nearby Galaxies at Apache Point Observatory." The Astrophysical Journal 798.1 (2015): 7. DOI: [10.1088/0004-637X/798/1/7](https://doi.org/10.1088/0004-637X/798/1/7)
2. **Capivara Methodology**: Souza, Rafael S. de, et al. "Capivara: A Spectral-Based Segmentation Method for IFU Data Cubes." arXiv preprint (2024). DOI: [10.48550/arXiv.2410.21962](https://arxiv.org/abs/2410.21962)
3. **Torch in R**: Paszke, Adam, et al. "PyTorch: An Imperative Style, High-Performance Deep Learning Library." Advances in Neural Information Processing Systems. 2019.

---
For more information, visit the [Cosmostatistics Initiative](https://cosmostatistics-initiative.org/) or check the [Capivara GitHub repository](https://github.com/RafaelSdeSouza/capivara).



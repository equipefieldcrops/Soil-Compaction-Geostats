# Soil-Compaction-Geostats
This document contains the complete geostatistical workflow in English, adapted for the paper:  "Assessing soil mechanical resistance at different depths and compaction mapping at field scale" (2024)

# Geostatistical Workflow for Soil Mechanical Resistance Mapping


This repository contains a complete geostatistical workflow adapted for the paper:


**"Assessing soil mechanical resistance at different depths and compaction mapping at field scale" (2024)**


The workflow performs:
- Exploratory data analysis
- Variogram estimation (Matheron, Pairwise-Robust, Cressie)
- Variogram model fitting
- Ordinary Kriging
- IDW interpolation
- Cross-validation (RMSE, ME)
- Raster export (GeoTIFF)


## Directory Structure

geostats_kriging_repository/ ├─ data/ ├─ scripts/ ├─ results/ ├─ workflow/ └─ README.md

## Running the Workflow

```r
source("scripts/kriging_workflow.R")



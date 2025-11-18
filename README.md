# Geostats Kriging Repository (Adapted for *Assessing Soil Mechanical Resistance at Different Depths and Compaction Mapping at Field Scale*)

This repository contains a **fully reproducible geostatistical workflow** in **English**, adapted for the paper:

**“Assessing soil mechanical resistance at different depths and compaction mapping at field scale” (Camargo et al. 2024)**

It includes:

* `scripts/kriging_workflow.R` — complete and commented R script
* `README.md` — full documentation in English
* `workflow/` — reproducibility setup (`renv`, Makefile)

---

## Recommended Directory Structure

```
geostats_kriging_repository/
├─ data/               # input datasets (shapefile, point data)
├─ scripts/            # R scripts
│  └─ kriging_workflow.R
├─ results/            # outputs (rasters, tables, figures)
├─ docs/               # reports, Rmd documents
├─ workflow/           # reproducibility files (Makefile, renv)
├─ README.md
└─ LICENSE
```

---

## `scripts/kriging_workflow.R`

```r
###############################################################################
# Geostatistical Analysis Workflow
# Adapted for the paper:
# "Assessing soil mechanical resistance at different depths and compaction
# mapping at field scale" (Camargo et al. 2024)
# DOI: dx.doi.org/10.1139/cjss-2024-0005
# Author: Mauricio Fornalski
# Description: complete workflow for variogram modelling, kriging interpolation
# and model validation.
###############################################################################

# ----------------------------- 1. Setup -------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)

packages <- c(
  "gstat", "sp", "sf", "raster", "ggplot2", "maptools",
  "geoR", "rgdal", "hydroGOF", "RandomFields", "readr"
)

installed <- rownames(installed.packages())
for (p in packages) if (!p %in% installed) install.packages(p)
lapply(packages, library, character.only = TRUE)

# Define project directories
root <- getwd()
data_dir    <- file.path(root, "data")
results_dir <- file.path(root, "results")
plots_dir   <- file.path(results_dir, "plots")
rasters_dir <- file.path(results_dir, "rasters")

# Create directories if missing
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rasters_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- 2. Load data --------------------------------
shp_path <- list.files(data_dir, pattern = "\.shp$", full.names = TRUE)
if (length(shp_path) == 0) stop("No shapefile found in data/")
shp <- rgdal::readOGR(shp_path[1])
plot(shp, main = "Study Area")

point_file <- list.files(data_dir, pattern = "(txt|csv)$", full.names = TRUE)
if (length(point_file) == 0) stop("No point data found in data/")
pt <- readr::read_table2(point_file[1])
stopifnot(all(c("X", "Y") %in% names(pt)))

coordinates(pt) <- ~X + Y
proj4string(pt) <- proj4string(shp)

# ----------------------------- 3. Create grid --------------------------------
cellsize <- 1
agr.grid <- spsample(shp, type = "regular", cellsize = cellsize, offset = c(cellsize, cellsize))
gridded(agr.grid) <- TRUE

grid_df <- as.data.frame(agr.grid)
coordinates(grid_df) <- ~ x1 + x2
gridded(grid_df) <- TRUE

# ----------------------------- 4. Exploratory --------------------------------
if (!"layer5" %in% names(pt)) stop("Variable 'layer5' missing from dataset.")
summary(pt$layer5)

png(file.path(plots_dir, "hist_layer5.png"), width = 800, height = 600)
hist(as.numeric(pt$layer5), main = "Histogram of Soil Mechanical Resistance", xlab = "layer5 (MPa)")
dev.off()

# ----------------------------- 5. Variogram ----------------------------------
vgm_matheron <- variogram(layer5 ~ 1, pt)
vgm_pr      <- variogram(layer5 ~ 1, pt, PR = TRUE)
vgm_cressie <- variogram(layer5 ~ 1, pt, cressie = TRUE)

init_model <- vgm(psill = var(as.numeric(pt$layer5)), model = "Sph", range = 400, nugget = 0)
fit_matheron <- fit.variogram(vgm_matheron, model = init_model)
fit_pr       <- fit.variogram(vgm_pr,       model = init_model)
fit_cressie  <- fit.variogram(vgm_cressie,  model = init_model)

png(file.path(plots_dir, "variograms.png"), width = 1200, height = 400)
par(mfrow = c(1,3))
plot(vgm_matheron, fit_matheron, main = "Matheron")
plot(vgm_pr,       fit_pr,       main = "Pairwise-Robust")
plot(vgm_cressie,  fit_cressie,  main = "Cressie")
dev.off()

chosen_fit <- fit_cressie
if (is.null(chosen_fit$psill)) chosen_fit <- fit_matheron

# ----------------------------- 6. Interpolation ------------------------------
idw_map <- idw(layer5 ~ 1, pt, grid_df)
kriged <- krige(layer5 ~ 1, pt, grid_df, model = chosen_fit)

write.table(as.data.frame(kriged), file = file.path(results_dir, "kriged_layer5.txt"), sep = "	", row.names = FALSE)
write.table(as.data.frame(idw_map), file = file.path(results_dir, "idw_layer5.txt"), sep = "	", row.names = FALSE)

# ----------------------------- 7. Validation ---------------------------
cv <- krige.cv(layer5 ~ 1, pt, model = chosen_fit, nfold = 5)
validation <- as.data.frame(cv)

png(file.path(plots_dir, "cv_residuals.png"), width = 800, height = 600)
bubble(cv, "residual", main = "Cross-validation Residuals")
dev.off()

rmse_value <- rmse(validation$observed, validation$var1.pred)
me_value   <- me(validation$observed, validation$var1.pred)
cat(sprintf("RMSE: %0.3f
ME: %0.3f
", rmse_value, me_value))

# ----------------------------- 8. Export rasters -----------------------------
gridded(kriged) <- TRUE
r_krig <- raster(kriged)
writeRaster(r_krig, filename = file.path(rasters_dir, "layer5_krig.tif"), format = "GTiff", overwrite = TRUE)

gridded(idw_map) <- TRUE
r_idw <- raster(idw_map)
writeRaster(r_idw, filename = file.path(rasters_dir, "layer5_idw.tif"), format = "GTiff", overwrite = TRUE)

save.image(file = file.path(results_dir, "session_workspace.RData"))
```

---

## `README.md`

```markdown
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
```

geostats_kriging_repository/
├─ data/
├─ scripts/
├─ results/
├─ workflow/
└─ README.md

````

## Running the Workflow

```r
source("scripts/kriging_workflow.R")
````

## Reproducibility

See the `workflow/` folder for the `Makefile` and instructions to use `renv`.

```

---

## `workflow/Makefile`

```

all: setup run

setup:
R -e "if(!require('renv')) install.packages('renv'); renv::init(bare=TRUE); renv::restore()"

run:
Rscript scripts/kriging_workflow.R

clean:
rm -rf results/*

```

---

If you'd like, I can also:
- insert the actual variable names from the PDF (soil resistance at different depths),
- generate an RMarkdown report,
- generate a ZIP file containing the full repository ready for GitHub.

Just let me know!

```




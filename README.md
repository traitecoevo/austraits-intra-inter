# Within- vs among-species trait variation along climate gradients in Australia

Analysis code using trait data from AusTraits and climate data from WorldClim.

## Workflow

1. **`Dataset_construction.qmd`** builds `data-processed/data_traits_climate.parquet` from AusTraits and WorldClim.
   The parquet file is committed, so this step is only needed to rebuild the dataset from scratch.
   It is skipped if the parquet file already exists, so delete it to force a rebuild.
2. **`Analysis.qmd`** fits the models and writes all figures and tables to `outputs/figures/` and `outputs/tables/`.
   It creates the `outputs/` folders itself.

## Data you need to download

`data/` is gitignored because of its size, so the files below must be downloaded into it.

| Data | Needed by | Where to put it |
|---|---|---|
| **WorldClim 2.1 bioclimatic variables, 30 arc-second** (Fick & Hijmans 2017). Download `wc2.1_30s_bio.zip` from <https://www.worldclim.org/data/worldclim21.html> ([direct link](https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip)) and unzip it. | `Dataset_construction.qmd` (BIO1, BIO12); `Analysis.qmd` (BIO1, BIO5, BIO10, BIO12) | `data/wc2.1_30s_bio/wc2.1_30s_bio_<n>.tif` |
| **AusTraits v7.0.0** (Falster et al. 2021, <https://doi.org/10.5281/zenodo.3568417>). Download the release's `.rds` file from Zenodo, or run `austraits::load_austraits(version = "7.0.0", path = "data")`. | `Dataset_construction.qmd` only | `data/austraits-7.0.0.rds` |

If you only run `Analysis.qmd` (using the committed parquet file), you need only the WorldClim files and are only used to run Figure S1, not for the models, so could be skipped unless Figure S1 is specifically to be rerun

## Data committed in `data-processed/`

- `data_traits_climate.parquet`: the analysis dataset (output of `Dataset_construction.qmd`).
- `australia.tif`: Australia mask used to crop the climate rasters.
- `datasets_with_leaf_area_edited.csv`: hand-annotated list of leaf-area datasets used to remove leaflet measurements of compound-leaved species.

## Model fits

The brms models in `Analysis.qmd` (`fit_model2d_bayesian_all_traits()` and `fit_model2d_bayesian_woody_all_traits()`) take a long time to fit.
They are saved to `outputs/model_fits/fit2d_bayesian/` and `outputs/model_fits/fit2d_bayesian_woody/` and read back by later chunks, so once fitted, those two chunks can be skipped on later runs.

## R packages

tidyverse, arrow, terra, tidyterra, sf, rnaturalearth, rnaturalearthdata, glmmTMB, brms, tmbstan, posterior, broom.mixed, easystats, modelbased, collapse, patchwork, ggsci, ggtext, see, scales, flextable, officer, austraits (only for `Dataset_construction.qmd`).

## References

- Falster, D. et al. (2021) AusTraits, a curated plant trait database for the Australian flora. *Scientific Data* 8, 254.
- Fick, S.E. & Hijmans, R.J. (2017) WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology* 37, 4302–4315.

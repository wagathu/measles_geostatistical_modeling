# Subnational inequities in measles vaccination coverage in Kenya

Code and model outputs for:

> *Subnational Inequities in Measles Vaccination Coverage in Kenya: A Bayesian Geostatistical Modelling Analysis.*

The study estimates coverage of the first (MCV1) and second (MCV2) doses of measles-containing vaccine in Kenya. It fits Bayesian geostatistical models (INLA-SPDE, via `inlabru`) to data from the 2022 Kenya Demographic and Health Survey (KDHS). Estimates are produced on a 5 × 5 km grid and aggregated to county and sub-county level using population weights.

## Repository structure

- **`code/`**: analysis code
  - `geostat_model.R`: the model script, covering model fitting, prediction on the 5 × 5 km grid, aggregation to county and sub-county level, and estimation of unvaccinated children
- **`data_outputs/`**: model estimates for MCV1 and MCV2 (each file below exists for both `mcv1_` and `mcv2_`)
  - `*_county_coverage.csv`: county-level coverage with 95% credible intervals
  - `*_subcounty_coverage.csv`: sub-county-level coverage with 95% credible intervals
  - `*_county_unvaccinated.csv`: estimated number of unvaccinated children per county, with 95% credible intervals
  - `*_subcounty_unvaccinated.csv`: estimated number of unvaccinated children per sub-county, with 95% credible intervals
  - `*_5km_mean.tif`: posterior mean coverage on the 5 × 5 km grid
  - `*_5km_uncertainty.tif`: width of the 95% credible interval on the 5 × 5 km grid
- **`README.md`**: this file
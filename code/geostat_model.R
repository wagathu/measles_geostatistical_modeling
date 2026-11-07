
# Importing packages ------------------------------------------------------

pacman::p_load(
  dplyr,
  raster,
  terra,
  data.table,
  sf,
  purrr,
  stringi,
  stringr,
  ggplot2,
  caret,
  recipes,
  ggridges,
  ggcorrplot,
  patchwork,
  INLA,
  fmesher,
  inlabru,
  ggspatial,
  lme4,
  RColorBrewer,
  pscl,
  parallel,
  doParallel,
  doSNOW
)
source('functions.R')

# Importing data ----------------------------------------------------------

mcv_2022 <- fread('data/county_mcv1.csv')
kenya_shapefile <- st_read('shapefiles/gadm41_KEN_0.shp', quiet = T)
county_shapefile2 <- st_read('../population_12_23/shapefiles/counties_without_lakes.shp',  quiet = T)
county_shapefile <- st_read(
  'shapefiles/gadm41_KEN_1.shp',  
  quiet = T, 
  query = "SELECT NAME_1 as county FROM gadm41_KEN_1"
)
subcounty_shapefile1 <- st_read('shapefiles/gadm41_KEN_2.shp', quiet = T)
cov <- fread("data/cov.csv")
all_df <- fread("data/all_df.csv")
pred1.1 <- fread("data/pred.csv")
covs <- brick("data/tif/tif_5KM.tif")
dhs_shp1 <- st_read(
  "dhs_shps/dhs_county.shp", 
  query = "select REGNAME as county from dhs_county", 
  quiet = T
)
tifs <- imap(c(
  list.files("data/dhs_tiff", full.names = T),
  list.files("data/tif_raw", full.names = T)
), ~ raster(.)) |>
  setNames(c(
    str_remove(list.files("data/dhs_tiff"), ".tif"),
    str_remove(list.files("data/tif_raw"), ".tif")
  ))
covs <- brick("data/tif/tif_5KM_KM.tif")

# Some initial cleaning ----------------------------------------------------

dhs_shp <- st_transform(dhs_shp1, crs = epsgKM(32737))

subcounty_shapefile <- data.table(subcounty_shapefile1) |> 
  _[, .(county = NAME_1, subcounty = NAME_2, geometry)] |>
  _[, subcounty := case_when(
    subcounty == "Webute West" ~ "Webuye West",
    TRUE ~ subcounty
  )] |>
  st_as_sf() |> 
  st_transform(st_crs(st_union(dhs_shp))) |> 
  st_intersection(st_union(dhs_shp))

county_shapefile_km <- county_shapefile |> 
  arrange(county) |> 
  st_transform(epsgKM(32737))

coords.clusters <- all_df[, .(lng, lat)] |>
  st_as_sf(coords = c('lng', 'lat'), crs = crs(kenya_shapefile)) |> 
  st_transform(crs = epsgKM(32737)) |> 
  st_coordinates() |> 
  data.table() |> 
  setNames(c('km_x', 'km_y'))

coords.shp <- all_df[, .(lng, lat)] |>
  st_as_sf(coords = c('lng', 'lat'), crs = crs(kenya_shapefile)) |> 
  st_transform(epsgKM(32737))
st_write(coords.shp, "shapefiles/coords.shp.shp", quiet = T)

kenya_shapefile2 <- kenya_shapefile |> 
  st_transform(epsgKM(32737))

coord.p <- pred1.1[, .(km_x, km_y)] |> as.matrix()
grid <- covs$prec

# The actual MCV1 County prevalence -----------------------------------------------------------

# Merging DHS estimates for mcv1 with shapefile
county_mcv_shp <- mcv_2022 |>
  _[, data.table::merge.data.table(
    .SD,
    dhs_shp,
    by = 'county',
  )] |>
  _[, let(prop = round(prop))] |> 
  _[, `:=`(
    mcv1_12_23 = case_when(
      prop >= 0 & prop <= 49 ~ "0-49%",
      prop > 49 & prop <= 59 ~ "50-59%",
      prop > 59 & prop <= 69 ~ "60-69%",
      prop > 69 & prop <= 79 ~ "70-79%",
      prop > 79 & prop <= 89 ~ "80-89%",
      prop > 89 & prop <= 94 ~ "90-94%",
      TRUE ~ "95% and above"
    ) |>
      factor(levels = c(
        "0-49%",
        "50-59%",
        "60-69%",
        "70-79%",
        "80-89%",
        "90-94%",
        "95% and above"
      ))
    )
    ] |> st_as_sf()

sf::st_write(
  county_mcv_shp |>
    dplyr::select(county, prop, mcv1_12_23) |> 
    st_as_sf(),
  "shapefiles/county_mcv_shp.shp", overwrite = T)

# Small arrangement of datas ----------------------------------------------

dm <-  all_df |> 
  mutate(lng2 = lng, lat2 = lat) |> 
  st_as_sf(coords = c("lng", "lat"), crs = crs(kenya_shapefile)) |>
  st_transform(crs = epsgKM(32737)) %>%
  mutate(
    km_x = st_coordinates(.)[, 1],
    km_y = st_coordinates(.)[, 2]
  ) |> 
  as.data.table() |> 
  _[, map_if(.SD, is.numeric, as.double)] |> 
  _[, !"geometry"] |> 
  mutate(
    IsUrban = fifelse(IsUrban == 1, "Yes", "No"),
    IsUrban = relevel(factor(IsUrban), ref = "No")
  ) |> 
  rename(lng = lng2, lat = lat2)
pop_1 <- pred1.1$pop_1

pred1 <- pred1.1 |>
  _[, map_if(.SD, is.numeric, as.double)] |> 
  mutate(
    IsUrban = fifelse(IsUrban == 1, "Yes", "No"),
    IsUrban = relevel(factor(IsUrban), ref = "No")
  ) |> 
  dplyr::select(-pop_1) 

dm2_recipe <- dm |>
  dplyr::select(!ends_with("_cat")) |> 
  mutate(IsUrban = relevel(factor(IsUrban), ref = "No")) |> 
  recipe(cluster_number + n +N + lng + lat + km_x + km_y + prop ~ .) |> 
  step_log(travel_hfa_malaria_atlas, travel_time_cities, ken_pop_density, facility_density, offset = 1e-3) |> 
  step_center(all_double_predictors()) |> 
  step_scale(all_double_predictors()) |> 
  step_dummy(all_nominal_predictors()) |> 
  prep()

dm2 <- juice(dm2_recipe) |> 
  rename(IsUrban = IsUrban_Yes)
fwrite(dm2, "data/dm2.csv", row.names = F)

pred <- pred1 |> 
  mutate(IsUrban = relevel(factor(IsUrban), ref = "No")) %>%
  bake(dm2_recipe, .) |> 
  rename(IsUrban = IsUrban_Yes)

any(apply(dm2, 2, \(z) any(is.infinite(z))))
any(apply(pred, 2, \(z) any(is.infinite(z))))

# Step-wise regression -----------------------------------------------------

glm_model <- glm(
  cbind(
    n,N-n
  ) ~ 
    dpt1 +
    literate_women +
    lowest_wealth_quantile_prop +
    pn3 +
    vacc_card +
    lst_day +
    prec +
    travel_hfa_malaria_atlas +
    travel_time_cities + 
    ken_pop_density +
    facility_density + 
    IsUrban,
  family = "binomial",
  data = dm2
)
summary(glm_model)

mode_step <- MASS::stepAIC(glm_model, direction = 'backward', scope = list(lower = ~IsUrban))
summary(mode_step)

# Glm model for the residuals and obtaining the variogram -----------------


# Fit GLMM with observation-level random effect
glmm_model <- glmer(
  cbind(n, N - n) ~ dpt1 + literate_women + lowest_wealth_quantile_prop +
    travel_hfa_malaria_atlas + ken_pop_density + IsUrban +
    (1 | cluster_number),
  family = binomial(link = "logit"),
  data = dm2
)

# Extract random effects
re <- ranef(glmm_model)$cluster_number[, 1]

# Compute variogram from random effects
variogram <- ggvario(
  coords = coords.clusters,
  data = re,
  bins = 25,
  maxdist = 250
)

glm_model <- glm(
  cbind(n, N - n) ~ dpt1 + literate_women + lowest_wealth_quantile_prop +
     travel_hfa_malaria_atlas + ken_pop_density + IsUrban,
  family = binomial(link = "logit"),
  data = dm2
)

resids <- resid(glm_model)

variogram <- ggvario(
  coords = coords.clusters,
  data =  resids,
  bins = 25,
  maxdist = 250,
  xlab = "Distance between point pairs (Km)",
  show_nbins = F
)  +
  theme(plot.title = element_text(color = "black", hjust = .5)) +
  ggtitle("Empirical variogram for MCV1") +
  labs(y = "Semivariance")

ggsave(
  plot = variogram,
  'plots/mcv1_variogram.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 2 * 5,
  units = "in"
)

# Creating the mesh and the SPDE model -------------------------------------------------------

# Bounding box of Kenya
bb <- st_bbox(kenya_shapefile2)
diameter <- bb["ymax"] - bb["ymin"] 

# Fuglstad's recommendation
r0_fuglstad <- diameter / 5

# Mesh
mesh <- fm_mesh_2d_inla(
  boundary = kenya_shapefile2,
  loc = coords.clusters,
  max.edge = c(20, 40),
  cutoff = 10
)

# Matern SPDE model object using inla.pcmatern
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 1.5,
  prior.range = c(r0_fuglstad, .5),   # P(r0_fuglstad < r0) = 0.5
  prior.sigma = c(1, .01)   # P(sigma > 1) = 0.01
)

# Creating the model components ---------------------------------------------------------------

cmp <- ~ 
  dpt1(dpt1, model = "linear") +
  literate_women(literate_women, model = "linear") +
  lowest_wealth_quantile_prop(lowest_wealth_quantile_prop, model = "linear") +
  travel_hfa_malaria_atlas(travel_hfa_malaria_atlas, model = "linear") +
  ken_pop_density(ken_pop_density, model = "linear") +
  IsUrban(IsUrban, model = "linear") +
  field(cbind(km_x, km_y), model = spde) +
  nugget(cluster_number, model = "iid", hyper = list(theta = list(prior = "pc.prec", param = c(1, .01)))) +
  Intercept(1)

# The INLA models ----------------------------------------------------------

# The formula
frm <-
  n ~
  dpt1 +
  literate_women +
  lowest_wealth_quantile_prop +
  travel_hfa_malaria_atlas +
  ken_pop_density +
  IsUrban +
  field +
  nugget +
  Intercept 

mod <- bru(
  components = cmp,
  formula = frm,
  family = "binomial",
  Ntrials = dm2$N,
  data = dm2,
  options = list(
    control.fixed = list(mean = 0, prec = 1, mean.intercept = 0, prec.intercept = 0.001),
    control.compute = list(config = TRUE, return.marginals.predictor = TRUE, cpo = T, waic = T, dic = T)
  )
)

# Some techqs -------------------------------------------------------------

#Parameters
coeff.reg <- tidy.inla(mod, exp = T) |>
  mutate(
    terms = str_replace_all(terms, "x1", "DPT1")
  ) |> 
  filter(terms != "b0", terms != "(Intercept)", terms != "Intercept")
fwrite(coeff.reg, "data/coeff.reg.csv", row.names = F)

forest_plot <- coeff.reg |> 
  mutate(
    terms = case_when(
      terms == 'vacc_card' ~ 'Vaccination card',
      terms == "ses_deprivation" ~ "Deprivation",
      terms == 'travel_time_cities' ~ 'Travel time to cities',
      terms == "lowest_wealth_quantile_prop" ~ "Proportion of population in the lowest wealth quintile",
      terms == 'travel_hfa_malaria_atlas' ~ 'Travel time to health facilties',
      terms == 'prec' ~ 'Precipitation',
      terms == 'pn3' ~ 'Pneumococcal 3',
      terms == 'night_light' ~ 'Night light',
      terms == 'lst_day' ~ 'Day land surface temperature',
      terms == 'literate_women' ~ 'Literate women',
      terms == 'ken_pop_density' ~ 'Population density',
      terms == 'dpt1' ~ 'DPT1',
      terms == "IsUrban" ~ "IsUrban:Yes",
      TRUE ~ terms
    ) |> 
      str_wrap(25)
  ) |> 
  ggplot(aes(x = mean, y = terms)) +
  geom_point(color = "#0072B2", size = 3) +
  geom_errorbar(
    aes(xmin = `0.025quant`, xmax = `0.975quant`),
    width = .15,
    linewidth = 1,
    color = "#0072B2"
  ) +
  geom_vline(xintercept = 1, color = "#009E73", linetype = 2, linewidth = 1) +
  theme_minimal() +
  theme(
    axis.title = element_text(color = "black", size = 18),
    axis.text = element_text(color = "black", size = 18)
  ) +
  labs(x = "Odds", y = "Variable")

ggsave(
  plot = forest_plot,
  filename  = "plots/forest_plot.png",
  dpi = 5e2,
  width = 2*5,
  height = 2*5,
  bg = NULL
)

# Sampling from the posterior distribution ------------------------------------------------

nsamp <- 1000

# Posterior sampling
linpred <- generate(
  mod,
  pred,
  formula = ~ {
    k <- 16 * sqrt(3) / (15 * pi)
    sigma2 <- 1 / Precision_for_nugget
    
    eta <- dpt1 + literate_women + lowest_wealth_quantile_prop +  IsUrban +
      travel_hfa_malaria_atlas + ken_pop_density + field + Intercept
    
    eta / sqrt(1 + k^2 * sigma2)
  },
  n.samples = nsamp,
  num.threads = "7:1"
)

std.error <- apply(linpred, 1, \(.) sd(.))
inv.linpred <- gtools::inv.logit(linpred)
pred.grid <- data.frame(t(apply(inv.linpred, 1, FUN = function(x){ c(mean(x), sd(x), quantile(x, probs = c(0.025,0.5,0.975)))})))
colnames(pred.grid) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
fitted.pred.mean <- pred.grid[["mean"]]
fitted.pred.sd <- pred.grid[["sd"]]
fitted.pred.median <- pred.grid[["0.5quant"]]
fitted.pred.low <- pred.grid[["0.025quant"]]
fitted.pred.up <- pred.grid[["0.975quant"]]
fitted.uncertain  <- fitted.pred.up - fitted.pred.low

#' Rasterizing the results
mcv1.mean <- rasterize(x = coord.p, y = grid, field = fitted.pred.mean)
plot(mcv1.mean)
writeRaster(mcv1.mean, 'rasters/mcv1_coverage.tif', overwrite = T)

#' Rasterizing the low CI
mcv1.lower <- rasterize(x = coord.p, y = grid, field = fitted.pred.low)

#' Rasterizing the upper ci
mcv1.upper <- rasterize(x = coord.p, y = grid, field = fitted.pred.up)

#' Rasterizing the uncertainty
mcv1.uncertain <- rasterize(x = coord.p, y = grid, field = fitted.uncertain)

#' Rasterizing the std error
mcv1.std_error <- rasterize(x = coord.p, y = grid, field = std.error)

#' The plot
rasters <- list(mcv1.lower, mcv1.mean, mcv1.upper) |>
  setNames(c(
    'Lower CI',
    'mcv1 probability',
    'Upper CI'
  ))
saveRDS(rasters, "../raster results/mcv1_results_2022.RDS")
par(mfrow = c(1,3))
map(rasters, plot)
par(mfrow = c(1,1))

# Obtaining predictions for the clusters -------------------------------

yy <- generate(
  mod,
  dm2,
  formula = ~ {
    eta <- dpt1 + literate_women + lowest_wealth_quantile_prop +
      travel_hfa_malaria_atlas + ken_pop_density + IsUrban + field +
      nugget + Intercept
    prob <- 1/(1 + exp(-eta))
    rbinom(length(prob), prob = prob, size = dm2$N)
  },
  n.samples = nsamp
)

prev_yy <- apply(yy, 2, \(z) z/dm2$N)
pred.grid_cluster <- data.frame(t(apply(prev_yy, 1, FUN = function(x){ c(mean(x), sd(x), quantile(x, probs = c(0.025,0.5,0.975)))})))
colnames(pred.grid_cluster) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant") 
rownames(pred.grid_cluster) <- NULL

pred_actual <- data.frame(pred.grid_cluster) |> 
  data.table() |> 
  mutate(
    actual = dm2$n/dm2$N,
    actual_cat = case_when(
      actual == 0 & actual < 0.1 ~ "0-0.1",
      actual >= 0.1 & actual < 0.2 ~ "0.1-0.2",
      actual >= 0.2 & actual < 0.3 ~ "0.2-0.3",
      actual >= 0.3 & actual < 0.4 ~ "0.3-0.4",
      actual >= 0.4 & actual < 0.5 ~ "0.4-0.5",
      actual >= 0.5 & actual < 0.6 ~ "0.5-0.6",
      actual >= 0.6 & actual < 0.7 ~ "0.5-0.6",
      actual >= 0.7 & actual < 0.8 ~ "0.7-0.8",
      actual >= 0.8 & actual < 0.9 ~ "0.8-0.9",
      actual >= 0.9 & actual <= 1 ~ "0.9-1"
    ),
    pred_cat = case_when(
      mean == 0 & mean < 0.1 ~ "0-0.1",
      mean >= 0.1 & mean < 0.2 ~ "0.1-0.2",
      mean >= 0.2 & mean < 0.3 ~ "0.2-0.3",
      mean >= 0.3 & mean < 0.4 ~ "0.3-0.4",
      mean >= 0.4 & mean < 0.5 ~ "0.4-0.5",
      mean >= 0.5 & mean < 0.6 ~ "0.5-0.6",
      mean >= 0.6 & mean < 0.7 ~ "0.5-0.6",
      mean >= 0.7 & mean < 0.8 ~ "0.7-0.8",
      mean >= 0.8 & mean < 0.9 ~ "0.8-0.9",
      mean >= 0.9 & mean <= 1 ~ "0.9-1"
    )
  ) |> 
  dplyr::summarise(
    actual_mean = mean(actual, na.rm = TRUE),
    predicted_mean = mean(mean, na.rm = TRUE),
    .by = actual_cat
  ) 

pred_actual_plot <- ggplot(pred_actual, aes(x = actual_mean, y = predicted_mean)) +
  geom_point(color = "blue", size = 3) +  # Points for each category
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "black"
  ) +  # Perfect calibration
  labs(x = "Mean actual)", y = "Mean predicted") +
  theme_minimal() +
  theme(
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 14),
  )

ggsave(
  plot = pred_actual_plot,
  filename = "plots/pred_actual_plot.png",
  dpi = 5e2,
  width = 2 * 5,
  height = 2 * 5,
  bg = NULL
)

prev_yy_plot <- prev_yy |>
  data.table() %>%
  melt(measure.vars = names(.)) |>
  ggplot(aes(x = value)) +  
  geom_density(aes(group = variable), colour = "grey60") +  
  geom_density(data = dm2, aes(x = n/N), colour = "black", linewidth = 1) +
  theme_minimal() +
  theme(
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 14),
  )

ggsave(
  plot = prev_yy_plot,
  filename = "plots/prev_yy_plot.png",
  dpi = 5e2,
  width = 2 * 5,
  height = 2 * 5,
  bg = NULL
)

# Obtaining variograms after modelling ------------------------------------

df_after <- yy |> 
  as.data.frame() |> 
  setNames(c(paste0("sample", 1:1000))) |> 
  data.table() %>%
  melt(measure.vars = names(.), value.name = "n", variable.name = 'sample') %>%
  split(.$sample)

# Register parallel backend
numCores <- detectCores() - 1
cl <- makeCluster(numCores, type = "SOCK")  # Use SOCK for doSNOW
registerDoSNOW(cl)

# Initialize progress bar
pb <- txtProgressBar(max = length(df_after), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run foreach loop
dfvario_list <- foreach(
  x = df_after, 
  y = names(df_after), 
  .packages = c("data.table", "geoR"),
  .options.snow = opts
) %dopar% {
  library(data.table)  # Ensure packages are loaded inside foreach
  library(geoR)
  library(INLA)
  
  x <- data.table(x)[, !'sample', with = FALSE]
  
  df <- cbind(x, data.table(dm2)[, !c("lng", "lat", "cluster_number", "prop", 'n'), with = FALSE])
  
  glm_model <- glm(
    cbind(n, N - n) ~ dpt1 + literate_women + lowest_wealth_quantile_prop +
      travel_hfa_malaria_atlas + ken_pop_density + IsUrban,
    family = binomial(link = "logit"),
    data = df
  )
  
  resids <- resid(glm_model)
  coords <- as.matrix(coords.clusters)
  min_dist <- min(dist(coords))
  bins <- 15
  maxdist <- max(dist(coords)) / 3
  uvec <- seq(min_dist, maxdist, length.out = bins)
  
  empvario <- variog(coords = coords, data = resids, uvec = uvec, messages = FALSE)
  
  empvario$v
}

# Close progress bar and cluster
close(pb)
stopCluster(cl)

dfEmp <- imap(dfvario_list, \(x,y){
  t(matrix(x)) |> 
    data.table()
}) |> 
  rbindlist()

# The actual Variogram
min_dist <- min(dist(coords.clusters))
bins <- 15
maxdist <- max(dist(coords.clusters)) / 3
uvec <- seq(min_dist, maxdist, length.out = bins)
empvario <- geoR::variog(coords = coords.clusters, data = resids, uvec = uvec, messages = FALSE)


dfEmp_summary <- data.table(
  variable = colnames(dfEmp),
  min = apply(dfEmp, 2, min),
  max = apply(dfEmp, 2, max),
  quant2.5 = apply(dfEmp, 2, quantile, probs = 0.025),
  quant97.5 = apply(dfEmp, 2, quantile, probs = 0.975),
  distance = empvario$u,
  empirical = empvario$v
)

empDiagnostic <- dfEmp_summary |> 
  ggplot(aes(x = distance, y = empirical)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#27AAE1", alpha = 0.3) +
  geom_line(aes(y = quant2.5), linetype = 2) +
  geom_line(aes(y = quant97.5), linetype = 2) +
  geom_point() +
  theme_minimal() +
  theme(
    axis.title = element_text(color = "black", size = 17),
    axis.text = element_text(color = "black", size = 17),
  ) +
  labs(y = "Semivariance", x = "Distance (KM)")

ggsave(
  plot = empDiagnostic,
  filename = "plots/empDiagnostic.png",
  dpi = 5e2,
  width = 2 * 5,
  height = 2 * 5,
  bg = NULL
)

# Obtaining the key statistic's -----------------------------------------------

#range for spatial RE
range_std <- mod$summary.hyperpar[1:2,][,1:5]
range_std[2, c(1,3,4,5)] <- (range_std[2, c(1,3,4,5)]) ^ 2

# variance for IID RE
var.ind <- inla.tmarginal(function(x)
  1 / x, mod$marginals.hyperpar[[3]])
var.iid <- inla.zmarginal(var.ind, silent = TRUE)
variance.iid <- c(var.iid$mean,
                  var.iid$sd,
                  var.iid$quant0.025,
                  var.iid$quant0.5,
                  var.iid$quant0.975)
param.all <- rbind( range_std, variance.iid)
#param.all <- rbind(range_std)
param.all <- param.all[, c(1,3,5)]
param.all <- apply(param.all, 2, \(x) round(x,2))
#rownames(param.all) <- c("Spatial range", "Spatial variance")
rownames(param.all) <- c("Spatial range", "Spatial variance", "Nugget")
write.csv(param.all, "data/parameter_output.csv")

# Cropping to DHS shapefile areas -------------------------------------------------------

temp_list <- list(mcv1.mean,
                  mcv1.lower,
                  mcv1.upper,
                  mcv1.uncertain,
                  mcv1.std_error
) |> 
  imap(\(x, y){
    x |> 
      crop(dhs_shp) |> 
      mask(dhs_shp)
  }) |> 
  setNames(c(
    'mcv1.mean',
    'mcv1.lower',
    'mcv1.upper',
    'mcv1.uncertain',
    "mcv1.std_error"
  ))

# Plotting the rasters --------------------------------------------

rast_ll <- ggplot() +
  geom_tile(data = as.data.frame(rast(temp_list$mcv1.lower) * 1e2, xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = dhs_shp,
    fill = NA
  ) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 100)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "") 


rast_mean <- ggplot() +
  geom_tile(data = as.data.frame(rast(temp_list$mcv1.mean) * 1e2, xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = dhs_shp,
    fill = NA
  ) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 100)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "")

rast_upper <- ggplot() +
  geom_tile(data = as.data.frame(rast(temp_list$mcv1.upper) * 1e2, xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = dhs_shp,
    fill = NA
  ) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 100)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "")

rr_all <- rast_ll + rast_mean + rast_upper
ggsave(
  plot = rr_all,
  'plots/rr_all.png',
  dpi = 5e2,
  height = 5,
  width = 3*5,
  bg = NULL,
  units = "in"
)

ggsave(
  plot = rast_mean,
  'plots/rast_mean.png',
  dpi = 5e2,
  height = 5,
  width = 5,
  bg = NULL,
  units = "in"
)


# Uncertainty
rast_uncertain <- ggplot() +
  geom_tile(data = as.data.frame(rast(temp_list$mcv1.uncertain) * 1e2, xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = dhs_shp,
    fill = NA
  ) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(0, 100)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "")

ggsave(
  plot = rast_uncertain,
  'plots/rast_uncertain.png',
  dpi = 5e2,
  height = 5,
  width = 5,
  bg = NULL,
  units = "in"
)

# Uncertainity plus points
rast_uncertain_points <- rast_uncertain +
  geom_point(data = coords.clusters,  aes(x = km_x, y = km_y))

ggsave(
  plot = rast_uncertain_points,
  'plots/rast_uncertain_points.png',
  dpi = 5e2,
  height = 5,
  width = 5,
  bg = NULL,
  units = "in"
)

# std_error
rast_std_error <- ggplot() +
  geom_tile(data = as.data.frame(rast(temp_list$mcv1.std_error), xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = dhs_shp,
    fill = NA
  ) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "")

ggsave(
  plot = rast_std_error,
  'plots/rast_std_error.png',
  dpi = 5e2,
  height = 5,
  width = 5,
  bg = NULL,
  units = "in"
)

# Std error and the points
rast_std_error_points <- rast_std_error +
  geom_point(data = coords.clusters,  aes(x = km_x, y = km_y))

ggsave(
  plot = rast_std_error_points,
  'plots/rast_std_error_points.png',
  dpi = 5e2,
  height = 5,
  width = 5,
  bg = NULL,
  units = "in"
)

# Population weighted for counties data -------------------------

# Population raster file
# Getting where each pixel falls within a county
cc <- as(as(dhs_shp, "Spatial"), "SpatialPolygons")
pred_locations <- st_as_sf(data.frame(coord.p), coords = c('km_x', 'km_y'), crs = crs(dhs_shp)) |>
  as('Spatial') |> 
  as('SpatialPoints') 
cc_id <- over(pred_locations, cc)

joint_data <- data.frame(
  ID = cc_id,
  pop = pop_1,
  inv.linpred
)

cc_weighted_all <- data.table(joint_weighted_prev(joint_data = joint_data)) |> 
  _[, let(ID = as.numeric(ID), sample_id = as.numeric(sample_id))] |> 
  _[, setorder(.SD, ID)] |>  
  dcast(...~sample_id, value.var = "prev")

apply(cc_weighted_all[, -1], 1, \(x) length(which(x >= 0.95)) / 1000)

cc_weighted <- t(apply(cc_weighted_all[, -1], 1, \(x) c(
  mean(x), 
  sd(x), 
  quantile(x, probs = c(0.025, 0.5, 0.975)),
  length(which(x >= 0.95)) / 1000,
  length(which(x >= 0.80)) / 1000
))) |> cbind(cc_weighted_all[, 1])

colnames(cc_weighted) <- c(
  "mean",
  "sd",
  "0.025quant",
  "0.5quant",
  "0.975quant",
  "thresh95",
  "thresh80",
  "ID"
)
cmat <- cc_weighted

# Population weighted for Sub-counties data -------------------------

subcounty_shapefile2 <- subcounty_shapefile %>%
  mutate(ID = 1:nrow(.)) |> 
  st_as_sf() |> 
  st_transform(crs = epsgKM(32737))

# Getting where each pixel falls within a sub county
sc <- as(as(subcounty_shapefile2, "Spatial"), "SpatialPolygons")
sub_pred_locations <- st_as_sf(data.frame(coord.p), coords = c('km_x', 'km_y'), crs = crs(subcounty_shapefile2)) |>
  as('Spatial') |> 
  as('SpatialPoints') 
sc_id <- over(sub_pred_locations, sc)
setdiff(unique(subcounty_shapefile2$ID), sc_id) |> length()

# Get centroids of missing subcounties and convert coord.p to sf
mc <- st_centroid(subcounty_shapefile2) %>% 
  filter(ID %in% setdiff(unique(subcounty_shapefile2$ID), sc_id))
coord_sf <- st_as_sf(data.frame(coord.p), coords = c("km_x", "km_y"), crs = st_crs(subcounty_shapefile2))

# Map nearest pixel indices to missing subcounty IDs
nearest_pixel_map <- data.frame(
  pixel_index = st_nearest_feature(mc, coord_sf),
  ID = mc$ID
)

# Take values from inv.linpred and pop_1 at those nearest pixels
missing <- data.frame(
  ID = nearest_pixel_map$ID,
  pop = 1,  # Set population = 1 to act as weight-neutral
  inv.linpred[nearest_pixel_map$pixel_index, , drop = FALSE]
)

sub_joint_data <- data.frame(
  ID = sc_id,
  pop = pop_1,
  inv.linpred
) |> 
  rbind(missing) |> 
  arrange(ID)
setdiff(unique(subcounty_shapefile2$ID), sub_joint_data$ID) |> length()

sc_weighted_all <- data.table(joint_weighted_prev(joint_data = sub_joint_data)) |> 
  _[, let(ID = as.numeric(ID), sample_id = as.numeric(sample_id))] |> 
  _[, setorder(.SD, ID)] |>  
  dcast(...~sample_id, value.var = "prev")

sc_weighted <- t(apply(sc_weighted_all[, -1], 1, \(x) c(
  mean(x), 
  sd(x), 
  quantile(x, probs = c(0.025, 0.5, 0.975)),
  length(which(x >= 0.95)) / 1000,
  length(which(x >= 0.80)) / 1000
))) |> cbind(sc_weighted_all[, 1])

colnames(sc_weighted) <- c(
  "mean",
  "sd",
  "0.025quant",
  "0.5quant",
  "0.975quant",
  "thresh95",
  "thresh80",
  "ID"
)
smat <- sc_weighted

# Joining weighted mcv data with shape file -----------------------------------------------------

# County
countyLevel <- setDT(cbind(dhs_shp, cmat)) |>
  _[, .(
    county,
    mcv1 = round(as.numeric(mean * 1e2), 1),
    lower_ci = X0.025quant * 1e2,
    upper_ci = X0.975quant * 1e2
  )] |>
  _[, mcv1_12_23 := case_when(
    mcv1 >= 0  & mcv1 < 50 ~ "0-49%",
    mcv1 >= 50 & mcv1 < 60 ~ "50-59%",
    mcv1 >= 60 & mcv1 < 70 ~ "60-69%",
    mcv1 >= 70 & mcv1 < 80 ~ "70-79%",
    mcv1 >= 80 & mcv1 < 90 ~ "80-89%",
    mcv1 >= 90 & mcv1 < 95 ~ "90-94%",
    TRUE ~ "95% and above"
  ) |>
    factor(levels = c(
      "0-49%",
      "50-59%",
      "60-69%",
      "70-79%",
      "80-89%",
      "90-94%",
      "95% and above"
    ))] |>
  data.table::merge.data.table(dhs_shp, by = "county") |>
  st_as_sf()

round((cor(county_mcv_shp$prop, countyLevel$mcv1)) ^ 2, 2) # The R-Squared
#st_write(countyLevel, 'shapefile_insights/countyLevel.shp', append = F)

# Subcounty
sub_shp <- subcounty_shapefile1 |>
  data.table() |>
  st_as_sf()

subcountyLevel <- setDT(merge(subcounty_shapefile2, smat, by = 'ID')) |>
  _[, .(county, subcounty, mcv1 = round(as.numeric(mean * 1e2)), geometry)] |>
  _[, mcv1_12_23 := case_when(
    mcv1 >= 0  & mcv1 < 50 ~ "0-49%",
    mcv1 >= 50 & mcv1 < 60 ~ "50-59%",
    mcv1 >= 60 & mcv1 < 70 ~ "60-69%",
    mcv1 >= 70 & mcv1 < 80 ~ "70-79%",
    mcv1 >= 80 & mcv1 < 90 ~ "80-89%",
    mcv1 >= 90 & mcv1 < 95 ~ "90-94%",
    TRUE ~ "95% and above"
  ) |>
    factor(levels = c(
      "0-49%",
      "50-59%",
      "60-69%",
      "70-79%",
      "80-89%",
      "90-94%",
      "95% and above"
    ))] |>
  st_as_sf()
#st_write(subcountyLevel, 'shapefile_insights/subcountyLevel.shp', append = F)

# Joining the threshold probability with the shapefile --------------------

# County
countyLevelThresh <- setDT(merge(
  data.table(dhs_shp)[, ID := 1:47],
  cc_weighted[, .(ID, thresh95, thresh80)],
  by = 'ID'
)) |>
  _[, .(county, thresh95, thresh80, geometry)] |>
  st_as_sf()
# st_write(countyLevelThresh,
#          'shapefile_insights/countyLevelThresh.shp',
#          quiet = T,
#          append = F)

# Subcounty
subcountyLevelThresh <- setDT(merge(
  subcounty_shapefile2,
  sc_weighted[, .(ID, thresh95_sub = thresh95, thresh80_sub = thresh80)],
  by = 'ID'
)) |>
  _[, .(county, subcounty, thresh95_sub, thresh80_sub, geometry)] |>
  st_as_sf()

# Plotting for mcv1 ----------------------------------------------------------------

#' Plotting for the MCV1 probability only
plt_mcv1_pred <- countyLevel |>
  ggplot() +
  geom_sf(aes(fill = mcv1_12_23), show.legend = T) +
  scale_fill_manual(
    values = brewer.pal(n = 7, "RdYlBu"),
    drop = F
  ) +
  theme_void() +
  theme(
    plot.title = element_text(color = "black", hjust = .5, size = 17),
    legend.key.size = unit(1, "lines"),
    legend.title = element_text(colour = "black", size = 13),
    legend.text = element_text(color = 'black', size = 13)
  ) +
  labs(fill = "") 
ggsave(
  plot = plt_mcv1_pred,
  'plots/plt_mcv1_pred.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 5,
  units = "in"
)

plt_mcv1_dhs <- county_mcv_shp |>
  ggplot() +
  geom_sf(aes(fill = mcv1_12_23), show.legend = TRUE) +
  scale_fill_manual(
    values = brewer.pal(n = 7, "RdYlBu"),
    drop = F
  ) +
  theme_void() +
  theme(
    plot.title = element_text(color = "black", hjust = .5, size = 14),
    legend.key.size = unit(1, "lines"),
    legend.title = element_text(colour = "black", size = 11),
    legend.text = element_text(color = 'black', size = 11)
  ) +
  labs(fill = "MCV1") +
  ggtitle('Actual MCV1')
plt_mcv1_dhs + plt_mcv1_pred

subcountyPred <- subcountyLevel |>
  st_as_sf() |>
  ggplot() +
  geom_sf(aes(fill = mcv1_12_23), show.legend = T) +
  scale_fill_manual(
    values = c(brewer.pal(n = 7, "RdYlBu")),
    drop = F
  ) +
  theme_void() +
  theme(
    plot.title = element_text(color = "black", hjust = .5, size = 17),
    legend.key.size = unit(1, "lines"),
    legend.title = element_text(colour = "black", size = 13),
    legend.text = element_text(color = 'black', size = 13)
  ) +
  labs(fill = "") +
  ggtitle('') +
  guides(col = 'none')
ggsave(
  plot = subcountyPred,
  'plots/subcountyPred.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 5,
  units = "in"
)

# Plotting for threshold --------------------------------------------------

plt_mcv1_thresh95 <- countyLevelThresh |>
  ggplot() +
  geom_sf(aes(fill = thresh95), show.legend = TRUE) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1, limits = c(0, 1)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "") 

ggsave(
  plot = plt_mcv1_thresh95,
  'plots/plt_mcv1_thresh95.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 5,
  units = "in"
)

plt_mcv1_thresh80 <- countyLevelThresh |>
  ggplot() +
  geom_sf(aes(fill = thresh80), show.legend = TRUE) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1, limits = c(0, 1)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "") 

ggsave(
  plot = plt_mcv1_thresh80,
  'plots/plt_mcv1_thresh80.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 5,
  units = "in"
)

# Subcounty level
subplt_mcv1_thresh95 <- subcountyLevelThresh |>
  st_as_sf() |>
  ggplot() +
  geom_sf(aes(fill = thresh95_sub), show.legend = TRUE) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1, limits = c(0, 1)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "") 

ggsave(
  plot = subplt_mcv1_thresh95,
  'plots/subplt_mcv1_thresh95.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 5,
  units = "in"
)

subplt_mcv1_thresh80 <- subcountyLevelThresh |>
  st_as_sf() |>
  ggplot() +
  geom_sf(aes(fill = thresh80_sub), show.legend = TRUE) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1, limits = c(0, 1)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "") 

ggsave(
  plot = subplt_mcv1_thresh80,
  'plots/subplt_mcv1_thresh80.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 5,
  units = "in"
)

# Plotting ridgeline plots ------------------------------------------------

ext <- data.table( cbind(data.table(dhs_shp)[, 1], cc_weighted_all[, -1])) |> 
  _[, melt(.SD, id.vars = "county", value.name = "mcv1")]

county_order <- ext[, .(median = median(mcv1)), by = county][order(median)]
ext[, county := factor(county, levels = county_order$county)]

# Plotting
ridgePlot <- ggplot(ext, aes(x = mcv1, y = county, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 5, rel_min_height = 0.02) +
  scale_x_continuous(expand = c(0, 0), labels = scales::label_percent()) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = rev(c(
    "#1a9641", "#a6d96a", "#fdae61", "#d7191c"
  )), guide = F) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(colour = "black", size = 14)
  ) +
  labs(x = "MCV1(%)", y = "County")

ggsave(
  plot = ridgePlot,
  filename = "plots/ridgePlot.png",
  dpi = 5e2,
  width = 4 * 5,
  height = 3 * 5,
  units = "in",
  bg = NULL
)

# Obtaining some statistics after predicting ------------------------------
# Printing mean prevalence and CI
round(mean(subcountyLevel$mcv1, na.rm = T), 2)
round(quantile(subcountyLevel$mcv1, probs = c(0.025, 0.5, 0.975), na.rm = T), 2)
data.table(subcountyLevel) |> 
  na.omit() |> 
  _[, .N, by = mcv1_12_23]

# Probability of exceedance at the pixel level -----------------------------------------------

#Threshold calculations - 10%, 2%
prob.95 <- apply(inv.linpred, 1, \(.) length(which(. >= 0.95))/nsamp) 
prob95 <- rasterize(x = coord.p, y = grid, field = prob.95)

prob.80 <- apply(inv.linpred, 1, \(.) length(which(. >= 0.80))/nsamp) 
prob80 <- rasterize(x = coord.p, y = grid, field = prob.80)

# Cropping to match DHS shapefile
probs <- list(
  prob95,
  prob80
) |>   
  imap(\(x, y){
    x |> 
      crop(dhs_shp) |> 
      mask(dhs_shp)
  }) |> 
  setNames(c("prob95", "prob80"))


# Plotting the raster files
rast_95 <- ggplot() +
  geom_tile(data = as.data.frame(rast(probs$prob95), xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = dhs_shp,
    fill = NA
  ) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 1)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "") 

ggsave(
  plot = rast_95,
  'plots/rast_95.png',
  dpi = 5e2,
  height = 5,
  width = 5,
  bg = NULL,
  units = "in"
)

rast_80 <- ggplot() +
  geom_tile(data = as.data.frame(rast(probs$prob80), xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = dhs_shp,
    fill = NA
  ) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(0, 1)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  guides(fill = guide_colorbar(keyheight = unit(2, "inches"))) +
  labs(fill = "") 

ggsave(
  plot = rast_80,
  'plots/rast_80.png',
  dpi = 5e2,
  height = 5,
  width = 5,
  bg = NULL,
  units = "in"
)

# Forest plot for sub-counties ---------------------------------------------

mm_c <- merge(data.table(countyLevel)[, .(county, lower_ci, upper_ci,
                                          pred_prev = mcv1, geometry)], 
              data.table(county_mcv_shp)[, .(county, actual_prev = prop)],
              by = c("county")
) |> 
  st_as_sf()
nrow(mm_c)

mm_c |>
  ggplot(aes(x = county)) +
  geom_point(aes(y = pred_prev, colour = 'Predicted')) +
  geom_point(aes(y = actual_prev, colour = "Actual")) +
  geom_errorbar(aes(
    ymin = lower_ci,
    ymax = upper_ci,
    y = pred_prev,
    colour = 'Predicted'
  )) +
  coord_flip() +
  scale_colour_manual(values = c("Predicted" = "#1f78b4", "Actual" = "#e31a1c")) +
  theme_classic()

# Which are within the CI?
mm_c |> 
  data.table() |> 
  dplyr::mutate(
    within_ci = case_when(
      actual_prev >= lower_ci & actual_prev <= upper_ci ~ "yes",
      TRUE ~ "no"
    )
  ) |> 
  ungroup() |> 
  group_by(within_ci) |> 
  dplyr::summarise(n = n())

# Getting the number of unvaccinated children at county level -----------------------------

temp_shp <- dhs_shp |>
  st_transform(crs = st_crs(covs$pop_1))
county_pop1 <-
  data.table(
    county = temp_shp$county,
    pop_1 = raster::extract(covs$pop_1, temp_shp, fun = sum, na.rm = T)[, 1]
  )
un_cc_weighted_all <- apply(cc_weighted_all[, -1], 2, \(z) (1 - z) * county_pop1$pop_1)

un_cc_weighted <- t(apply(un_cc_weighted_all[, -1], 1, \(x) c(
  mean(x), 
  sd(x), 
  quantile(x, probs = c(0.025, 0.5, 0.975))))
) %>%
  apply(., 2, round) |>
  cbind(cc_weighted_all[, 1])

colnames(un_cc_weighted) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "ID")
county_nomcv1 <- cbind(dhs_shp, un_cc_weighted) |> 
  mutate(
    pop_1 = round(county_pop1$pop_1),
    prop_unvacc = mean/pop_1
  )

# Getting the number of unvaccinated children at sub-county level -----------------------------

temp_shp2 <- subcounty_shapefile |>
  st_transform(crs = st_crs(covs$pop_1))

subcounty_pop1 <-
  data.table(
    county = subcounty_shapefile$county,
    subcounty = subcounty_shapefile$subcounty,
    pop_1 = raster::extract(covs$pop_1, temp_shp2, fun = sum, na.rm = T)[, 1]
  )
un_sc_weighted_all <- apply(sc_weighted_all[, -1], 2, \(z) (1 - z) * subcounty_pop1$pop_1)
un_sc_weighted <- t(apply(un_sc_weighted_all[, -1], 1, \(x) c(
  mean(x), 
  sd(x), 
  quantile(x, probs = c(0.025, 0.5, 0.975))))
) %>%
  apply(., 2, round) |>
  cbind(sc_weighted_all[, 1])

colnames(un_sc_weighted) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "ID")
subcounty_nomcv1 <- un_sc_weighted[, ID := as.integer(ID)] |> 
  data.table::merge.data.table(data.table(subcounty_shapefile2), by = 'ID') |> 
  mutate(
    pop_1 = round(subcounty_pop1$pop_1),
    prop_unvacc = mean/pop_1
  )

# Plotting for the unvaccinated at county level -------------------------------------------
# 
# # Continous-number
plot_unvacc_cont <- county_nomcv1 |>
  st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = mean)) +
  scale_fill_distiller(
    palette = "RdYlBu",
    direction = -1,
    labels = scales::label_comma()
  ) +
  annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale() +
  theme_void() +
  theme(
    plot.title = element_text(color = "black", hjust = .5, size = 17),
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    legend.key.size = unit(1, "lines"),
    legend.title = element_text(colour = "black", size = 13),
    legend.text = element_text(color = 'black', size = 13)
  ) +
  labs(fill = "")
ggsave(
  plot = plot_unvacc_cont,
  'plots/plot_unvacc_cont.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 5,
  units = "in"
)
# 
# # Categorical
# plot_unvacc_cat <-unvacc_shp |> 
#   ggplot() +
#   geom_sf(aes(fill = prop_unvacc_cat)) +
#   scale_fill_manual(
#     values = c(
#       "0-9%" = '#4575b4',
#       "10-19%" = "#91bfdb",
#       "20-29%" = "#ffffbf",
#       "30-39%" = '#fee090',
#       "40-49%" = '#fc8d59',
#       "50-59%" = '#d73027',
#       "60-70%" = '#a50026'
#     )
#     
#     ) +
#   annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
#   annotation_scale() +
#   theme_void() +
#   theme(
#     plot.title = element_text(color = "black", hjust = .5, size = 17),
#     # legend.position = "bottom",
#     # legend.direction = "horizontal",
#     legend.key.size = unit(1, "lines"),
#     legend.title = element_text(colour = "black", size = 13),
#     legend.text = element_text(color = 'black', size = 13)
#   ) +
#   labs(fill = "") 
# ggsave(
#   plot = plot_unvacc_cat,
#   'plots/plot_unvacc_cat.png',
#   bg = NULL,
#   dpi = 5e2,
#   height = 5,
#   width = 5,
#   units = "in"
# )
# 
# 
# Plotting for the unvaccinated at subcounty level -------------------------------------------
# # Continous-number
subplot_unvacc_cont <- subcounty_nomcv1 |>
  st_as_sf() |> 
  ggplot() +
  geom_sf(aes(fill = mean)) +
  scale_fill_distiller(
    palette = "RdYlBu",
    labels = scales::label_comma(),
    direction = -1
  ) +
  annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale() +
  theme_void() +
  theme(
    plot.title = element_text(color = "black", hjust = .5, size = 17),
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    legend.key.size = unit(1, "lines"),
    legend.title = element_text(colour = "black", size = 13),
    legend.text = element_text(color = 'black', size = 13)
  ) +
  labs(fill = "")

ggsave(
  plot = subplot_unvacc_cont,
  'plots/subplot_unvacc_cont.png',
  bg = NULL,
  dpi = 5e2,
  height = 5,
  width = 5,
  units = "in"
)
# 
# # Categorical
# subplot_unvacc_cat <-sub_unvacc_shp |> 
#   ggplot() +
#   geom_sf(aes(fill = prop_unvacc_cat)) +
#   scale_fill_manual(
#     values = c(
#       "0-9%" = '#4575b4',
#       "10-19%" = "#91bfdb",
#       "20-29%" = "#ffffbf",
#       "30-39%" = '#fee090',
#       "40-49%" = '#fc8d59',
#       "50-59%" = '#d73027',
#       "60-70%" = '#a50026',
#       "Above 70%" = '#FF4500'
#     )
#     
#     ) +
#   annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
#   annotation_scale() +
#   theme_void() +
#   theme(
#     plot.title = element_text(color = "black", hjust = .5, size = 17),
#     # legend.position = "bottom",
#     # legend.direction = "horizontal",
#     legend.key.size = unit(1, "lines"),
#     legend.title = element_text(colour = "black", size = 13),
#     legend.text = element_text(color = 'black', size = 13)
#   ) +
#   labs(fill = "") 
# ggsave(
#   plot = subplot_unvacc_cat,
#   'plots/subplot_unvacc_cat.png',
#   bg = NULL,
#   dpi = 5e2,
#   height = 5,
#   width = 5,
#   units = "in"
# )

#beepr::beep(4)

# Modified non-randomized probability integral transform --------------------------------------

# Sourcing files and logging
# source("07_anpit_functions.R")
# log_file <- "logs/nrpit_mcv1.log"
# dir.create("logs", showWarnings = FALSE)

# Reps
reps <- 10

# Running the model
# nrpit_results <- map(1:reps, \(i) {
#   set.seed(i * 42)
#   nrpit_results <- validate(
#     dm2 = dm2,
#     kenya_shapefile2 = kenya_shapefile2,
#     p = 0.3,
#     u_grid = seq(0, 1, by = 0.1),
#     n_sim = 1000,
#     n_env_sim = 1000,
#     log_file = log_file
#   )
# })

nrpit_results <- readRDS("data/nrpit_mcv1.RDS")

# Averaging the results
all_anpit <- imap(nrpit_results, \(x, i) {
  x$anpit[, rep := i]
}) |> rbindlist()

# Averaged AnPIT + averaged envelope
anpit_avg <- all_anpit[, .(
  anpit = mean(anpit),
  lower = mean(lower),
  upper = mean(upper)
), by = u]

# Mean CPRS
mean_crps <- mean(map_dbl(nrpit_results, \(x) x$anpit$crps[1]))

# Saving the results
# saveRDS(nrpit_results, "data/nrpit_mcv1.RDS")

# Plotting for AnPIT --------------------------------------------------------------------------

average_anpit_plot <- anpit_avg |> 
  ggplot(aes(x = u, y = anpit)) +
  geom_abline(slope = 1, intercept = 0, color = "#d6604d", linewidth = 0.5) +
  geom_line(aes(y = lower), color = "black", lty = 2) +
  geom_line(linewidth = 0.8, color = "#2166ac") +
  geom_line(aes(y = upper), color = "black", lty = 2) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    x = "u",
    y = "AnPIT(u)",
    caption = "AnPIT - Average non-randomized Probability Integral Transform"
  ) +
  theme_minimal(base_size = 11) +
  coord_equal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position   = "bottom",
    strip.text = element_text(color = "black")
  )

ggsave(
  plot = average_anpit_plot,
  filename = "plots/average_anpit_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

ind_anpit_plot <- all_anpit |>
  ggplot(aes(x = u, y = anpit, group = rep)) +
  geom_abline(slope = 1, intercept = 0, color = "#d6604d", linewidth = 0.5) +
  geom_line(color = "#2166ac", alpha = 0.6, linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "u", y = "AnPIT(u)") +
  theme_minimal() +
  coord_equal()

ggsave(
  plot = ind_anpit_plot,
  filename = "plots/ind_anpit_plot.png",
  dpi = 5e2,
  width = 1 * 5,
  height = 1 * 5,
  units = "in",
  bg = "white"
)

# Saving work space ---------------------------------------------------------------------------

save.image("workspaces/inla_mcv1.RData", safe = T)


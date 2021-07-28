# Convert 2000 geography characteristics to 2010 equivalents and interpolate

library(doSNOW)
library(parallel)

library(data.table)
library(progress)

cores <- 40

# Load data ---------------------------------------------------------------

nhgis_2000 <- fread("/media/qnap3/Covariates/nhgis/output/2000_zcta.csv.gz"
                    )[, `:=` (geoid = sprintf("%05.f", as.numeric(GEOID)),
                              source = "census",
                              weight = 4)
                      ][, `:=` (med_household_income = median_household_income,
                                pct_non_hispanic_white = pct_nonhispanic_white,
                                pct_non_hispanic_black = pct_nonhispanic_black,
                                pct_non_hispanic_asian = pct_nonhispanic_asian)]

nhgis_2010 <- fread("/media/qnap3/Covariates/nhgis/output/2010_zcta.csv.gz"
                    )[, `:=` (geoid = sprintf("%05.f", as.numeric(GEOID)),
                              source = "census",
                              weight = 4)
                      ][, `:=` (pct_non_hispanic_white = pct_nonhispanic_white,
                                pct_non_hispanic_black = pct_nonhispanic_black,
                                pct_non_hispanic_asian = pct_nonhispanic_asian)]

acs5 <- fread("/media/qnap3/Covariates/census-ses-covariates/outputs/tables-us/time_series_2009_to_2019_zip_codes_long.csv.gz"
              )[, `:=` (geoid = substr(GEOID, 8, 8+5),
                        source = "acs5",
                        weight = 2)
                ][, `:=` (population = population.x,
                          houses = n_housing_units)]

tdw_pass1 <- fread("output/tdw1-bd-areas.csv"
                   )[, `:=` (ZCTA5CE00 = sprintf("%05.f", as.numeric(ZCTA5CE00)),
                             ZCTA5CE10 = sprintf("%05.f", as.numeric(ZCTA5CE10)),
                             tdw1 = tdw)]

tdw_pass2 <- fread("output/tdw2-buffer-areas.csv",
                   )[, `:=` (ZCTA5CE00 = sprintf("%05.f", as.numeric(ZCTA5CE00)),
                             ZCTA5CE10 = sprintf("%05.f", as.numeric(ZCTA5CE10)),
                             tdw2 = tdw)]

tdw_pass3 <- fread("output/tdw3-raw-areas.csv",
                   )[, `:=` (ZCTA5CE00 = sprintf("%05.f", as.numeric(ZCTA5CE00)),
                             ZCTA5CE10 = sprintf("%05.f", as.numeric(ZCTA5CE10)),
                             tdw3 = tdw)]

matched <- fread("output/matches.csv")[, `:=` (ZCTA5CE00 = sprintf("%05.f", as.numeric(ZCTA5CE00)),
                                               ZCTA5CE10 = sprintf("%05.f", as.numeric(ZCTA5CE10)))]

# Merge TDWs into one data.table
tdw <- Reduce(
  function(left, right) return(left[right, on = list(ZCTA5CE10, ZCTA5CE00)]),
  lapply(
    list(tdw_pass1, tdw_pass2, tdw_pass3),
    function(dt) return(dt[, !"tdw"])
  )
)

# 2000 geographies with missing 2010 geographies: pull the population from 2010
tdw <- rbindlist(
  list(
    tdw,
    nhgis_2010[matched[missing == "ZCTA5CE10"][, geoid := ZCTA5CE10],
               on = "geoid"
               ][, .(ZCTA5CE00,
                     ZCTA5CE10,
                     matched_population = population,
                     matched_houses = houses)
                 ]
  ),
  fill = TRUE
)
  
# 2010 geographies with missing 2000 geographies: leave a marker indicating that
# we will later pull the population directly
tdw <- rbindlist(
  list(
    tdw,
    matched[missing == "ZCTA5CE00", .(ZCTA5CE00,
                                      ZCTA5CE10,
                                      matched_2010_to_2000 = TRUE)]
  ),
  fill = TRUE
)
tdw[, matched_2010_to_2000 := ifelse(is.na(matched_2010_to_2000), FALSE, matched_2010_to_2000)]

# Convert 2000 Census to 2010 geography equivalent ------------------------

imputable_vars <- setdiff(
  intersect(names(nhgis_2000), names(acs5)),
  c("source", "houses", "weight", "population", "GEOID", "geoid")
)
imputable_population_vars <- imputable_vars[
  !grepl("housing", imputable_vars) & !grepl("heating", imputable_vars)
]
imputable_housing_vars <- imputable_vars[
  (grepl("housing", imputable_vars) | grepl("heating", imputable_vars))
]

## Initial setup ----

# Join TDWs
tdw_joined <- nhgis_2000[tdw, on = c("geoid" = "ZCTA5CE00")]

# Join information about which pass TDW to use
tdw_joined <- tdw[, .(tdw1_ok = !any(is.na(tdw1)),
                      tdw2_ok = !any(is.na(tdw2)),
                      tdw3_ok = !any(is.na(tdw3))),
                  by = "ZCTA5CE10"
                  ][tdw_joined, on = "ZCTA5CE10"
                    ][, tdw := fcase(tdw1_ok == TRUE, tdw1,
                                     tdw2_ok == TRUE, tdw2,
                                     tdw3_ok == TRUE, tdw3)]

# Calculate population and houses by TDW
nhgis_2000_2010equiv <- tdw_joined[, .(# Population
                                       population = first(fcase(
                                         matched_2010_to_2000 == TRUE, as.double(population),
                                         !is.na(matched_population), as.double(matched_population),
                                         tdw1_ok == TRUE, sum(population * tdw1),
                                         tdw2_ok == TRUE, sum(population * tdw2),
                                         tdw3_ok == TRUE, sum(population * tdw3)
                                       )),
                                       
                                       # Houses
                                       houses = first(fcase(
                                         matched_2010_to_2000 == TRUE, as.double(houses),
                                         !is.na(matched_houses), as.double(matched_houses),
                                         tdw1_ok == TRUE, sum(houses * tdw1),
                                         tdw2_ok == TRUE, sum(houses * tdw2),
                                         tdw3_ok == TRUE, sum(houses * tdw3)
                                       )),
                                       
                                       # Interpolation method
                                       method = as.factor(first(fcase(
                                         matched_2010_to_2000 == TRUE, "matched2010to2000",
                                         !is.na(matched_houses), "matched2000to2010",
                                         tdw1_ok == TRUE, "TDWpass1",
                                         tdw2_ok == TRUE, "TDWpass2",
                                         tdw3_ok == TRUE, "TDWpass3"
                                       )))),
                                   by = "ZCTA5CE10"]

### Diagnostics ----
if (FALSE) {
  library(dplyr)
  library(mapview)
  
  tiger_2010 <- st_read("output/zctas/zcta5_2000.gpkg")
  
  table(nhgis_2000_2010equiv$method, exclude = "no")
  round(table(nhgis_2000_2010equiv$method, exclude = "no") / nrow(nhgis_2000_2010equiv) * 100, 2)
  
  # Where / what are the missing geographies?
  missing_geographies <- tiger_2010 %>%
    filter(ZCTA5CE10 %in% nhgis_2000_2010equiv[is.na(method)]$ZCTA5CE10) %>%
    left_join(
      transmute(nhgis_2010, population, ZCTA5CE10 = geoid),
      by = "ZCTA5CE10"
    )
  
  missing_geographies %>%
    select(ZCTA5CE10, population) %>%
    st_drop_geometry()
  
  mapview(missing_geographies)
}

## Interpolate population-related variables ----

nhgis_2010_matched <- nhgis_2010[, ZCTA5CE10 := geoid
                                 ][nhgis_2000_2010equiv[method == "matched2010to2000"], on = "ZCTA5CE10"]
acs_only_variables <- c("ZCTA5CE10", setdiff(imputable_vars, names(nhgis_2010_matched)))
nhgis_2010_matched <- acs5[, ZCTA5CE10 := geoid
                           ][year == 2011, ..acs_only_variables # FIXME: we only had 2011-2019 data, 2009-2010 was missing
                            ][nhgis_2010_matched, on = "ZCTA5CE10"]

bar <- progress_bar$new(
  "Interpolating 2000 -> 2010 Census population variables :current/:total (:percent) [:bar] eta :eta",
  total = length(imputable_population_vars)
)
for (variable in imputable_population_vars) {
  tdw_interpolated <- rbindlist(list(
    # Original TDW: hat y_t = sum_s(w_s * y_s)
    # TDW with precalculated weights: hat y_t = sum_s(w_s * y_s * z_s) / sum(w_s * z_s)
    tdw_joined[matched_2010_to_2000 == FALSE,
               setNames(
                 list(
                   # Need to convert integers to doubles to avoid overflow
                   sum(as.double(get(variable)) * population * tdw, na.rm = TRUE) / sum(population * tdw, na.rm = TRUE)
                 ),
                 variable
               ),
               by = "ZCTA5CE10"],
    
    # For matched 2010 -> 2000, just pull directly from 2010
    nhgis_2010_matched[, setNames(list(ZCTA5CE10, get(variable)), c("ZCTA5CE10", variable))]
  ))
  nhgis_2000_2010equiv <- nhgis_2000_2010equiv[tdw_interpolated, on = "ZCTA5CE10"]
  bar$tick()
}

## Interpolate housing-related variables ----
bar <- progress_bar$new(
  "Interpolating 2000 -> 2010 Census housing variables :current/:total (:percent) [:bar] eta :eta",
  total = length(imputable_housing_vars)
)
for (variable in imputable_housing_vars) {
  tdw_interpolated <- rbindlist(list(
    # Original TDW: hat y_t = sum_s(w_s * y_s)
    # TDW with precalculated weights: hat y_t = sum_s(w_s * y_s * z_s) / sum(w_s * z_s)
    tdw_joined[matched_2010_to_2000 == FALSE,
               setNames(
                 list(
                   # Need to convert integers to doubles to avoid overflow
                   sum(as.double(get(variable)) * houses * tdw, na.rm = TRUE) / sum(houses * tdw, na.rm = TRUE)
                 ),
                 variable
               ),
               by = "ZCTA5CE10"],
    
    # For matched 2010 -> 2000, just pull directly from 2010
    nhgis_2010_matched[, setNames(list(ZCTA5CE10, get(variable)), c("ZCTA5CE10", variable))]
  ))
  nhgis_2000_2010equiv <- nhgis_2000_2010equiv[tdw_interpolated, on = "ZCTA5CE10"]
  bar$tick()
}

## Diagnostics ----

# Unreachable code so doesn't run in non-interactive mode
if (FALSE) {
  # Check variables in ZCTAs that did not change much between the years
  test_zcta <- "02115" # Longwood + Northeastern
  round(setNames(
    as.data.frame(cbind(
      t(nhgis_2000[geoid == test_zcta, ..imputable_vars]),
      t(nhgis_2000_2010equiv[ZCTA5CE10 == test_zcta, ..imputable_vars])
    )),
    c("2000", "2010 interpolated")
  ), 10)
}

# Linear interpolation ----------------------------------------------------

## Merge data frames ----

# Last-minute normalization
setnames(nhgis_2000_2010equiv, "ZCTA5CE10", "geoid")
nhgis_2000_2010equiv[, `:=` (year = 2000,
                             source = "census",
                             weight = 4)]

all_columns <- setdiff(names(nhgis_2000_2010equiv), "method")
combined <- rbindlist(list(nhgis_2000_2010equiv, acs5[, ..all_columns]), fill = TRUE)

# For some reason, this gives us an error about the geoid index being invalid,
# so we need to rebuild the data.table
nrow(combined)
combined <- as.data.table(as.data.frame(combined))
nrow(combined)

## Interpolation ----

variables_to_interpolate <- c("population", "houses", imputable_vars)
interpolate_census <- function(target_geoid) {
  input <- combined[geoid == target_geoid]
  result <- data.table(
    geoid = target_geoid,
    year = 2000:2019
  )
  for (variable in variables_to_interpolate) {
    if (all(is.na(input[[variable]]))) {
      result[[variable]] <- NA
    } else {
      model <- lm(
        as.formula(sprintf("%s ~ year", variable)),
        data = input,
        weights = weight
      )
      result[[variable]] <- predict(model, newdata = result)
    }
  }
  return(result)
}

geoids <- sort(unique(combined$geoid))
bar <- progress_bar$new(
  "Interpolating Census variables :current/:total (:percent) [:bar] eta :eta",
  total = length(geoids)
)

# Single-core implementation
# interpolated <- rbindlist(lapply(
#   geoids,
#   function(geoid) {
#     bar$tick()
#     return(interpolate_census(geoid))
#   }
# ))

# Parallel implementation
message(sprintf("Instantiating cluster with %s cores", cores))
cluster <- makeSOCKcluster(cores)
registerDoSNOW(cluster)
on.exit(stopCluster(cluster))
interpolated <- foreach(
  geoid = geoids,
  .packages = c("progress", "data.table"),
  .options.snow = list(progress = function(n) bar$tick())
) %dopar% {
  return(interpolate_census(geoid))
}
stopCluster(cluster)
interpolated <- rbindlist(interpolated)[, source := "imputed"]
fwrite(interpolated, "output/interpolated.csv")

## Diagnostics ----

library(ggplot2)

# Unreachable code so doesn't run in non-interactive mode
if (FALSE) {
  # National correlations
  # DT[i, j, by] -> "combined" are prefixed with "i."
  comparison <- interpolated[combined, on = c("geoid", "year")]
  correlations_df <- data.table(
    variable = variables_to_interpolate,
    overall = sapply(
      variables_to_interpolate,
      function(variable) {
        variable_source <- paste0("i.", variable)
        return(cor(
          comparison[[variable]],
          comparison[[variable_source]],
          use = "complete.obs"
        ))
      }
    )
  )
  for (target_year in unique(comparison$year)) {
    comparison_subset <- comparison[year == target_year]
    temp <- data.table(
      variable = variables_to_interpolate,
      target_year = sapply(
        variables_to_interpolate,
        function(variable) {
          variable_source <- paste0("i.", variable)
          return(cor(
            comparison_subset[[variable]],
            comparison_subset[[variable_source]],
            use = "complete.obs"
          ))
        }
      )
    )
    setnames(temp, "target_year", as.character(target_year))
    correlations_df <- correlations_df[temp, on = "variable"]
  }
  View(correlations_df[, lapply(.SD, round, digits = 4), by = variable])
  fwrite(correlations_df, "diagnostics/overall_correlations.csv")
  
  # National scatter plot
  yearly_means_interpolated <- interpolated[, !c("geoid", "source")
                                            ][, lapply(.SD, mean, na.rm = TRUE), by = "year"
                                              ][, source := "imputed"]
  yearly_means <- combined[, !c("geoid", "source")
                           ][, lapply(.SD, mean, na.rm = TRUE), by = "year"
                             ][, source := ifelse(year == 2000, "census", "acs5")]
  ggplot() +
    aes(x = year, y = pct_heating_solar, color = source) +
    geom_point(data = yearly_means_interpolated) +
    geom_point(data = yearly_means)
  
  # Single-ZCTA
  target_geoid <- "02115"
  target_geoid <- sample(unique(comparison$geoid), 1)
  comparison_subset <- comparison[geoid == target_geoid]
  correlations_df_onegeoid <- data.table(
    variable = variables_to_interpolate,
    correlation = lapply(
      variables_to_interpolate,
      function(variable) {
        variable_source <- paste0("i.", variable)
        return(cor(
          comparison_subset[[variable]],
          comparison_subset[[variable_source]],
          use = "complete.obs"
        ))
      }
    )
  )
  correlations_df_onegeoid[]
  fwrite(correlations_df_onegeoid, sprintf("diagnostics/correlations_%s.csv", target_geoid))
  ggplot() +
    aes(x = year, y = population, color = source) +
    geom_point(data = combined[geoid == target_geoid]) +
    geom_point(data = interpolated[geoid == target_geoid])
  ggplot() +
    aes(x = year, y = pct_female, color = source) +
    geom_point(data = combined[geoid == target_geoid]) +
    geom_point(data = interpolated[geoid == target_geoid])
  ggplot() +
    aes(x = year, y = pct_black, color = source) +
    geom_point(data = combined[geoid == target_geoid]) +
    geom_point(data = interpolated[geoid == target_geoid])
  ggplot() +
    aes(x = year, y = pct_travel_60_to_89_min, color = source) +
    geom_point(data = combined[geoid == target_geoid]) +
    geom_point(data = interpolated[geoid == target_geoid])
}
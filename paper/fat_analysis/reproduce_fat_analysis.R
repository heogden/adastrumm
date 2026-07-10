# Reproduce the figures from the body-fat data analysis in the AdaStruMM paper.
# This script should be run from the root of the adastrumm repository

suppressPackageStartupMessages({
  library(adastrumm)
  library(dplyr)
  library(ggplot2)
  library(nlme)
  library(tidyr)
  library(tikzDevice)
})


analysis_directory <- file.path("paper", "fat_analysis")
figure_directory <- file.path(analysis_directory, "figures")

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

figure_path <- function(filename) {
  file.path(figure_directory, filename)
}

write_tikz <- function(plot, filename, width, height) {
  tikzDevice::tikz(
    file = figure_path(filename),
    width = width,
    height = height
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(NULL)
}

# Load and prepare the data -------------------------------------------------

data("fat", package = "ALA")

fat_data <- fat %>%
  as_tibble() %>%
  transmute(
    c = as.integer(as.factor(id)),
    x = time.menarche,
    y = percent.fat,
    bf4 = x * (x > 0)
  )

# Fit the comparison model from Section 8.8 of Fitzmaurice et al.
fit_fitzmaurice <- nlme::lme(
  fixed = y ~ x + bf4,
  random = ~ x + bf4 | c,
  data = fat_data
)

# Fit the AdaStruMM model. The fitting procedure is deterministic for fixed
# data, package versions and numerical libraries, so no random seed is needed
# before this call.
fit_ada <- adastrumm::fit_adastrumm(fat_data, nbasis = 20)

# Predictions and uncertainty samples --------------------------------------

x_grid <- seq(
  from = min(fat_data$x),
  to = max(fat_data$x),
  length.out = 100L
)
cluster_ids <- sort(unique(fat_data$c))

prediction_data <- tidyr::crossing(
  x = x_grid,
  c = cluster_ids
) %>%
  arrange(x, c)

n_samples <- 1000L

# find_samples() uses a fixed seed for each sample, so these samples are
# reproducible for a fixed fitted model.
samples <- adastrumm::find_samples(fit_ada, n_samples)

predict_for_samples <- function(samples, model, newdata, derivative = FALSE) {
  vapply(
    samples,
    adastrumm::predict_y_given_sample,
    FUN.VALUE = numeric(nrow(newdata)),
    mod = model,
    newdata = newdata,
    deriv = derivative
  )
}

mean_samples <- predict_for_samples(
  samples = samples,
  model = fit_ada,
  newdata = prediction_data
)

derivative_samples <- predict_for_samples(
  samples = samples,
  model = fit_ada,
  newdata = prediction_data,
  derivative = TRUE
)

subject_predictions <- prediction_data %>%
  mutate(
    estimate = adastrumm::predict_adastrumm(
      fit_ada,
      newdata = prediction_data
    ),
    lower = apply(mean_samples, 1L, stats::quantile, probs = 0.025),
    upper = apply(mean_samples, 1L, stats::quantile, probs = 0.975)
  )

# Subject-specific fitted curves -------------------------------------------

clusters_to_plot <- head(cluster_ids, 20L)
fat_subset <- fat_data %>%
  filter(c %in% clusters_to_plot)

plot_fitted_curves <- subject_predictions %>%
  filter(c %in% clusters_to_plot) %>%
  ggplot(aes(x = x)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2
  ) +
  geom_line(aes(y = estimate)) +
  geom_point(aes(y = y), data = fat_subset) +
  facet_wrap(vars(c)) +
  coord_cartesian(xlim = range(fat_subset$x)) +
  xlab("Time relative to menarche (years)") +
  ylab("Percent body fat")

write_tikz(
  plot = plot_fitted_curves,
  filename = "fat_fitted_curves.tex",
  width = 6,
  height = 3.5
)

# Population-average trajectory --------------------------------------------

n_clusters <- length(cluster_ids)
n_x <- length(x_grid)

# prediction_data is ordered with clusters varying fastest within each x.
mean_sample_array <- array(
  mean_samples,
  dim = c(n_clusters, n_x, n_samples)
)

population_mean_samples <- apply(
  mean_sample_array,
  MARGIN = c(2L, 3L),
  FUN = mean
)

population_mean_data <- tibble(
  x = x_grid,
  estimate = rowMeans(population_mean_samples),
  lower = apply(
    population_mean_samples,
    1L,
    stats::quantile,
    probs = 0.025
  ),
  upper = apply(
    population_mean_samples,
    1L,
    stats::quantile,
    probs = 0.975
  )
)

fitzmaurice_coefficients <- nlme::fixef(fit_fitzmaurice)
fitzmaurice_intercept <- unname(fitzmaurice_coefficients[["(Intercept)"]])
fitzmaurice_slope_before_zero <- unname(fitzmaurice_coefficients[["x"]])
fitzmaurice_slope_after_zero <-
  fitzmaurice_slope_before_zero +
  unname(fitzmaurice_coefficients[["bf4"]])

plot_population_mean <- population_mean_data %>%
  ggplot(aes(x = x)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2
  ) +
  geom_line(aes(y = estimate)) +
  annotate(
    "segment",
    x = -4,
    xend = 0,
    y = fitzmaurice_intercept - 4 * fitzmaurice_slope_before_zero,
    yend = fitzmaurice_intercept,
    linetype = "dashed"
  ) +
  annotate(
    "segment",
    x = 0,
    xend = 4,
    y = fitzmaurice_intercept,
    yend = fitzmaurice_intercept + 4 * fitzmaurice_slope_after_zero,
    linetype = "dashed"
  ) +
  coord_cartesian(xlim = c(-4, 4)) +
  xlab("Time relative to menarche (years)") +
  ylab("Average percent body fat")

write_tikz(
  plot = plot_population_mean,
  filename = "fat_pa.tex",
  width = 4,
  height = 3
)

# Distribution of fitted subject-specific derivatives ----------------------

derivative_sample_array <- array(
  derivative_samples,
  dim = c(n_clusters, n_x, n_samples)
)

# Average over uncertainty samples for each subject and value of x, then
# summarise the fitted derivatives across subjects.
mean_derivative_by_subject <- apply(
  derivative_sample_array,
  MARGIN = c(1L, 2L),
  FUN = mean
)

derivative_quantiles <- apply(
  mean_derivative_by_subject,
  MARGIN = 2L,
  FUN = stats::quantile,
  probs = c(0.05, 0.25, 0.50, 0.75, 0.95)
)

derivative_data <- tibble(
  x = x_grid,
  q5 = derivative_quantiles[1L, ],
  q25 = derivative_quantiles[2L, ],
  estimate = derivative_quantiles[3L, ],
  q75 = derivative_quantiles[4L, ],
  q95 = derivative_quantiles[5L, ]
)

plot_derivatives <- derivative_data %>%
  filter(x >= -4, x <= 4) %>%
  ggplot(aes(x = x)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.4) +
  geom_line(aes(y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Time relative to menarche (years)") +
  ylab("Rate of change of percent body fat")

write_tikz(
  plot = plot_derivatives,
  filename = "fat_estimated_derivs.tex",
  width = 4,
  height = 3
)

message("Figures written to: ", figure_directory)

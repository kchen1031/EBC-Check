library(ebcalibration)

bundle <- get_data("narasimhan", "surgery")

result <- eb_calibration_check(
  X_obs              = bundle$data$s,
  prior_sampler      = bundle$prior_sampler,
  likelihood_sampler = bundle$likelihood_sampler,
  feature_transform  = bundle$feature_transform,
  B = 500, seed = 1
)
print(result)
plot(result)
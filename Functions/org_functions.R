

gen_trials <- function(blocks=3,numinputs=3,shuffle=FALSE) {
  if (shuffle) {
    examples <- as.vector(apply(replicate(blocks,seq(1, numinputs)), 2,
                                sample, numinputs))
  } else{
    examples <- as.vector(replicate(blocks, seq(1, numinputs)))
  }
}



normal_predictive_distribution <- function(mu_samples,
                                           sigma_samples,
                                           N_obs) {
  map2_dfr(mu_samples, sigma_samples, function(mu, sigma) {
    tibble(
      trialn = seq_len(N_obs),
      t_pred = rnorm(N_obs, mu, sigma)
    )
  }, .id = "iter") %>%
    # .id is always a string and
    # needs to be converted to a number
    mutate(iter = as.numeric(iter))
}

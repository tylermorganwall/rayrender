# Independent reference: per-step RGB mixture free flights with implicit capture,
# and a separate single-channel analog walk. No production proposal code is used.
sss_reference_slab = function(sigma_a, sigma_s, n = 2000L, seed = 194L) {
  set.seed(seed)
  sigma_a = rep(sigma_a, length.out = 3L)
  sigma_s = rep(sigma_s, length.out = 3L)
  extinction = sigma_a + sigma_s
  result = matrix(0, n, 6L)
  for (i in seq_len(n)) {
    z = 0
    mu = 1
    weight = rep(1, 3L)
    repeat {
      endpoint = if (mu > 0) (1 - z) / mu else -z / mu
      channel = sample.int(3L, 1L)
      rate = extinction[channel]
      flight = if (rate == 0) Inf else rexp(1L, rate)
      collision = flight < endpoint
      distance = if (collision) flight else endpoint
      transmission = exp(-extinction * distance)
      if (!collision) {
        weight = weight * transmission / mean(transmission)
        offset = if (mu > 0) 3L else 0L
        result[i, offset + seq_len(3L)] = weight
        break
      }
      density = mean(extinction * transmission)
      weight = weight * transmission * sigma_s / density
      if (!any(weight > 0)) {
        break
      }
      z = z + mu * flight
      mu = runif(1L, -1, 1)
      # This reference uses its own roulette, independent of production depth.
      survival = min(0.99, max(weight))
      if (runif(1L) >= survival) {
        break
      }
      weight = weight / survival
    }
  }
  list(mean = colMeans(result), se = apply(result, 2, stats::sd) / sqrt(n))
}

sss_analog_slab = function(absorption, scattering, n = 3000L, seed = 912L) {
  set.seed(seed)
  result = integer(n)
  extinction = absorption + scattering
  for (i in seq_len(n)) {
    z = 0
    mu = 1
    repeat {
      endpoint = if (mu > 0) (1 - z) / mu else -z / mu
      flight = if (extinction == 0) Inf else rexp(1L, extinction)
      if (flight >= endpoint) {
        result[i] = if (mu > 0) 2L else 1L
        break
      }
      if (runif(1L) < absorption / extinction) {
        break
      }
      z = z + mu * flight
      mu = runif(1L, -1, 1)
    }
  }
  probabilities = tabulate(result + 1L, nbins = 3L) / n
  list(mean = probabilities, se = sqrt(probabilities * (1 - probabilities) / n))
}

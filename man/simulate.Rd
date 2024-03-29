% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulate}
\alias{simulate}
\title{Simulate a utilization distribution and track from a Langevin diffusion model}
\usage{
simulate(
  xlim = c(-2, 2),
  ylim = c(-2, 2),
  boundary_sharpness = 2.5,
  boundary_limit = 0.8,
  cv_code = 1,
  cv_pars = c(2, 0.5, 6.5),
  pred_loc_delta = 0.1,
  nt = 100,
  nping = 40,
  gamma = 0.08,
  ping_tau = 0.01 * c(1, 1),
  ping_cor = 0.3,
  loc_class_probs = c(G = 0.026, `3` = 0.04, `2` = 0.035, `1` = 0.02, `0` = 0.065, A =
    0.126, B = 0.688),
  seed
)
}
\arguments{
\item{xlim, }{ylim The approximate bounding box for movement}

\item{boundary_sharpness}{How penalized should movement outside the bounding box be?}

\item{boundary_limit}{A scaling factor applied to xlim and ylim to determine where the penalty be applied.}

\item{cv_code}{0 = exponential (don't use), 1 = gaussian, 2 = matern, 3 = matern32}

\item{cv_pars}{c(marginal std. dev., range, smoothness)}

\item{pred_loc_delta}{Raster cell size for interpolated utilization distribution raster. If equal to NA, predictions won't be made.}

\item{nt}{How many true movement locations should there be?}

\item{nping}{How many location pings should there be?}

\item{gamma}{Speed parameter.}

\item{ping_tau}{Std. dev.s for ping observation error.}

\item{ping_cor}{Correlation between ping observation error coordinates.}

\item{loc_class_probs}{Named probability factor giving location quality class frequencies.}

\item{seed}{Optional simulation seed.}
}
\value{
A list:
\itemize{
\item nn_graph The output of make_nn_graph
\item field A stars object with the simulated utilization distribution
\item pred_field An upscaled / interpolated copy of field
\item track An sf data.frame with the true location track and field gradient values
\item pings An sf data.frame with the observation pings
\item tmap A list of tmap objects for the above simulated values
\item data A copy of data passed to TMB
\item para A copy of parameters passed to TMB
}
}
\description{
Simulate a utilization distribution and track from a Langevin diffusion model
}

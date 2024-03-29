% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting_subroutines.R
\name{fit_utilization_distribution}
\alias{fit_utilization_distribution}
\title{Fit a utilization distribution using a filtered movement track estimate}
\usage{
fit_utilization_distribution(
  filtered_locations,
  cv_code = 1,
  max.edge = 1,
  ...
)
}
\arguments{
\item{filtered_locations}{The output of fit_rw, e.g. a list
\itemize{
\item pings An sf data.frame with columns
\itemize{
\item t Time for the location ping
\item q Location quality class for the location ping
\item geom Point geometries giving the observed location pings
}
\item track An sf data.frame with columns
\itemize{
\item x_se, y_se Standard errors for the x and y coordinate, respectively
\item t Time for the location estimate
\item geom Point geometries giving the estimated locations
}
\item parameters parameter estimates
}}

\item{cv_code}{0 = exponential (don't use), 1 = gaussian, 2 = matern, 3 = matern32}

\item{max.edge}{The maximum edge length for the mesh}

\item{...}{Additional arguments to pass to make_starve_graph}
}
\value{
A list with the following elements:
\itemize{
\item filtered_locations: A copy of the filtered_locations argument
\item graph: The mesh used for the random field
\item track_graph: The parents used to connect the track and field
\item opt: The output of nlminb
\item sdr: The output of sdreport
\item cv_code: The covariance function code
\item mesh_predictions: A data.frame containing the field predictions for the mesh.
}
}
\description{
Fit a utilization distribution using a filtered movement track estimate
}

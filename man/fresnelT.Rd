% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gisaxs_model.R
\name{fresnelT}
\alias{fresnelT}
\title{Fresnel transmission coeficent for rough interface}
\usage{
fresnelT(ai, delta, beta, surfRms = 0, lambda = 0.154)
}
\arguments{
\item{ai}{angle of incidence in [rad.]}

\item{delta}{Re part of index of refraction}

\item{beta}{Im part of index of refraction}

\item{surfRms}{RMS surface roughness in [nm]}

\item{lambda}{X-ray wavelength}
}
\value{
ft trnamsission coeficient
}
\description{
calculate Fresnel coeficient for transmission at rough surface
}
\examples{
fresnelT(ai = 0.5, delta = 7e-5, beta = 5e-7, surfRms = 0)
}
\author{
K. Juraic
}

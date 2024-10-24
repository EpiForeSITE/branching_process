% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pGenSwitch1.R
\name{pGenSwitch1}
\alias{pGenSwitch1}
\title{Probability that one initial case leads to an outbreak lasting
less than g generations of transmission.}
\usage{
pGenSwitch1(gMax, R0, k0, Rc, kc)
}
\arguments{
\item{gMax}{maximum number of generations}

\item{R0}{basic reproduction number: mean of negative binomial offspring distribution from generation one}

\item{k0}{dispersion of negative binomial offspring distribution from generation one}

\item{Rc}{control reproduction number: mean of negative binomial offspring distribution from generation two plus}

\item{kc}{dispersion of negative binomial offspring distribution from generation two plus}
}
\description{
Probability that one initial case leads to an outbreak lasting
less than g generations of transmission.
}
\author{
Damon Toth
}

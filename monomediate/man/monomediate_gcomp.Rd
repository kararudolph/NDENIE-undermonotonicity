% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monomediate_gcomp.R
\name{monomediate_gcomp}
\alias{monomediate_gcomp}
\title{Parametric mediation estimator under monotonicity assumption}
\usage{
monomediate_gcomp(data, npsem, y_form, z_form, m_form, y_family)
}
\arguments{
\item{data}{[\code{data.frame}]\cr
A \code{data.frame} containing all necessary variables
for the estimation problem.}

\item{npsem}{[\code{R6(Npsem)}]\cr
An \code{Npsem} object}

\item{y_form}{[\code{formula}]\cr
A formula for the outcome model}

\item{z_form}{[\code{formula}]\cr
A formula for the "Z" model}

\item{m_form}{[\code{formula}]\cr
A \bold{RIGHT-HAND SIDE ONLY} formula for the mediator model}

\item{y_family}{[\code{character(1)}]\cr
Outcome variable type (i.e., "gaussian", "binomial").}
}
\value{
Natural (in)direct effect estimates
}
\description{
Parametric mediation estimator under monotonicity assumption
}
\examples{
npsem <- Npsem$new(c("w1", "w2", "w3"), "a", "z", c("m1", "m2"), "y")
y_form <- y ~ w1 + w2 + w3 + z + m1 + m2
z_form <- z ~ w1 + w2 + w3 + a
m_form <- ~ w3 + z
monomediate_gcomp(multiple_m, npsem, y_form, z_form, m_form, "binomial")
}

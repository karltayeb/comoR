% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twococomo.R
\name{pointnormalsusie}
\alias{pointnormalsusie}
\title{Point Normal SuSiE
Fits Point Normal SuSiE, observations are drawn from a two component mixture
(point-mass at zero, normal).
The prior log odds of being non-zero are a linear function of covariates X
and the effects have a SuSiE prior}
\usage{
pointnormalsusie(
  X,
  betahat,
  se,
  Z = NULL,
  L = min(10, ncol(X)),
  f1_params = list(),
  estimate_f1 = FALSE,
  prior_mean = 0,
  prior_variance = 1,
  prior_weights = NULL,
  intercept = TRUE,
  estimate_prior_variance = TRUE,
  s_init = NULL,
  coverage = 0.95,
  min_abs_corr = 0.5,
  max_iter = 100,
  tol = 0.001,
  verbose = FALSE,
  n_purity = 100
)
}
\arguments{
\item{X}{a n x p matrix of covariates}

\item{betahat}{an n vector of observations}

\item{se}{a n vector of standard errors}

\item{Z}{fixed effect covaraites (including intercept). If null just a n x 1 matrix of ones}

\item{L}{number of single effects}

\item{f1_params}{list of parameters to initialize f1, must be a subset of those accepted by the constructor for f1}

\item{estimate_f1}{boolean whether or not to estimate the parameters of f1}

\item{prior_mean}{the prior mean of each non-zero element of b. Either a scalar or vector of length L.}

\item{prior_variance}{the prior variance of each non-zero element of b. Either a scalar or vector of length L. If \code{estimate_prior_variance=TRUE} the value provides an initial estimate of prior variances}

\item{prior_weights}{prior probability of selecting each column of X, vector of length p summing to one, or an L x p matrix}

\item{s_init}{a logistic susie object to initialize with, NOTE if non-null, we ignore \code{prior_mean}, \code{prior_variance}, and \code{prior_weights}}

\item{returns}{a fit Binomial SuSiE model, which is compatable with summary functions from \code{susieR} package}
}
\description{
Point Normal SuSiE
Fits Point Normal SuSiE, observations are drawn from a two component mixture
(point-mass at zero, normal).
The prior log odds of being non-zero are a linear function of covariates X
and the effects have a SuSiE prior
}

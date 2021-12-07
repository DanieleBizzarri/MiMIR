\name{covid_betas}
\alias{covid_betas}
\docType{data}
\title{
Parameters to calculate the susceptibility to CODIV-19 score
}
\description{
Abbreviatinos of the metabolites, metabolites and coefficients of the score
}
\usage{data("covid_betas")}
\format{
  A data frame with 25 observations on the following 3 variables.
  \describe{
    \item{\code{Abbreviation}}{a character vector}
    \item{\code{Metabolite}}{a character vector}
    \item{\code{Beta_value}}{a numeric vector}
  }
}
\source{
https://elifesciences.org/articles/63033
}
\examples{
data(covid_betas)
## maybe str(covid_betas) ; plot(covid_betas) ...
}
\keyword{datasets}

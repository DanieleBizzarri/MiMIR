\name{mort_betas}
\alias{mort_betas}
\docType{data}
\title{
Parameters to calculate the mortality score by Deelen et al.
}
\description{
Abbreviatinos of the metabolites, metabolites and coefficients of the score
}
\usage{data("mort_betas")}
\format{
  A data frame with 14 observations on the following 3 variables.
  \describe{
    \item{\code{Abbreviation}}{a character vector}
    \item{\code{Metabolite}}{a character vector}
    \item{\code{Beta_value}}{a numeric vector}
  }
}

\source{
https://www.nature.com/articles/s41467-019-11311-9
}
\examples{
data(mort_betas)
## maybe str(mort_betas) ; plot(mort_betas) ...
}
\keyword{datasets}

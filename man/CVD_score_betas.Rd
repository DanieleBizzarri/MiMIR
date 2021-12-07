\name{CVD_score_betas}
\alias{CVD_score_betas}
\docType{data}
\title{
Parameters to calculate the Cardiovascular disease score
}
\description{
Abbreviatinos of the metabolites, metabolites and coefficients of the score
}
\usage{data("CVD_score_betas")}
\format{
  A data frame with 12 observations on the following 3 variables.
  \describe{
    \item{\code{Abbreviation}}{a character vector}
    \item{\code{Metabolite}}{a character vector}
    \item{\code{Beta_value}}{a numeric vector}
  }
}
\source{
https://pubmed.ncbi.nlm.nih.gov/25573147/
}
\examples{
data(CVD_score_betas)
## maybe str(CVD_score_betas) ; plot(CVD_score_betas) ...
}
\keyword{datasets}

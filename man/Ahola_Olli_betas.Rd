\name{Ahola_Olli_betas}
\alias{Ahola_Olli_betas}
\docType{data}
\title{
Parameters to calculate the Ahola Olli type 2 diabetes score
}
\description{
Abbreviatinos of the metabolites, metabolites and coefficients of the score
}
\usage{data("Ahola_Olli_betas")}
\format{
  A data frame with 7 observations on the following 3 variables.
  \describe{
    \item{\code{Abbreviation}}{a character vector}
    \item{\code{Metabolite}}{a character vector}
    \item{\code{Beta_value}}{a numeric vector}
  }
}
\source{
https://pubmed.ncbi.nlm.nih.gov/31584131/
}
\examples{
data(Ahola_Olli_betas)
## maybe str(Ahola_Olli_betas) ; plot(Ahola_Olli_betas) ...
}
\keyword{datasets}

\name{metabo_names_translator}
\alias{metabo_names_translator}
\docType{data}
\title{
Metabolites names translators
}
\description{
This data.frame is used to translate alternative names of Nightingale dataset to the metabolites names used in BBMRI-nl
}
\usage{data("metabo_names_translator")}
\format{
  A data frame with 228 observations on the following 9 variables.
  \describe{
    \item{\code{abbreviation}}{a character vector}
    \item{\code{machine_readable_name}}{a character vector}
    \item{\code{name}}{a character vector}
    \item{\code{description}}{a character vector}
    \item{\code{alternative_names}}{a list vector}
    \item{\code{group}}{a character vector}
    \item{\code{subgroup}}{a character vector}
    \item{\code{unit}}{a character vector}
    \item{\code{BBMRI_names}}{a character vector}
  }
}
\details{
This file contains: 
abbreviations of the metabolites names,
machine_readable_name are the same abbreviations without spaces or other non readable characters
name is the full name of the metabolites
description what the metabolic feature represents
alternative_names other names used by Nightingale Health to indicate the same metabolites
group what group of variables it is part of
subgroup of the variables
unit measure of the metabolites
BBMMRI_names which are the names used in BBMRI for those variables
}
\source{
This file was first created in the package ggforestplot and then modified to include the names used in BBMRI-nl.
}
\references{
https://nightingalehealth.github.io/ggforestplot/articles/ggforestplot.html
}
\examples{
data(metabo_names_translator)
}
\keyword{datasets}

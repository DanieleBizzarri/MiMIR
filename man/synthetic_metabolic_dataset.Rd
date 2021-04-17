\name{synthetic_metabolic_dataset}
\alias{synthetic_metabolic_dataset}
\docType{data}
\title{
Synthetic example Nightingal metabolomics dataset made for the package MetaboRiSc
}
\description{
Synthetic example dataset made with synthpop
}
\usage{data("synthetic_metabolic_dataset")}
\format{
  The format is:
 num [1:5000, 1:228] 1.14e-12 2.11e-10 3.18e-11 1.05e-10 4.38e-10 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:5000] "Dummy1" "Dummy2" "Dummy3" "Dummy4" ...
  ..$ : chr [1:228] "xxl_vldl_p" "xxl_vldl_l" "xxl_vldl_pl" "xxl_vldl_c" ...
}
\examples{
data(synthetic_metabolic_dataset)
## maybe str(synthetic_metabolic_dataset) ; plot(synthetic_metabolic_dataset) ...
}
\keyword{datasets}

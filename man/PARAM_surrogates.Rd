\name{PARAM_surrogates}
\alias{PARAM_surrogates}
\docType{data}
\title{
Parameters to calculate the surrogate variables
}
\description{
This list contains all the informations to calculate the surrogate clinical variables by Bizzarri et all.
}
\usage{data("PARAM_surrogates")}
\format{
  The format is:
List of 6
 $ MET         : chr [1:56] "ala" "gln" "his" "phe" ...
 $ Nmax_miss   : num 1
 $ Nmax_zero   : num 1
 $ mean        :'data.frame':	1 obs. of  56 variables:
  ..$ ala      : num 0.315
  ..$ gln      : num 0.454
  ..$ his      : num 0.0581
  ..$ phe      : num 0.0472
  ..$ tyr      : num 0.0604
  ..$ ile      : num 0.0519
  ..$ leu      : num 0.0622
  ..$ val      : num 0.156
  ..$ glc      : num 4.68
  ..$ lac      : num 1.29
  ..$ cit      : num 0.103
  ..$ ace      : num 0.0455
  ..$ acace    : num 0.0423
  ..$ crea     : num 0.0709
  ..$ alb      : num 0.0872
  ..$ gp       : num 1.34
  ..$ m_vldl_l : num 0.647
  ..$ s_vldl_l : num 0.68
  ..$ xs_vldl_l: num 0.554
  ..$ idl_l    : num 1.09
  ..$ l_ldl_l  : num 1.21
  ..$ m_ldl_l  : num 0.691
  ..$ s_ldl_l  : num 0.446
  ..$ m_hdl_l  : num 0.825
  ..$ s_hdl_l  : num 1.03
  ..$ idl_c    : num 0.695
  ..$ serum_c  : num 4.47
  ..$ vldl_c   : num 0.838
  ..$ ldl_c    : num 1.55
  ..$ hdl_c    : num 1.38
  ..$ hdl2_c   : num 0.907
  ..$ hdl3_c   : num 0.476
  ..$ vldl_d   : num 36.8
  ..$ ldl_d    : num 23.6
  ..$ hdl_d    : num 9.96
  ..$ serum_tg : num 1.38
  ..$ totpg    : num 1.92
  ..$ pc       : num 1.99
  ..$ sm       : num 0.459
  ..$ totcho   : num 2.33
  ..$ apoa1    : num 1.58
  ..$ apob     : num 0.963
  ..$ totfa    : num 11.6
  ..$ dha      : num 0.141
  ..$ la       : num 3.05
  ..$ faw3     : num 0.416
  ..$ faw6     : num 3.83
  ..$ pufa     : num 4.25
  ..$ mufa     : num 2.97
  ..$ sfa      : num 4.36
  ..$ faw3_fa  : num 3.61
  ..$ faw6_fa  : num 33.4
  ..$ pufa_fa  : num 37
  ..$ mufa_fa  : num 25.4
  ..$ sfa_fa   : num 37.7
  ..$ unsatdeg : num 1.21
 $ sd          :'data.frame':	1 obs. of  56 variables:
  ..$ ala      : num 0.0794
  ..$ gln      : num 0.0761
  ..$ his      : num 0.014
  ..$ phe      : num 0.011
  ..$ tyr      : num 0.015
  ..$ ile      : num 0.0193
  ..$ leu      : num 0.0163
  ..$ val      : num 0.037
  ..$ glc      : num 1.56
  ..$ lac      : num 1.06
  ..$ cit      : num 0.0324
  ..$ ace      : num 0.0203
  ..$ acace    : num 0.0324
  ..$ crea     : num 0.0181
  ..$ alb      : num 0.00688
  ..$ gp       : num 0.205
  ..$ m_vldl_l : num 0.37
  ..$ s_vldl_l : num 0.236
  ..$ xs_vldl_l: num 0.131
  ..$ idl_l    : num 0.282
  ..$ l_ldl_l  : num 0.353
  ..$ m_ldl_l  : num 0.216
  ..$ s_ldl_l  : num 0.131
  ..$ m_hdl_l  : num 0.172
  ..$ s_hdl_l  : num 0.123
  ..$ idl_c    : num 0.195
  ..$ serum_c  : num 0.986
  ..$ vldl_c   : num 0.289
  ..$ ldl_c    : num 0.524
  ..$ hdl_c    : num 0.333
  ..$ hdl2_c   : num 0.31
  ..$ hdl3_c   : num 0.0625
  ..$ vldl_d   : num 1.35
  ..$ ldl_d    : num 0.134
  ..$ hdl_d    : num 0.259
  ..$ serum_tg : num 0.661
  ..$ totpg    : num 0.36
  ..$ pc       : num 0.366
  ..$ sm       : num 0.087
  ..$ totcho   : num 0.381
  ..$ apoa1    : num 0.202
  ..$ apob     : num 0.22
  ..$ totfa    : num 2.44
  ..$ dha      : num 0.0548
  ..$ la       : num 0.684
  ..$ faw3     : num 0.136
  ..$ faw6     : num 0.756
  ..$ pufa     : num 0.831
  ..$ mufa     : num 0.899
  ..$ sfa      : num 0.951
  ..$ faw3_fa  : num 0.948
  ..$ faw6_fa  : num 3.53
  ..$ pufa_fa  : num 3.63
  ..$ mufa_fa  : num 3.67
  ..$ sfa_fa   : num 1.89
  ..$ unsatdeg : num 0.0722
 $ models_betas: num [1:19, 1:57] -0.103 -2.397 -1.33 -0.165 -0.546 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:19] "s_sex" "s_diabetes" "s_lipidmed" "s_blood_pressure_lowering_med" ...
  .. ..$ : chr [1:57] "(Intercept)" "ala" "gln" "his" ...
}
\details{
 MET=Metabolites that are necessary for the logistic regression. Nmax_miss and Nmax_zero= the max number of missing and zeros raccommended for each sample. mean and sd= mean and standard deviation of each metabolites in BBMRI-NL. model_betas= a matrix with coefficients for the pre-trained models
}
\source{
reference to a publication or URL from which the data were obtained
}
\references{
possibly secondary sources and usages
}
\examples{
data(PARAM_surrogates)
## maybe str(PARAM_surrogates) ; plot(PARAM_surrogates) ...
}
\keyword{datasets}

#' acc_LOBOV
#' 
#' Accuracy of the Leave One Biobank Out Validation of the surrogate metabolic-modesl performed in BBMRI-nl
#'
#' @docType data
#'
#' @usage data("acc_LOBOV")
#'
#' @examples
#' data("acc_LOBOV")
#' 
#' @details 
#' Dataframe containing the accuracy obtained during the Leave One Biobank Out Validation of the surrogate metabolic-modesl  in BBMRI-nl.
#' 
#' @references
#' The method is described in:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' 
#' 
"acc_LOBOV"

#' T2D-score Betas
#' 
#' The coefficients used to compute the T2Diabetes score by Ahola Olli.
#'
#' @docType data
#'
#' @usage data("Ahola_Olli_betas")
#'
#' @examples
#' data("Ahola_Olli_betas")
#' 
#' @details 
#' Dataframe containing the abbreviation of the metabolites, the metabolites names and finally the Coefficients to compute the T2Diabetes score
#' 
#' @references
#' Ahola-Olli,A.V. et al. (2019) Circulating metabolites and the risk of type 2 diabetes: a prospective study of 11,896 young adults from four Finnish cohorts. Diabetologia, 62, 2298-2309, <doi:10.1007/s00125-019-05001-w>
#' 
"Ahola_Olli_betas"


#' COVID-score betas
#' 
#' The coefficients used to compute the COVID score by Nightingale Health UK Biobank Initiative et al.
#'
#' @docType data
#'
#' @usage data("covid_betas")
#'
#' @examples
#' data("covid_betas")
#' 
#' @details 
#' Dataframe containing the abbreviation of the metabolites, the metabolites names and finally the Coefficients to compute the COVID score
#' 
#' @references
#' Nightingale Health UK Biobank Initiative et al. (2021) Metabolic biomarker profiling for identification of susceptibility to severe pneumonia and COVID-19 in the general population. eLife, 10, e63033, <doi:10.7554/eLife.63033>
#' 
"covid_betas"

#' CVD-score betas
#' 
#' The coefficients used to compute the CVD score by Wurtz et al.
#'
#' @docType data
#'
#' @usage data("CVD_score_betas")
#'
#' @examples
#' data("CVD_score_betas")
#' 
#' @details 
#' Dataframe containing the abbreviation of the metabolites, the metabolites names and finally the Coefficients to compute the COVID score
#' 
#' @references
#' Wurtz,P. et al. (2015) Metabolite profiling and cardiovascular event risk: a prospective study of 3 population-based cohorts. Circulation, 131, 774-785, <doi:10.1161/CIRCULATIONAHA.114.013116>
#' 
"CVD_score_betas"


#' Mortality score betas
#' 
#' The coefficients used to compute the mortality score by Deelen et al.
#'
#' @docType data
#'
#' @usage data("mort_betas")
#'
#' @examples
#' data("mort_betas")
#' 
#' @details 
#' Dataframe containing the abbreviation of the metabolites, the metabolites names and finally the Coefficients to compute the mortality score
#' 
#' @references
#' Deelen,J. et al. (2019) A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals. Nature Communications, 10, 1-8, <doi:10.1038/s41467-019-11311-9>
#' 
"mort_betas"

#' PARAMETERS MetaboAge
#' 
#' The coefficients used to compute the MetaboAge by van den Akker et al.
#'
#' @docType data
#'
#' @usage data("PARAM_metaboAge")
#'
#' @examples
#' data("PARAM_metaboAge")
#' 
#' @details 
#' List containing all the information to pre-process and compute the MetaboAge.
#' 
#' @references
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541-547, <doi:10.1161/CIRCGEN.119.002610>
#' 
"PARAM_metaboAge"

#' PARAMETERS surrogates
#' 
#' The coefficients used to compute the metabolomics-based surrogate clinical variables by Bizzarri et al.
#'
#' @docType data
#'
#' @usage data("PARAM_surrogates")
#'
#' @examples
#' data("PARAM_surrogates")
#' 
#' @details 
#' List containing all the information to pre-process and compute the surrogate clinical variables.
#' 
#' @references
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' 
"PARAM_surrogates"

#' BBMRI_hist_scaled
#' 
#' Z-scaled distributions of the Nightingale Health metabolic features in BBMRI-nl
#'
#' @docType data
#'
#' @usage data("BBMRI_hist_scaled")
#'
#' @examples
#' data("BBMRI_hist_scaled")
#' 
#' @details 
#' List containing the histograms of the scaled metabolomics-features in BBMRI-nl.
#' 
#' 
"BBMRI_hist_scaled"

#' BBMRI_hist
#' 
#' Distributions of the Nightingale Health metabolic features in BBMRI-nl
#'
#' @docType data
#'
#' @usage data("BBMRI_hist")
#'
#' @examples
#' data("BBMRI_hist")
#' 
#' @details 
#' List containing the histograms of the metabolomics-features in BBMRI-nl.
#' 
"BBMRI_hist"

#' c21
#' 
#' Colors attributed to each metabolomics-based model in MiMIR
#'
#' @docType data
#'
#' @usage data("c21")
#'
#' @examples
#' data("c21")
#' 
"c21"


#' metabolomics feature nomenclatures
#' 
#' Translator of the names of the metabolomics-features to the ones used in BBMRI-nl
#'
#' @docType data
#'
#' @usage data("metabo_names_translator")
#'
#' @examples
#' data("metabo_names_translator")
#' 
#' @references 
#' This is a list originally created for the package ggforestplot and modified ad-hoc for our package 
#' (https://nightingalehealth.github.io/ggforestplot/articles/index.html).
#' 
"metabo_names_translator"

#' metabolomics feature subsets
#' 
#' List containing all the subset of the metabolomics-based features used for our models
#'
#' @docType data
#'
#' @usage data("metabolites_subsets")
#'
#' @examples
#' data("metabolites_subsets")
#' 
#' @references 
#' The selection of metabolic features available is the one selected by the papers:
#' Deelen,J. et al. (2019) A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals. Nature Communications, 10, 1-8, <doi:10.1038/s41467-019-11311-9>
#' Ahola-Olli,A.V. et al. (2019) Circulating metabolites and the risk of type 2 diabetes: a prospective study of 11,896 young adults from four Finnish cohorts. Diabetologia, 62, 2298-2309, <doi:10.1007/s00125-019-05001-w>
#' Wurtz,P. et al. (2015) Metabolite profiling and cardiovascular event risk: a prospective study of 3 population-based cohorts. Circulation, 131, 774-785, <doi:10.1161/CIRCULATIONAHA.114.013116>
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541-547, <doi:10.1161/CIRCGEN.119.002610>
#' 
"metabolites_subsets"

#' phenotypic features names
#' 
#' List containing all the subsets of phenotypics variables used in the app
#'
#' @docType data
#'
#' @usage data("phenotypes_names")
#'
#' @examples
#' data("phenotypes_names")
#'
"phenotypes_names"

#' synthetic metabolomics dataset
#' 
#' Data.frame containing a synthetic dataset of the Nightingale Metabolomics dataset created with the package synthpop from the LLS_PAROFF dataset.
#'
#' @docType data
#'
#' @usage data("synthetic_metabolic_dataset")
#'
#' @examples
#' data("synthetic_metabolic_dataset")
#'
#' @references 
#' M. Schoenmaker et al., 'Evidence of genetic enrichment for exceptional survival using a family approach: the Leiden Longevity Study', Eur. J. Hum. Genet., vol. 14, no. 1, Art. no. 1, Jan. 2006, <doi:10.1038/sj.ejhg.5201508>
#' B. Nowok, G. M. Raab, and C. Dibben, 'synthpop: Bespoke Creation of Synthetic Data in R', J. Stat. Softw., vol. 74, no. 1, Art. no. 1, Oct. 2016, <doi:10.18637/jss.v074.i11>
#'
"synthetic_metabolic_dataset"

#' synthetic metabolomics dataset
#' 
#' Data.frame containing a synthetic dataset of phenotypic dataset created with the package synthpop from the LLS_PAROFF dataset.
#'
#' @docType data
#'
#' @usage data("synthetic_metabolic_dataset")
#'
#' @examples
#' data("synthetic_metabolic_dataset")
#'
#' @references
#' M. Schoenmaker et al., 'Evidence of genetic enrichment for exceptional survival using a family approach: the Leiden Longevity Study', Eur. J. Hum. Genet., vol. 14, no. 1, Art. no. 1, Jan. 2006, <doi:10.1038/sj.ejhg.5201508>
#' B. Nowok, G. M. Raab, and C. Dibben, 'synthpop: Bespoke Creation of Synthetic Data in R', J. Stat. Softw., vol. 74, no. 1, Art. no. 1, Oct. 2016, <doi:10.18637/jss.v074.i11>
#'
"synthetic_phenotypic_dataset"



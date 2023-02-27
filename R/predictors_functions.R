#####################
### Libraries #######
#####################
#Import packages
#' @import foreach
#' @importFrom purrr map_chr
#' @importFrom purrr map_lgl
#' @importFrom caret createDataPartition
#' @importFrom stats model.matrix
#' @importFrom stats na.omit
#' @importFrom matrixStats colMedians
#' @importFrom matrixStats colSds
#' @import shiny
#' @import survival
#' @import survminer
#' @import dplyr
#' @importFrom stats as.formula  
#' @importFrom stats coefficients
#' @importFrom stats p.adjust
#' @importFrom graphics par
NULL

###################################
### DEFINITIONS / CONSTANTS #######
###################################
#Global variables
utils::globalVariables(c("bin_phenotypes", "i", "y", "mortScore", 
                         "surro","AUC","outcome", "Uploaded", "biobank","MiMIR::BBMRI_translator",
                         "BBMRI_hist_scaled","BBMRI_hist", "ord", "pval.adj",
                         "met_name", "prepped_dat", "CVD_score", "MiMIR::metabo_names_translator",
                         "PARAM_metaboAge","PARAM_surrogates",
                         "mort_betas", "Ahola_Olli_betas", "CVD_score_betas", "MiMIR::metabolites_subsets", "covid_betas",
                         "pheno_names", "out_list", "bin_names", "bin_surro", "MiMIR::c21","dropdownMenuOutput"))

######################
## Upload functions ##
######################
#' find_BBMRI_names
#' 
#' Function to translate Nightingale metabolomics alternative metabolite names to the ones used in BBMRI-nl
#'
#' @param names vector of strings with the metabolic features names to be translated
#' @return data.frame with the uploaded metabolites names on the first column and the BBMRI names on the second column.
#' 
#' @examples
#' library(MiMIR)
#' library(purrr)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' #Find the metabolites names used in BBMRI-nl
#' nam<-find_BBMRI_names(colnames(metabolic_measures)) 
#' 
#' @references 
#' This is a function originally created for the package ggforestplot and modified ad hoc for our package 
#' (https://nightingalehealth.github.io/ggforestplot/articles/index.html).
#' 
#' @export
#' 
find_BBMRI_names<-function(names){
  names<-tolower(names)
  new_names <- names %>% purrr::map_chr(function(id) {
    # Look through the alternative_ids
    hits <-
      purrr::map_lgl(
        MiMIR::metabo_names_translator$alternative_names,
        ~ id %in% .
      )
    
    # If one unambiguous hit, return it.
    if (sum(hits) == 1L) {
      return(MiMIR::metabo_names_translator$BBMRI_names[hits])
      # If not found, give a warning and pass through the input.
    } else {
      warning("Biomarker not found: ", id, call. = FALSE)
      return(id)
    } 
  })
  n<-data.frame(uploaded=names,BBMRI_names=new_names)
  return(n)
}

########################################
## Metabolomics-based score functions ##
########################################
#' prep_met_for_scores
#' 
#' Helper function to pre-process the Nightingale Health metabolomics data-set before applying the mortality, Type-2-diabetes and CVD scores.
#'
#' @param dat numeric data-frame with Nightingale-metabolomics
#' @param featID vector of strings with the names of metabolic features included in the score selected
#' @param plusone logical to determine if a value of 1.0 should be added to all metabolic features (TRUE) or only to the ones featuring zeros before log-transforming (FALSE)
#' @param quiet logical to suppress the messages in the console
#' @return The Nightingale-metabolomics data-frame after pre-processing (checked for zeros, zscale and log-transformed) according to what has been done by the authors of the original papers.
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the Nightingale metabolomics dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' #Prepare the metabolic features fo the mortality score
#' prepped_met <- prep_met_for_scores(metabolic_measures,featID=MiMIR::mort_betas$Abbreviation)
#' 
#' @references 
#' This function is constructed to be able to follow the pre-processing steps described in:
#' Deelen,J. et al. (2019) A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals. Nature Communications, 10, 1-8, <doi:10.1038/s41467-019-11311-9>.
#'  
#' Ahola-Olli,A.V. et al. (2019) Circulating metabolites and the risk of type 2 diabetes: a prospective study of 11,896 young adults from four Finnish cohorts. Diabetologia, 62, 2298-2309, <doi:10.1007/s00125-019-05001-w>
#' 
#' Wurtz,P. et al. (2015) Metabolite profiling and cardiovascular event risk: a prospective study of 3 population-based cohorts. Circulation, 131, 774-785, <doi:10.1161/CIRCULATIONAHA.114.013116>
#' 
#' @seealso 
#' comp.mort_score, mort_betas, comp.T2D_Ahola_Olli, comp.CVD_score
#' 
#' @export
#' 
prep_met_for_scores <- function(dat, featID, plusone=FALSE, quiet=FALSE){
  if(!quiet){
    cat("|| Preparing data ... \n")
  }
  ## 1. Check for zeroes:
  if(plusone){
    if(!quiet){
    cat("| Adding 1.0 to all the metabolites\n")
    }
    dat <- dat[,featID] + 1
  }else{
    to_fix <- names(which(colSums(dat[,featID]==0,na.rm=TRUE)>0))
    if(length(to_fix)>0){
      if(!quiet){
        cat(paste0("| Adding 1.0 to metabolites featuring zero's: '",paste(to_fix,collapse="', '"),"'\n"))
      }
      dat[,to_fix] <- dat[,to_fix] + 1
    } else {
      if(!quiet){
        cat(paste0("| No metabolites found featuring zero's \n"))
      }
    }
  }
  ## 2. Scale:
  dat[,featID] <- scale(log(dat[,featID]),center=TRUE,scale=TRUE)
  if(!quiet){
    cat("| Perform log transform & scaling to zero mean and unity sd .. ")
  }
  
  if(!quiet){
    cat("Done!\n")
  }
  return(dat)
}

#' comp.mort_score
#' 
#' Function to compute the mortality score made by Deelen et al. on Nightingale metabolomics data-set.
#'
#' @param dat numeric data-frame with Nightingale-metabolomics
#' @param betas data.frame containing the coefficients used for the regression of the mortality score
#' @param quiet logical to suppress the messages in the console
#' @return data-frame containing the value of the mortality score on the uploaded data-set
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' #Prepare the metabolic features fo the mortality score
#' mortScore<-comp.mort_score(metabolic_measures,quiet=TRUE)
#' 
#' @details 
#' This multivariate model predicts all-cause mortality at 5 or 10 years better than clinical variables normally associated with mortality. 
#' It is constituted of 14 metabolic features quantified by Nightingale Health. 
#' It was originally trained using a stepwise Cox regression analysis in a meta-analysis on 12 cohorts composed by 44,168 individuals.
#' 
#' @references 
#' This function is constructed to be able to apply the mortality score as described in:
#' Deelen,J. et al. (2019) A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals. Nature Communications, 10, 1-8, <doi:10.1038/s41467-019-11311-9>
#' 
#' @seealso 
#' prep_met_for_scores, mort_betas, comp.T2D_Ahola_Olli, comp.CVD_score
#'  
#' @export
#' 
comp.mort_score <- function(dat, betas=mort_betas, quiet=FALSE){
  ## 1. Prepare data:
  if(!quiet){
    cat("=== Computing mortality score === \n")
  }
  prepped_dat <- prep_met_for_scores(dat,featID=betas$Abbreviation, plusone=F, quiet=quiet)
  ## 2. Compute:
  if(!quiet){
    cat("| Computing score .. \n")
  }
  mortScore <- as.data.frame(as.matrix(prepped_dat[,betas$Abbreviation]) %*% betas$Beta_value)
  colnames(mortScore)<-"mortScore"
  rownames(mortScore)<-rownames(dat)
  
  if(!quiet){
    cat("Done!\n\n")
  }
  return(mortScore)
}

#' comp.T2D_Ahola_Olli
#' 
#' Function to compute the T2D score made by Ahola Olli et al. on Nightingale metabolomics data-set.
#'
#' @param met numeric data-frame with Nightingale-metabolomics
#' @param phen data-frame containing phenotypic information of the samples (in particular: sex, age, BMI and the clinically measured glucose)
#' @param betas The betas of the linear regression composing the T2D-score
#' @param quiet logical to suppress the messages in the console
#' @return data-frame containing the value of the T2D-score on the uploaded data-set
#' 
#' @details 
#' This metabolomics-based score is associated with incident Type 2 Diabetes, made by Ahola-Olli et al. 
#' It is constructed using phe, l_vldl_ce_percentage and l_hdl_fc quantified by Nightingale Health,
#' and some phenotypic information: sex, age, BMI, fasting glucose.
#' It was trained using a stepwise logistic regression on 3 cohorts.
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the dataset
#' met <- synthetic_metabolic_dataset
#' phen<-synthetic_phenotypic_dataset
#' #Prepare the metabolic features fo the mortality score
#' T2Dscore<-comp.T2D_Ahola_Olli(met= met, phen=phen,betas=MiMIR::Ahola_Olli_betas, quiet=TRUE)
#' 
#' @references 
#' This function is constructed to be able to apply the T2D-score as described in:
#' Ahola-Olli,A.V. et al. (2019) Circulating metabolites and the risk of type 2 diabetes: a prospective study of 11,896 young adults from four Finnish cohorts. Diabetologia, 62, 2298-2309, <doi:10.1007/s00125-019-05001-w>
#' 
#' @seealso 
#' prep_met_for_scores, Ahola_Olli_betas, comp.mort_score, comp.CVD_score
#' @export
#' 
comp.T2D_Ahola_Olli<- function(met, phen, betas, quiet=FALSE){
  ## 1. Prepare data:
  if(!quiet){
    cat("=== Computing Ahola Olli T2D score === \n")
  }
  prepped_met <- prep_met_for_scores(met,featID=betas$Abbreviation[1:3],plusone=T,quiet=quiet)
  ## 2. Compute:
  if(!quiet){
    cat("| Computing score .. \n")
  }
  prepped_data<-cbind(met,phen)
  if(length(which(colnames(prepped_data) %in% betas$Abbreviation))==7){
    if(any(colSums(is.na(prepped_data[,betas$Abbreviation])) == dim(prepped_data)[1])){
      miss_phen<-names(which(colSums(is.na(prepped_data[,betas$Abbreviation])) == dim(prepped_data)[1]))
      return(paste("It was not possible to calculate the score because all the values of", paste(miss_phen, collapse=", "), "are missing"))
    }else{
      T2DScore <- as.data.frame(as.matrix(prepped_data[,betas$Abbreviation]) %*% betas$Beta_value)
      colnames(T2DScore)<-"T2Dscore"
      rownames(T2DScore)<-rownames(met)
      
      if(!quiet){
        cat("Done!\n\n")
      }
      return(T2DScore)
    }
  }else{
    i<-colnames(prepped_data)[which(colnames(prepped_data) %in% betas$Abbreviation)] 
    miss_phen<-setdiff(betas$Abbreviation, i)
    return(paste("It was not possible to calculate the score because of the missing variable(s):", paste(miss_phen, collapse=", ")))
  }
}


#' comp.CVD_score
#' 
#' Function to compute CVD-score made by Peter Wurtz et al. made by Deelen et al. on Nightingale metabolomics data-set.
#'
#' @param met numeric data-frame with Nightingale-metabolomics
#' @param phen data-frame containing phenotypic information of the samples (specifically: sex, systolic_blood_pressure, current_smoking, diabetes, blood_pressure_lowering_med, lipidmed, totchol, and hdlchol)
#' @param betas The betas of the linear regression composing the CVD-score
#' @param quiet logical to suppress the messages in the console
#' @return data-frame containing the value of the CVD-score on the uploaded data-set
#' 
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the dataset
#' met <- synthetic_metabolic_dataset
#' phen<-synthetic_phenotypic_dataset
#' #Prepare the metabolic features fo the mortality score
#' CVDscore<-comp.CVD_score(met= met, phen=phen, betas=MiMIR::CVD_score_betas, quiet=TRUE)
#' 
#' 
#' @references 
#' This function is constructed to be able to apply the CVD-score as described in:
#' Wurtz,P. et al. (2015) Metabolite profiling and cardiovascular event risk: a prospective study of 3 population-based cohorts. Circulation, 131, 774-785, <doi:10.1161/CIRCULATIONAHA.114.013116>
#' 
#' @seealso 
#' prep_met_for_scores, CVD_score_betas, comp.T2D_Ahola_Olli, comp.mort_score
#'  
#' @export
#' 
comp.CVD_score <- function(met, phen, betas, quiet=FALSE){
  ## 1. Prepare data:
  if(!quiet){
    cat("=== Computing CVD score === \n")
  }
  prepped_met <- prep_met_for_scores(met,featID=betas$Abbreviation[1:4],plusone=FALSE,quiet=quiet)
  ## 2. Compute:
  if(!quiet){
    cat("| Computing score .. \n")
  }
  prepped_data<-cbind(met,phen)
  if(length(which(colnames(prepped_data)%in% betas$Abbreviation))==12){
    if(any(colSums(is.na(prepped_data[,betas$Abbreviation])) == dim(prepped_data)[1])){
      miss_phen<-names(which(colSums(is.na(prepped_data[,betas$Abbreviation])) == dim(prepped_data)[1]))
      return(paste("It was not possible to calculate the score because all the values of", paste(miss_phen, collapse=", "), "are missing"))
    }else{
      CVD_score <- as.data.frame(as.matrix(prepped_data[,betas$Abbreviation]) %*% betas$Beta_value)
      colnames(CVD_score)<-"CVD_score"
      rownames(CVD_score)<-rownames(met)
      
      if(!quiet){
        cat("Done!\n\n")
      }
      return(CVD_score)
    }
  } else{
    i<-colnames(prepped_data)[which(colnames(prepped_data) %in% betas$Abbreviation)] 
    miss_phen<-setdiff(betas$Abbreviation, i)
    return(paste("It was not possible to calculate the score because of the missing variable(s):", paste(miss_phen, collapse=", ")))
  }
  
}


#' prep_data_COVID_score
#' 
#' Helper function to pre-process the Nightingale Health metabolomics data-set before applying the COVID score.
#'
#' @param dat numeric data-frame with Nightingale-metabolomics
#' @param featID vector of strings with the names of metabolic features included in the COVID-score
#' @param quiet logical to suppress the messages in the console
#' @return The Nightingale-metabolomics data-frame after pre-processing (checked for zeros, z-scaled and log-transformed) according to what has been done by the authors of the original papers.
#' 
#' @examples
#' require(MiMIR)
#' require(matrixStats)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' #Prepare the metabolic features fo the mortality score
#' prepped_met <- prep_data_COVID_score(dat=metabolic_measures)
#' 
#' 
#' @references 
#' This function is constructed to be able to follow the pre-processing steps described in:
#' Nightingale Health UK Biobank Initiative et al. (2021) Metabolic biomarker profiling for identification of susceptibility to severe pneumonia and COVID-19 in the general population. eLife, 10, e63033, <doi:10.7554/eLife.63033>
#' 
#' @seealso 
#' prep_met_for_scores, covid_betas, comp_covid_score
#' 
#' @export
#'
prep_data_COVID_score <- function(dat, featID=c("gp","dha","crea","mufa", "apob_apoa1","tyr","ile","sfa_fa","glc","lac","faw6_faw3",
                                               "phe", "serum_c", "faw6_fa","ala","pufa","glycine","his","pufa_fa","val","leu",
                                               "alb","faw3","ldl_c","serum_tg"), quiet=FALSE){
  if(!quiet){
    cat("|| Preparing data ... \n")
  }
  
  if(length(which(colnames(dat)=="faw6_faw3"))==0){
    dat$faw6_faw3<-dat$faw6/dat$faw3
  }
  MEANS=matrixStats::colMedians(as.matrix(dat[,featID]),na.rm=T)
  names(MEANS)<-featID
  SDS=matrixStats::colSds(as.matrix(dat[,featID]),na.rm=T)
  names(SDS)<-featID
  # 1. Subset samples on SD:
  dat <- subset_samples_sd_surrogates(x=as.matrix(dat[,featID]), MEAN=MEANS, SD=SDS,N=4, quiet=quiet)
  
  ## 2. Scale:
  dat[,featID] <- scale(dat[,featID],center=TRUE,scale=TRUE)
  if(!quiet){
    cat("| z-scaling to zero mean and unity sd .. ")
  }
  
  if(!quiet){
    cat("Done!\n")
  }
  return(dat)
}


#' comp_covid_score
#' 
#' Function to compute the COVID severity score made by Nightingale Health UK Biobank Initiative et al. on Nightingale metabolomics data-set.
#'
#' @param dat numeric data-frame with Nightingale-metabolomics
#' @param betas data.frame containing the coefficients used for the regression of the COVID-score
#' @param quiet logical to suppress the messages in the console
#' @return data-frame containing the value of the COVID-score on the uploaded data-set
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' 
#' #Compute the mortality score
#' mortScore<-comp_covid_score(dat=metabolic_measures, quiet=TRUE)
#' 
#' @details 
#' Multivariate model predicting the risk of severe COVID-19 infection. 
#' It is based on 37 metabolic features and trained using LASSO regression on 52,573 samples from the UK-biobanks.
#' 
#' @references 
#' This function is constructed to be able to apply the COVID-score as described in:
#' Nightingale Health UK Biobank Initiative et al. (2021) Metabolic biomarker profiling for identification of susceptibility to severe pneumonia and COVID-19 in the general population. eLife, 10, e63033, <doi:10.7554/eLife.63033>
#' 
#' @seealso 
#' prep_data_COVID_score, covid_betas, comp.mort_score
#'  
#' @export
#' 
comp_covid_score <- function(dat, betas=MiMIR::covid_betas ,quiet=FALSE){
  ## 1. Prepare data:
  if(!quiet){
    cat("=== Computing mortality score === \n")
  }
  prepped_data <- prep_data_COVID_score(dat,featID=betas$Abbreviation,quiet=quiet)
  ## 2. Compute:
  if(!quiet){
    cat("| Computing score .. \n")
  }
  
  
  if(length(which(colnames(prepped_data)%in% betas$Abbreviation))==25){
    if(any(colSums(is.na(prepped_data[,betas$Abbreviation])) == dim(prepped_data)[1])){
      miss_phen<-names(which(colSums(is.na(prepped_data[,betas$Abbreviation])) == dim(prepped_data)[1]))
      return(paste("It was not possible to calculate the score because all the values of", paste(miss_phen, collapse=", "), "are missing"))
    }else{
      covidScore<-as.data.frame(as.matrix(prepped_data[,betas$Abbreviation]) %*% betas$Beta_value)
      colnames(covidScore)<-"covidScore"
      rownames(covidScore)<-rownames(prepped_data)
      
      if(!quiet){
        cat("Done!\n\n")
      }
      return(covidScore)
    }
  } else{
    i<-colnames(prepped_data)[which(colnames(prepped_data) %in% betas$Abbreviation)] 
    miss_phen<-setdiff(betas$Abbreviation, i)
    return(paste("It was not possible to calculate the score because of the missing variable(s):", paste(miss_phen, collapse=", ")))
  }
  
  if(!quiet){
    cat("Done!\n\n")
  }
  return(covidScore)
}


#########################
## MetaboAge functions ##
#########################
#' QCprep
#' 
#' Helper function to pre-process the Nightingale Health metabolomics data-set before applying the MetaboAge score by van den Akker et al.
#'
#' @param mat numeric data-frame NH-metabolomics matrix.
#' @param PARAM_metaboAge list containing all the parameters to compute the metaboAge (metabolic features list,BBMRI-nl means and SDs of the metabolic features, and coefficients)
#' @param quiet logical to suppress the messages in the console
#' @param Nmax_zero numberic value indicating the maximum number of zeros allowed per sample (Number suggested=1)
#' @param Nmax_miss numberic value indicating the maximum number of missing values allowed per sample (Number suggested=1)
#' @return Nightingale-metabolomics data-frame after pre-processing (checked for zeros, missing values, samples>5SD from the BBMRI-mean, imputing the missing values and z-scaled)
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' 
#' #Pre-process the metabolic features
#' prepped_met<-QCprep(as.matrix(metabolic_measures[,metabolites_subsets$MET63]), PARAM_metaboAge)
#' 
#' @references 
#' This function is constructed to be able to follow the pre-processing steps described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541-547, <doi:10.1161/CIRCULATIONAHA.114.013116>
#' 
#' @seealso 
#' apply.fit
#' 
#' @export
#'
QCprep<-function(mat, PARAM_metaboAge, quiet=TRUE, Nmax_zero=1, Nmax_miss=1){
  # 0. Start:
  if(!quiet){
    cat(report.dim(mat,header="Start"))
  }
  # 1. Subset required metabolites:
  mat <- subset_metabolites_overlap(mat,metabos=PARAM_metaboAge$MET,quiet=quiet)
  # 2. Subset samples on missingness:
  mat <- subset_samples_miss(mat, Nmax=Nmax_miss, quiet=quiet)
  # 3. Subset samples on zeros:
  mat <- subset_samples_zero(mat, Nmax=Nmax_zero, quiet=quiet)
  # 4. Subset samples on SD:
  mat <- subset_samples_sd(as.matrix(mat), MEAN=PARAM_metaboAge$logMEAN, SD=PARAM_metaboAge$logSD, quiet=quiet)
  # 5. Perform scaling:
  mat <- apply.scale(mat, MEAN=PARAM_metaboAge$MEAN, SD=PARAM_metaboAge$SD)
  if(!quiet){
    cat("| Performing scaling ... ")
    cat(" DONE!\n")
  }
  # 6. Perform imputation:
  mat <- impute_miss(mat)
  if(!quiet){
    cat("| Imputation ... ")
    cat(" DONE!\n")
  }
  return(mat)
}

#' apply.fit
#' 
#' Function to compute the MetaboAge score made by van den Akker et al. on Nightingale metabolomics data-set.
#'
#' @param mat numeric data-frame with Nightingale-metabolomics
#' @param FIT The betas of the linear regression composing the MetaboAge by van den Akker et al.
#' @return data-frame containing the value of the MetaboAge by van den Akker et al.
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' #Pre-process the metabolic features
#' prepped_met<-QCprep(as.matrix(metabolic_measures[,metabolites_subsets$MET63]), PARAM_metaboAge)
#' #Apply the metaboAge
#' metaboAge<-apply.fit(prepped_met, FIT=PARAM_metaboAge$FIT_COEF)
#' 
#' @details 
#' Multivariate model indicating the biological age of an individual, based on 56 metabolic features. 
#' It was trained using a linear regression in BBMRI-nl, a Consortium of 28 cohorts comprising ~25,000 individuals.
#' 
#' @references 
#' This function is constructed to be able to apply the metaboAge as described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541-547, <doi:10.1161/CIRCULATIONAHA.114.013116>
#' 
#' @seealso 
#' QCprep, subset_metabolites_overlap, subset_samples_miss, subset_samples_zero, subset_samples_sd, impute_miss, apply.scale,report.dim
#'  
#' @keywords internal
#'  
#' @export
#' 
apply.fit<-function(mat,FIT){
  # Resort:
  BETA <- FIT[colnames(mat)]
  INTC <- FIT[1]
  # Predict age:
  AGE <- data.frame(MetaboAge=as.vector(mat %*% BETA) + as.vector(INTC),stringsAsFactors=FALSE)
  rownames(AGE)<-rownames(mat)
  return(AGE)
}
######################################
## Metabolomics-Surrogate functions ##
######################################
#' BMI_LDL_eGFR
#' 
#' #' Function created to calculate: 1) BMI using height and weight; 2) LDL cholesterol using HDL cholesterol, triglycerides, totchol;
#' 3) eGFR creatinine levels, sex and age.
#'
#' @param phenotypes data.frame containing height and weight, HDL cholesterol, triglycerides, totchol, sex and age
#' @param metabo_measures numeric data-frame with Nightingale metabolomics quantifications containing creatinine levels (crea)
#' @return phenotypes data.frame with the addition of BMI, LDL cholesterol and eGFR
#' @export
#'
#' @examples
#' library(MiMIR)
#' 
#' #load the dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' phenotypes <- synthetic_phenotypic_dataset
#' #Calculate BMI, LDL cholesterol and eGFR
#' phenotypes<-BMI_LDL_eGFR(phenotypes, metabolic_measures)
#' 
#' @references 
#' This function is constructed to calculate BMI, LDL cholesterol and eGFR as in the following papers:
#' 
#' BMI: Flint AJ, Rexrode KM, Hu FB, Glynn RJ, Caspard H, Manson JE et al. Body mass index, waist circumference, and risk of coronary heart disease: a prospective study among men and women. Obes Res Clin Pract 2010; 4: e171-e181, <doi:10.1016/j.orcp.2010.01.001>
#' 
#' LDL-cholesterol: Friedewald WT, Levy RI, Fredrickson DS. Estimation of the Concentration of Low-Density Lipoprotein Cholesterol in Plasma, Without Use of the Preparative Ultracentrifuge. Clin Chem 1972; 18: 499-502, <doi.org/10.1093/clinchem/18.6.499>
#' 
#' eGFR: Carrero Juan Jesus, Andersson Franko Mikael, Obergfell Achim, Gabrielsen Anders, Jernberg Tomas. hsCRP Level and the Risk of Death or Recurrent Cardiovascular Events in Patients With Myocardial Infarction: a Healthcare-Based Study. J Am Heart Assoc 2019; 8: e012638, <doi: 10.1161/JAHA.119.012638>
#' 
#' @export
#'
BMI_LDL_eGFR<-function(phenotypes, metabo_measures){
  #BMI calculation
  if(!c("BMI") %in% colnames(phenotypes)){
    if(all(c("weight","height") %in% colnames(phenotypes))){
      phenotypes<-cbind(phenotypes, data.frame(BMI=phenotypes$weight/(phenotypes$height/100)^2))
    }
  }
  
  #ldl-chol calculation
  if(!c("ldl_chol") %in% colnames(phenotypes)){
    if(all(c("hdlchol","triglycerides", "totchol") %in% colnames(phenotypes))){
    phenotypes<-cbind(phenotypes, data.frame(
      ldl_chol= (phenotypes$totchol-(phenotypes$hdlchol + (0.45 * phenotypes$triglycerides)))))
    phenotypes[which(phenotypes$triglycerides>=4.52),"ldl_chol"]<-NA
    }
  }
  
  #eGFR calculation
  if(!c("eGFR") %in% colnames(phenotypes)){
    if((all(c("sex","age") %in% colnames(phenotypes))) & 
       c("crea") %in% colnames(metabo_measures)){
      crea<-metabo_measures$crea/0.0884
      eGFR<-rep(NA, dim(phenotypes)[1])
      male_ind<-which(phenotypes$sex==1)
      m_indmin<-which(crea[male_ind]<=0.9)
      m_indplus<-which(crea[male_ind]>0.9)
      female_ind<-which(phenotypes$sex==0)
      f_indmin<-which(crea[female_ind]<=0.7)
      f_indplus<-which(crea[female_ind]>0.7)
      eGFR[male_ind[m_indmin]]<-141*((crea[male_ind[m_indmin]]/0.9)^(-0.411))*((0.993)^phenotypes[male_ind[m_indmin],"age"])
      eGFR[male_ind[m_indplus]]<-141*((crea[male_ind[m_indplus]]/0.9)^(-1.209))*((0.993)^phenotypes[male_ind[m_indplus],"age"])
      eGFR[female_ind[f_indmin]]<-144*((crea[female_ind[f_indmin]]/0.7)^(-0.329))*((0.993)^phenotypes[female_ind[f_indmin],"age"])
      eGFR[female_ind[f_indplus]]<-144*((crea[female_ind[f_indplus]]/0.7)^(-1.209))*((0.993)^phenotypes[female_ind[f_indplus],"age"])
      phenotypes<-cbind(phenotypes,eGFR)
      }
    }
  return(phenotypes)
}

#' binarize_all_pheno
#' 
#' Helper function created to binarize the phenotypes used to calculate the metabolomics based surrogate made by Bizzarri et al.
#'
#' @param data phenotypes data.frame containing some of the following variables (with the same namenclature):
#' "sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking",
#' "metabolic_syndrome", "alcohol_consumption", "age","BMI", "ln_hscrp","waist_circumference",
#' "weight","height", "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb" 
#' @return The phenotypic variables binarized following the thresholds in in the metabolomics surrogates made by by Bizzarri et al.
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the phenotypes dataset
#' phenotypes <- synthetic_phenotypic_dataset
#' #Calculate BMI, LDL cholesterol and eGFR
#' binarized_phenotypes<-binarize_all_pheno(phenotypes)
#' 
#' @details 
#' Bizzarri et al. built multivariate models,using 56 metabolic features quantified by Nightingale, to predict the 19 binary characteristics of an individual. 
#' The binary variables are: sex, diabetes status, metabolic syndrome status, lipid medication usage, blood pressure lowering medication,
#' current smoking, alcohol consumption, high age, middle age, low age, high hsCRP, high triglycerides, high ldl cholesterol,
#' high total cholesterol, low hdl cholesterol, low eGFR, low white blood cells, low hemoglobin levels.
#' 
#' @references 
#' This function was made to binarize the variables following the same rules indicated in the article:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' 
#' @seealso 
#' pheno_barplots
#'  
#' @export
#'
binarize_all_pheno<-function(data){
  available_pheno<-intersect(colnames(data),c("sex","diabetes","lipidmed","current_smoking", "blood_pressure_lowering_med", "alcohol_consumption", "metabolic_syndrome"))
  binarized_phenotypes<-cbind(data[,available_pheno],
                              matrix(data=NA,nrow=dim(data)[1],ncol=12))
  colnames(binarized_phenotypes)<-c(available_pheno,
                                    "high_age", "middle_age","low_age","high_hscrp", "high_triglycerides","low_hdlchol","high_ldl_chol","high_totchol",
                                    "low_eGFR","obesity", "low_wbc","low_hgb")
  #gender
  if(c("sex") %in% colnames(data)){
    if(!is.null(which(data$sex=="male"))){
      binarized_phenotypes[which(binarized_phenotypes$sex=="male"),"sex"]<-1
      binarized_phenotypes[which(binarized_phenotypes$sex=="female"),"sex"]<-0
    }
    binarized_phenotypes$sex<-as.factor(binarized_phenotypes$sex)
  }
  #diabetes
  if(c("diabetes") %in% colnames(data)){
    if(!is.null(which(data$diabetes=="TRUE"))){
      binarized_phenotypes[which(binarized_phenotypes$diabetes=="TRUE"),"diabetes"]<-1
      binarized_phenotypes[which(binarized_phenotypes$diabetes=="FALSE"),"diabetes"]<-0
    }
    binarized_phenotypes$diabetes<-as.factor(binarized_phenotypes$diabetes)
  }
  #lipidmed
  if(c("lipidmed") %in% colnames(data)){
    if(!is.null(which(data$lipidmed=="statins"))){
      binarized_phenotypes[which(binarized_phenotypes$lipidmed=="statins"),"lipidmed"]<-1
      binarized_phenotypes[!which(binarized_phenotypes$lipidmed=="statins"),"lipidmed"]<-0
      binarized_phenotypes[which(is.na(data$lipidmed)),"lipidmed"]<-NA
    }
    
    if(length(levels(factor(data$lipidmed)))>2){
      binarized_phenotypes[which(binarized_phenotypes$lipidmed==2),"lipidmed"]<-NA
    }
    binarized_phenotypes$lipidmed<-as.factor(binarized_phenotypes$lipidmed)
  }
  #current smoking
  if(c("current_smoking") %in% colnames(data)){
    if(!is.null(which(data$current_smoking=="TRUE"))){
      binarized_phenotypes[which(binarized_phenotypes$current_smoking=="TRUE"),"current_smoking"]<-1
      binarized_phenotypes[which(binarized_phenotypes$current_smoking=="FALSE"),"current_smoking"]<-0
    }
    binarized_phenotypes$current_smoking<-as.factor(binarized_phenotypes$current_smoking)
  }
  #blood_pressure_lowering_med
  if(c("blood_pressure_lowering_med") %in% colnames(data)){
    if(!is.null(which(data$blood_pressure_lowering_med=="TRUE"))){
      binarized_phenotypes[which(binarized_phenotypes$blood_pressure_lowering_med=="TRUE"),"blood_pressure_lowering_med"]<-1
      binarized_phenotypes[which(binarized_phenotypes$blood_pressure_lowering_med=="FALSE"),"blood_pressure_lowering_med"]<-0
    }
    binarized_phenotypes$blood_pressure_lowering_med<-as.factor(binarized_phenotypes$blood_pressure_lowering_med)
  }
  #alcohol_consumption
  if(c("alcohol_consumption") %in% colnames(data)){
    if(!is.null(which(data$alcohol_consumption=="TRUE"))){
      binarized_phenotypes[which(binarized_phenotypes$alcohol_consumption=="TRUE"),"alcohol_consumption"]<-1
      binarized_phenotypes[which(binarized_phenotypes$alcohol_consumption=="FALSE"),"alcohol_consumption"]<-0
    }
    binarized_phenotypes$alcohol_consumption<-as.factor(binarized_phenotypes$alcohol_consumption)
  }
  #metabolic_syndrome
  if(c("metabolic_syndrome") %in% colnames(data)){
    if(!is.null(which(data$metabolic_syndrome=="TRUE"))){
      binarized_phenotypes[which(binarized_phenotypes$metabolic_syndrome=="TRUE"),"metabolic_syndrome"]<-1
      binarized_phenotypes[which(binarized_phenotypes$metabolic_syndrome=="FALSE"),"metabolic_syndrome"]<-0
    }
    binarized_phenotypes$metabolic_syndrome<-as.factor(binarized_phenotypes$metabolic_syndrome)
  }
  
  #hscrp
  if(c("hscrp") %in% colnames(data)){
    binarized_phenotypes[which(data$hscrp<3),"high_hscrp"]<-0
    binarized_phenotypes[which(data$hscrp>=3),"high_hscrp"]<-1
    binarized_phenotypes$high_hscrp<-as.factor(binarized_phenotypes$high_hscrp)
  }else if(c("ln_hscrp") %in% colnames(data)){
    hscrp<-exp(data[,"ln_hscrp"])
    binarized_phenotypes[which(hscrp<3),"high_hscrp"]<-0
    binarized_phenotypes[which(hscrp>=3),"high_hscrp"]<-1
    binarized_phenotypes$high_hscrp<-as.factor(binarized_phenotypes$high_hscrp)
  }
  
  #obesity
  if(all(c("BMI","waist_circumference","sex") %in% colnames(data))){
    binarized_phenotypes[which(data$BMI>=30 & data$waist_circumference>=100 & data$sex==1),"obesity"]<-1
    binarized_phenotypes[which(data$BMI>=30 & data$waist_circumference>=93 & data$sex==0),"obesity"]<-1
    binarized_phenotypes[which(data$BMI<30 | data$waist_circumference<100 & data$sex==1),"obesity"]<-0
    binarized_phenotypes[which(data$BMI<30 | data$waist_circumference<93 & data$sex==0),"obesity"]<-0
    binarized_phenotypes[which(is.na(data$BMI) | is.na(data$waist_circumference) | is.na(data$sex)),"obesity"]<-NA
    binarized_phenotypes$obesity<-as.factor(binarized_phenotypes$obesity)
  }
  #high age
  if(all(c("age") %in% colnames(data))){
    binarized_phenotypes[which(data$age>=65),"high_age"]<-1
    binarized_phenotypes[which(data$age<65),"high_age"]<-0
    binarized_phenotypes$high_age<-as.factor(binarized_phenotypes$high_age)
  }
  #middle age
  if(all(c("age") %in% colnames(data))){
    binarized_phenotypes[which(data$age>=45 & data$age<65),"middle_age"]<-1
    binarized_phenotypes[which(data$age<45),"middle_age"]<-0
    binarized_phenotypes[which(data$age>=65),"middle_age"]<-0
    binarized_phenotypes$middle_age<-as.factor(binarized_phenotypes$middle_age)
  }
  #low age
  if(all(c("age") %in% colnames(data))){
    binarized_phenotypes[which(data$age<=45),"low_age"]<-1
    binarized_phenotypes[which(data$age>45),"low_age"]<-0
    binarized_phenotypes$low_age<-as.factor(binarized_phenotypes$low_age)
  }
  #high triglycerides
  if(all(c("triglycerides") %in% colnames(data))){
    binarized_phenotypes[which(data$triglycerides>=2.3),"high_triglycerides"]<-1
    binarized_phenotypes[which(data$triglycerides<2.3),"high_triglycerides"]<-0
    binarized_phenotypes$high_triglycerides<-as.factor(binarized_phenotypes$high_triglycerides)
  }
  #high ldl cholesterol
  if(all(c("ldl_chol") %in% colnames(data))){
    binarized_phenotypes[which(data$ldl_chol>=4.1),"high_ldl_chol"]<-1
    binarized_phenotypes[which(data$ldl_chol<4.1),"high_ldl_chol"]<-0
    binarized_phenotypes$high_ldl_chol<-as.factor(binarized_phenotypes$high_ldl_chol)
  }
  #low hdl cholesterol
  if(all(c("hdlchol") %in% colnames(data))){
    binarized_phenotypes[which(data$hdlchol<=1.3),"low_hdlchol"]<-1
    binarized_phenotypes[which(data$hdlchol>1.3),"low_hdlchol"]<-0
    binarized_phenotypes$low_hdlchol<-as.factor(binarized_phenotypes$low_hdlchol)
  }
  #high total cholesterol
  if(all(c("totchol") %in% colnames(data))){
    binarized_phenotypes[which(data$totchol>=6.2),"high_totchol"]<-1
    binarized_phenotypes[which(data$totchol<6.2),"high_totchol"]<-0
    binarized_phenotypes$high_totchol<-as.factor(binarized_phenotypes$high_totchol)
  }
  #low eGFR
  if(all(c("eGFR") %in% colnames(data))){
    binarized_phenotypes[which(data$eGFR<=60),"low_eGFR"]<-1
    binarized_phenotypes[which(data$eGFR>60),"low_eGFR"]<-0
    binarized_phenotypes$low_eGFR<-as.factor(binarized_phenotypes$low_eGFR)               
  }
  #low wbc
  if(all(c("wbc") %in% colnames(data))){
    binarized_phenotypes[which(data$wbc<=4.5),"low_wbc"]<-1
    binarized_phenotypes[which(data$wbc>4.5),"low_wbc"]<-0
    binarized_phenotypes$low_wbc<-as.factor(binarized_phenotypes$low_wbc)
  }
  #low hgb
  if(all(c("hgb","sex") %in% colnames(data))){
    binarized_phenotypes[which(data$hgb<=8.67 & data$sex==1),"low_hgb"]<-1
    binarized_phenotypes[which(data$hgb<=7.62 & data$sex==0),"low_hgb"]<-1
    binarized_phenotypes[which(data$hgb>8.67 & data$sex==1),"low_hgb"]<-0
    binarized_phenotypes[which(data$hgb>7.62 & data$sex==0),"low_hgb"]<-0
    binarized_phenotypes$low_hgb<-as.factor(binarized_phenotypes$low_hgb)         
  }
  #factorize all variables
  binarized_phenotypes <- data.frame(sapply(binarized_phenotypes[,colnames(binarized_phenotypes)] , factor))
  
  rownames(binarized_phenotypes)<-rownames(data)
  return(binarized_phenotypes)
}

#' QCprep_surrogates
#' 
#' Helper function to pre-process the Nightingale Health metabolomics data-set before applying metabolomics-based surrogates by Bizzarri et al.
#'
#' @param mat numeric data-frame Nightingale metabolomics matrix.
#' @param PARAM_surrogates is a list holding the parameters to compute the surrogates
#' @param Nmax_zero numeric value indicating the maximum number of zeros allowed per sample (Number suggested=1)
#' @param Nmax_miss numeric value indicating the maximum number of missing values allowed per sample (Number suggested=1)
#' @param quiet logical to suppress the messages in the console
#' @return Nightingale-metabolomics data-frame after pre-processing (checked for zeros, missing values, samples>5SD from the BBMRI-mean, imputing the missing values and z-scaled)
#' 
#' @examples
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' #Pre-process the metabolic features
#' prepped_met<-QCprep_surrogates(as.matrix(metabolic_measures), MiMIR::PARAM_surrogates)
#' 
#' 
#' @details 
#' Bizzarri et al. built multivariate models,using 56 metabolic features quantified by Nightingale, to predict the 19 binary characteristics of an individual. 
#' The binary variables are: sex, diabetes status, metabolic syndrome status, lipid medication usage, blood pressure lowering medication,
#' current smoking, alcohol consumption, high age, middle age, low age, high hsCRP, high triglycerides, high ldl cholesterol,
#' high total cholesterol, low hdl cholesterol, low eGFR, low white blood cells, low hemoglobin levels.
#' 
#' @references 
#' This function was made to vidualize the binarized variables calculated following the rules indicated in the article:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' 
#' @seealso 
#' binarize_all_pheno
#' 
#' @export
#'
QCprep_surrogates<-function(mat,PARAM_surrogates, Nmax_miss=1, Nmax_zero=1,quiet=FALSE){
  ## mat is the input NH matrix; it may contain a mixture of flags and metabolites in the columns.
  ## PARAM is a list holding the parameters of the pipeline
  # 0. Start:
  if(!quiet){
    cat(report.dim(mat,header="Start"))
  }
  # 1. Subset required metabolites:
  mat <- subset_metabolites_overlap(mat,metabos=PARAM_surrogates$MET,quiet=quiet)
  # 2. Subset samples on missingness:
  mat <- subset_samples_miss(mat,Nmax=Nmax_miss,quiet=quiet)
  # 3. Subset samples on zeros:
  mat <- subset_samples_zero(mat,Nmax=Nmax_zero,quiet=quiet)
  # 4. Subset samples on SD:
  mat <- subset_samples_sd_surrogates(as.matrix(mat),MEAN=as.matrix(PARAM_surrogates$mean),
                                  SD=as.matrix(PARAM_surrogates$sd),quiet=quiet)
  sample_names<-rownames(mat)
  
  # 5. Perform scaling:
  mat <- sapply(PARAM_surrogates$MET,function(x) (mat[,x]-as.numeric(PARAM_surrogates$mean[x]))/as.numeric(PARAM_surrogates$sd[x]))
  rownames(mat)<-sample_names
  if(!quiet){
    cat("| Performing scaling ... ")
    cat(" DONE!\n")
  }
  
  # 6. Perform imputation:
  mat <- impute_miss(mat)
  if(!quiet){
  cat("| Imputation ... ")
  cat(" DONE!\n")
  }
  return(mat)
}


#' calculate_surrogate_scores
#' 
#' Function to compute the surrogate scores by Bizzarri et al. from  the Nightingale metabolomics matrix
#'
#' @param met numeric data-frame with Nightingale-metabolomics
#' @param pheno phenotypic data.frame including this clinical variables (with the same nomenclature): "sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking",
#' "metabolic_syndrome", "alcohol_consumption", "age","BMI", "ln_hscrp","waist_circumference",
#' "weight","height", "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb" 
#' @param PARAM_surrogates list containing the parameters to compute the metabolomics-based surrogates
#' @param bin_names vector of strings containing the names of the binary variables
#' @param Nmax_zero numeric value indicating the maximum number of zeros allowed per sample (Number suggested=1)
#' @param Nmax_miss numeric value indicating the maximum number of missing values allowed per sample (Number suggested=1)
#' @param quiet logical to suppress the messages in the console
#' @param post logical to indicate if the function should calculate the posterior probabilities
#' @param roc logical to plot ROC curves for the metabolomics surrogate (available only for the phenotypes included)
#' @return if pheno is not available: list with the surrogates and the Nightingale metabolomics matrix after QC.
#' if pheno is available: list with the surrogates, ROC curves, phenotypes, binarized phenotypes and the Nightingale metabolomics matrix after QC,
#' 
#' @examples
#' require(MiMIR)
#' require(foreach)
#' require(pROC)
#' require(foreach)
#' 
#' #load dataset
#' m <- synthetic_metabolic_dataset
#' p <- synthetic_phenotypic_dataset
#' #Apply the surrogates
#' sur<-calculate_surrogate_scores(met=m,pheno=p,MiMIR::PARAM_surrogates,bin_names=c("sex","diabetes"))
#' 
#' @details 
#' Bizzarri et al. built multivariate models,using 56 metabolic features quantified by Nightingale, to predict the 19 binary characteristics of an individual. 
#' The binary variables are: sex, diabetes status, metabolic syndrome status, lipid medication usage, blood pressure lowering medication,
#' current smoking, alcohol consumption, high age, middle age, low age, high hsCRP, high triglycerides, high ldl cholesterol,
#' high total cholesterol, low hdl cholesterol, low eGFR, low white blood cells, low hemoglobin levels.
#' 
#' @references 
#' This function was made to vidualize the binarized variables calculated following the rules indicated in the article:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' 
#' @seealso 
#' QCprep_surrogates
#'  
#' @export
#' 
calculate_surrogate_scores <- function(met, pheno, PARAM_surrogates, bin_names=c("sex","diabetes"), Nmax_miss=1,Nmax_zero=1, post=TRUE, roc=FALSE,quiet=FALSE){
  bin_surro=paste0("s_",bin_names)
  #QC of the metabolites
  metabo_measures<-QCprep_surrogates(as.matrix(met[,MiMIR::metabolites_subsets$MET63]), 
                                     PARAM_surrogates, quiet=quiet,Nmax_miss=Nmax_miss,Nmax_zero=Nmax_zero)
  
  if(roc){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar)) 
    
    phenotypes<-pheno[which(rownames(pheno)%in%rownames(metabo_measures)),]
    #Calculating the binarized surrogates
    bin_pheno<-binarize_all_pheno(phenotypes)
    all(rownames(bin_pheno)==rownames(metabo_measures))
    
    #Surrogates Calculation
    surrogates<-foreach::foreach(i=bin_names, .combine="cbind") %do% {
      pred<-apply.fit_surro(as.matrix(metabo_measures), 
                            PARAM_surrogates$models_betas[paste0("s_",i),])
    }
    #ROC curves of the available clinical variables
    ind<-sapply(bin_names, function(x) 
      !all(is.na(bin_pheno[,x])) && length(levels(bin_pheno[,x]))==2
    )
    graphics::par(mfrow=c(5,4))
    surro_images<-foreach::foreach(i=bin_names[ind], .combine="cbind") %do% {
      a<-data.frame(out=factor(bin_pheno[,i]), metabo_measures)
      colnames(a)[1]<-"out"
      pred<-predictions_surrogates(PARAM_surrogates$models_betas[paste0("s_",i),], data=a,title_img=i)
      return(pred$predictions_eval)
    }
  }else{
      #Surrogates Calculation
      surrogates<-foreach::foreach(i=bin_names, .combine="cbind") %do% {
        pred<-apply.fit_surro(as.matrix(metabo_measures),PARAM_surrogates$models_betas[paste0("s_",i),], post=post)
      }
  }
  colnames(surrogates)<-bin_surro
  rownames(surrogates)<-rownames(metabo_measures)
  #Outputs
  if(roc){
    return(list(dat=list(phenotypes=phenotypes,metabo_measures=metabo_measures,bin_phenotypes=bin_phenotypes),surrogates=surrogates, roc_curves=surro_images))
  }else{
    return(list(metabo_measures=met,surrogates=surrogates))
  }
}

####################
## Plot functions ##
####################
#' pheno_barplots
#' 
#' #' Function created to binarize the phenotypes used to calculate the metabolomics based surrogate made by Bizzarri et al.
#'
#' @param bin_phenotypes phenotypes data.frame containing some of the following variables (with the same namenclature):
#' "sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking",
#' "metabolic_syndrome", "alcohol_consumption", "age","BMI", "ln_hscrp","waist_circumference",
#' "weight","height", "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb" 
#' @return The phenotypic variables binarized following the thresholds in in the metabolomics surrogates made by by Bizzarri et al.
#' 
#' @examples
#' require(MiMIR)
#' require(foreach)
#' 
#' #load the phenotypes dataset
#' phenotypes <- synthetic_phenotypic_dataset
#' 
#' #Calculate BMI, LDL cholesterol and eGFR
#' binarized_phenotypes<-binarize_all_pheno(phenotypes)
#' #Plot the variables
#' pheno_barplots(binarized_phenotypes)
#' 
#' 
#' @details 
#' Bizzarri et al. built multivariate models,using 56 metabolic features quantified by Nightingale, to predict the 19 binary characteristics of an individual. 
#' The binary variables are: sex, diabetes status, metabolic syndrome status, lipid medication usage, blood pressure lowering medication,
#' current smoking, alcohol consumption, high age, middle age, low age, high hsCRP, high triglycerides, high ldl cholesterol,
#' high total cholesterol, low hdl cholesterol, low eGFR, low white blood cells, low hemoglobin levels.
#' 
#' @references 
#' This function was made to vidualize the binarized variables calculated following the rules indicated in the article:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' 
#' @seealso 
#' binarize_all_pheno
#'  
#' @export
#'
pheno_barplots<-function(bin_phenotypes){
  molten_pheno <- foreach::foreach(i=colnames(bin_phenotypes), .combine = "rbind") %do% {data.frame(ID=rownames(bin_phenotypes), Variable=i, Value=bin_phenotypes[,i])}
  
  tab<-as.data.frame(stats::xtabs( ~ Variable + Value, molten_pheno))
  tab1<-data.frame(Variable=tab[which(tab$Value==1),"Variable"],
                   true=tab[which(tab$Value==1),"Freq"],
                   false=tab[which(tab$Value==0),"Freq"])
  
  fig <- plotly::plot_ly(tab1, x = ~Variable, y = ~true, type = 'bar', name="1", marker = list(color = 'rgb(49,130,189)'))
  fig <- fig %>% plotly::add_trace(y = ~false, name = '0', marker = list(color = 'rgb(204,204,204)'))
  fig <- fig %>% plotly::layout(xaxis = list(title = "", tickangle = -45),
                                yaxis = list(title = ""),
                                margin = list(b = 100),
                                barmode = 'group')
  fig
}


#' plot_na_heatmap
#' 
#' Function plotting information about missing & zero values on the indicated matrix.
#'
#' @param dat The matrix or data.frame
#' @return Plot with a central heatmap and two histogram on the sides
#
#' @examples
#' library(graphics)
#' library(MiMIR)
#' 
#' #load the metabolites dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' #Plot the missing values in the metabolomics matrix
#' plot_na_heatmap(metabolic_measures)
#' 
#' @details
#' This heatmap indicates the available values in grey and missing or zeros in white. 
#' On the sides two bar plots on the sides, one showing the missingn or zero values 
#' per row and another to show the missing or zeroes per column.
#'  
#' @export
#'
plot_na_heatmap  <- function(dat){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  
  Dummi <- matrix(NA,ncol=ncol(dat),nrow=nrow(dat))
  Dummi[which(!is.na(dat==0))] <- 1
  graphics::layout(mat=matrix(c(1,2,3,4),ncol=2,byrow=TRUE),widths=c(4,1),heights=c(1,4))
  graphics::par(xaxs = "i")  # THIS IS EVIL MAGIC!
  graphics::par(yaxs = "i")
  graphics::par(mar=c(0,3,1,0))
  XCOUNT <- nrow(Dummi)-colSums(Dummi,na.rm=TRUE)
  XPERC  <- XCOUNT/nrow(Dummi)*100
  YCOUNT <- ncol(Dummi)-rowSums(Dummi,na.rm=TRUE)
  YPERC  <- YCOUNT/ncol(Dummi)*100
  if(max(XPERC)<10){
    ylim <- c(0,10)
  } else {
    ylim <- c(0,max(pretty(XPERC)))
  }
  graphics::par(mar=c(0.5,13,2.5,0))
  graphics::barplot(XPERC,axes=TRUE,main="Missingness Plot",col="lightblue",border="lightblue",
          las=2,cex.main=3,ylim=ylim,cex.axis=1.2,font.axis=1.5)
  graphics::par(mar=c(0,0,0,0))
  graphics::plot.new()
  graphics::plot.window(xlim=c(0,1),ylim=c(0,1))
  graphics::legend("topright",fill=c("white","grey30"),legend=c("missing","value"))
  graphics::par(mar=c(2,13,0,0))
  graphics::image(t(Dummi[nrow(Dummi):1,]),axes=F, col=c("grey30"),xlim=c(0,1),ylim=c(0,1))
  graphics::mtext(text=rownames(dat), side=2, line=0.3, at=c(0.985,seq(0.945,0.055, length=dim(dat)[1]-2),0.03),
        las=1, cex=0.8)
  graphics::mtext(text="samples",side=1,font=1.5,cex=1.5, line = 0.6)
  graphics::par(mar=c(2,0,0,0))
  if(max(YPERC)<10){
    xlim <- c(0,10)
  } else {
    xlim <- c(0,max(pretty(YPERC)))
  }
  graphics::barplot(rev(YPERC),horiz=TRUE,axes=TRUE,col="lightblue",border="lightblue",cex.axis=1,font.axis=2)
}

#' roc_surro
#' 
#' Function that creates a ROC curve of the selected metabolic surrogates as a plotly image
#'
#' @param surrogates numeric data.frame of metabolomics-based surrogate values by Bizzarri et al.
#' @param bin_phenotypes logic data.frame of binarized phenotypes
#' @param x_name vector of strings with the names of the selected binary phenotypes for the roc
#' @return plotly image with the ROC curves for one or more selected variables
#' 
#' @examples
#' require(pROC)
#' require(plotly)
#' require(foreach)
#' require(MiMIR)
#' 
#' #load the dataset
#' met <- synthetic_metabolic_dataset
#' phen<- synthetic_phenotypic_dataset
#' 
#' #Calculating the binarized surrogates
#' b_phen<-binarize_all_pheno(phen)
#' #Apply a surrogate models and plot the ROC curve
#' surr<-calculate_surrogate_scores(met, phen, MiMIR::PARAM_surrogates, colnames(b_phen))
#' #Plot the ROC curves
#' roc_surro(surr$surrogates, b_phen, "sex")
#' 
#' @export
#'
roc_surro<-function(surrogates, bin_phenotypes, x_name){
  ## Plots ##
  #Layouts
  axis_font <- list(
    family = "Arial",
    size = 80
  )
  tick_font <- list(
    family = "Arial",
    size = 6
  )
  title_font <- list(
    family = "Arial",
    size = 22,
    margin=10
  )
  
  available<-colnames(bin_phenotypes)[!colSums(is.na(bin_phenotypes))==nrow(bin_phenotypes)]
  non_available<-colnames(bin_phenotypes)[colSums(is.na(bin_phenotypes))==nrow(bin_phenotypes)]
  selected<-available[which(available %in% x_name)]
  sel_not_available<-non_available[which(non_available %in% x_name)]
  
  surro<-paste0("s_",selected[1])
  
  ROC_curve<-pROC::roc(bin_phenotypes[,selected[1]], as.numeric(surrogates[,paste0("s_",selected[1])]), plot=F, col=MiMIR::c21[surro], quiet = TRUE,
                 lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
  
  roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
  pl <- plotly::plot_ly(roc_df, x= ~Specificity, y= ~Sensitivity, 
                type = 'scatter', mode = 'lines', line = list(color= MiMIR::c21[surro],width = 4),
                name=paste0(surro, ":\nAUC=",round(as.numeric(ROC_curve$auc),digits = 3))) %>%
    plotly::add_segments(x = 1, xend =0, y = 0, yend = 1, line = list(dash = "dash", width= 1, color="black"),
                 name="random:\nAUC=0.5")
  
  if(length(selected)>1){
    surro_images<-foreach::foreach(i=2:length(selected), .combine="cbind") %do% {
      surro<-paste0("s_",selected[i])
  
      ROC_curve<-pROC::roc(bin_phenotypes[,selected[i]], as.numeric(surrogates[,paste0("s_",selected[i])]), plot=F, col=MiMIR::c21[surro], quiet = TRUE,
                     lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
      roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
      
      pl<-pl %>% 
        plotly::add_trace(x= roc_df$Specificity, y= roc_df$Sensitivity, 
                           type = 'scatter', mode = 'lines', line = list(color= MiMIR::c21[surro],width = 4),
                           name=paste0(selected[i], ":\nAUC=",round(as.numeric(ROC_curve$auc),digits = 3)))
    }
  }
  
    if(length(sel_not_available)>0){
      if(length(sel_not_available)>1){
        ann<-paste("Not available:\n",noquote(paste(sel_not_available, collapse = ',\n')))
      }else{ann<-paste("Not available:\n",sel_not_available)}
      
      pl<-pl %>% 
        plotly::layout(
        title = list(text=paste0('<b> ROC curves <b>'),font=title_font, y = 0.98),
        xaxis = list(autorange = "reversed", 
                     title = paste("<b> Specificity <b>")),
        yaxis = list (title = paste("<b> Sensitivity <b>")),
        annotations = list(
          text = ann,
          xref = "paper",
          yref = "paper",
          xanchor = c("center"),
          align = "center",
          x = c(1.12),
          y = c(0.5),
          showarrow=F,
          font=list(size=12))
      )
    }else{
      pl <- pl %>% 
        plotly::layout(
        title = list(text=paste0('<b> ROC curves <b>'),font=title_font, y = 0.98),
        xaxis = list(autorange = "reversed", 
                     title = paste("<b> Specificity <b>")),
        yaxis = list (title = paste("<b> Sensitivity <b>"))
      )
    }
    return(pl)
}

#' roc_surro_subplots
#' 
#' Function that plots the ROCs of the surrogates of all the available surrogate models as plotly sub-plots
#'
#' @param surrogates numeric data.frame containing the surrogate values by Bizzarri et al.
#' @param bin_phenotypes numeric data.frame with the binarized phenotypes output of binarize_all_pheno
#' @return plotly image with all the ROCs for all the available clinical variables
#' 
#' @examples
#' library(pROC)
#' library(plotly)
#' library(MiMIR)
#' 
#' #load the dataset
#' met <- synthetic_metabolic_dataset
#' phen<- synthetic_phenotypic_dataset
#' 
#' #Calculating the binarized surrogates
#' b_phen<-binarize_all_pheno(phen)
#' #Apply a surrogate models and plot the ROC curve
#' surr<-calculate_surrogate_scores(met, phen, MiMIR::PARAM_surrogates, colnames(b_phen))
#' 
#' roc_surro_subplots(surr$surrogates, b_phen)
#' 
#' 
#' @export
#'
roc_surro_subplots<-function(surrogates, bin_phenotypes){
  
  #Layouts
  axis_font <- list(
    family = "Arial",
    size = 80
  )
  tick_font <- list(
    family = "Arial",
    size = 6
  )
  title_font <- list(
    family = "Arial",
    size = 22,
    margin=10
  )
  
  available<-colnames(bin_phenotypes)[!colSums(is.na(bin_phenotypes))==nrow(bin_phenotypes)]
  surro_images<-lapply(available, function(i){
    surro<-paste0("s_",i)
    ROC_curve<-pROC::roc(bin_phenotypes[,i], as.numeric(surrogates[,surro]), plot=F, col=MiMIR::c21[surro], quiet = TRUE,
                   lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
    roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
    
    pl <- plotly::plot_ly(roc_df, x= ~Specificity, y= ~Sensitivity, 
                  type = 'scatter', mode = 'lines', line = list(color= MiMIR::c21[surro],width = 4)) %>%
      plotly::add_segments(x = 1, xend =0, y = 0, yend = 1, line = list(dash = "dash", width= 1, color="black")) %>% 
      plotly::layout(
        xaxis = list(autorange = "reversed", tickfont = tick_font),
        yaxis = list (tickfont = tick_font),
        margin = list(l = 10, r = 10, b = 7, t = 7),
        autosize=FALSE,
        width=770,
        height=500,
        showlegend = F,
        annotations = list(
          text = list(paste0('<b>',i,'<b>, AUC=<b>',round(as.numeric(ROC_curve$auc),digits = 2)),
                      '<b>Specificity<b>','<b>Sensitivity<b>'),
          xref = "paper",
          yref = "paper",
          yanchor = c("bottom","top","bottom"),
          xanchor = c("center","center","center"),
          align = "center",
          x = c(0.5,0.5,-0.07),
          y = c(0.923,-0.02, 0.2),
          showarrow=F,
          textangle=c(0,0,-90),
          font=list(size=7.5)
        )
      )
  })
  return(plotly::subplot(surro_images, nrows = 5, shareX = F,shareY = F, titleX = F, titleY = F, which_layout = 1))
}

#' ttest_surrogates
#' 
#' Function that calculates a t-test and a plotly image of the selected surrogates
#' @param surrogates numeric data.frame containing the surrogate values by Bizzarri et al.
#' @param bin_phenotypes numeric data.frame with the binarized phenotypes output of binarize_all_pheno
#' @return plotly image with all the ROCs for all the available clinical variables
#' 
#' @examples
#' require(pROC)
#' require(plotly)
#' require(MiMIR)
#' require(foreach)
#' 
#' #load the dataset
#' m <- synthetic_metabolic_dataset
#' p <- synthetic_phenotypic_dataset
#' 
#' #Calculating the binarized surrogates
#' b_p<-binarize_all_pheno(p)
#' #Apply a surrogate models and plot the ROC curve
#' surr<-calculate_surrogate_scores(met=m, pheno=p, MiMIR::PARAM_surrogates, bin_names=colnames(b_p))
#' ttest_surrogates(surr$surrogates, b_p)
#' 
#' @details 
#' Barplot and T-test indicating if the surrogate variables could split accordingly the real value of the binary clinical variables. 
#' @export
#'
ttest_surrogates<-function(surrogates,bin_phenotypes){
  available<-colnames(bin_phenotypes)[!colSums(is.na(bin_phenotypes))==nrow(bin_phenotypes)]
  
  Surrogates<-foreach::foreach(i=colnames(bin_phenotypes), .combine="rbind") %do% {
    surro<-paste0("s_",i)
    comp<-data.frame(ID=rownames(surrogates), 
                     value=bin_phenotypes[rownames(surrogates),i],
                     surrogate=surrogates[,surro], variable=i, 
                     stringsAsFactors = F)
  }
  
  tt_values<-foreach::foreach(i=available, .combine='rbind') %do%{
    comp<-Surrogates[which(Surrogates$variable==i),]
    if(length(which(is.na(comp[,"value"])))>0){
      comp<-comp[-which(is.na(comp[,"value"])),]
    }
    comp$value<-factor(comp$value)
    tt<-stats::t.test(comp[which(comp$value==1),"surrogate"],comp[which(comp$value==0),"surrogate"], var.equal = TRUE )
    ptext<-paste(ifelse(tt$p.value >= 0.05, "pval>0.05", 
                        ifelse( 0.01<= tt$p.value && tt$p.value <  0.05, "pval<0.05", 
                                ifelse(0.001<= tt$p.value && tt$p.value< 0.01, "pval<0.01", 
                                       ifelse(tt$p.value< 0.001, "pval<0.001")))))
    
    return(ptext)
  }
  
  rownames(tt_values)<-available
  
  for (i in available){
    Surrogates[which(Surrogates$variable==i),"variable"]<-paste0(i,", ",tt_values[i,])
  }
  Surrogates$value<-as.numeric(Surrogates$value)-1
  Surrogates[which(is.na(Surrogates),arr.ind=T)]<-NaN
  Surrogates$value<-factor(Surrogates$value)
  pl_surro<-suppressWarnings(plotly::plot_ly(Surrogates,
              x = ~variable,
              y = ~surrogate,
              color = ~value,
              type = "box",
              colors = c("#377EB8","#E41A1C", "grey")
              ) %>% 
                plotly::layout(boxmode = "group",
           title = list(text="<b>Surrogates' distributions split for their original values<b>",y = 0.98),
           xaxis = list(title = "<b>Clinical Variables<b>",
                        zeroline = FALSE),
           yaxis = list(title = "<b>Surrogate values<b>",
                        zeroline = FALSE)
           ))
  
  return(pl_surro)
}

#' LOBOV_accuracies
#' 
#' Function created to visualize the accuracies in the current dataset compared to the
#' accuracies in the Leave One Biobank Out Validation in Bizzarri et al.
#'
#' @param surrogates numeric data.frame containing the surrogate values by Bizzarri et al.
#' @param bin_phenotypes numeric data.frame with the binarized phenotypes output of binarize_all_pheno
#' @param bin_pheno_available vector of strings with the available phenotypes
#' @param acc_LOBOV accuracy of LOBOV calculated in Bizzarri et al.
#' @return Boxplot with the accuracies of the LOBOV
#' 
#' @examples
#' require(pROC)
#' require(plotly)
#' require(MiMIR)
#' require(foreach)
#' require(ggplot2)
#' 
#' #load the dataset
#' m <- synthetic_metabolic_dataset
#' p<- synthetic_phenotypic_dataset
#' 
#' #Calculating the binarized surrogates
#' b_p<-binarize_all_pheno(p)
#' #Apply a surrogate models and plot the ROC curve
#' sur<-calculate_surrogate_scores(m, p, MiMIR::PARAM_surrogates, bin_names=colnames(b_p))
#' p_avail<-colnames(b_p)[c(1:5)]
#' LOBOV_accuracies(sur$surrogates, b_p, p_avail, MiMIR::acc_LOBOV)
#' 
#' @details 
#' Comparison of the AUCs of the surrogates in the updated dataset and the 
#' results of the Leave One Biobank Out Validation made in BBMRI-nl.
#'
#' @references 
#' This function was made to vidualize the binarized variables calculated following the rules indicated in the article:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' 
#' @export
#'
LOBOV_accuracies<-function(surrogates, bin_phenotypes, bin_pheno_available, acc_LOBOV){
  AUCs<- foreach::foreach(i=1:length(bin_pheno_available), .combine="rbind") %do%{
    ROC_curve<-pROC::roc(bin_phenotypes[,bin_pheno_available[i]], as.numeric(surrogates[,paste0("s_",bin_pheno_available[i])]), plot=F, 
                         col=MiMIR::c21[surro], quiet = TRUE,
                         lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
    pROC::auc(ROC_curve)
  }
  
  rownames(AUCs)<-bin_pheno_available
  colnames(AUCs)<-"AUC"
  
  acc_bin<-foreach::foreach(i= 1:length(bin_pheno_available), .combine = "rbind") %do%{
    a<-acc_LOBOV[[bin_pheno_available[i]]][-c(dim(acc_LOBOV[[bin_pheno_available[i]]])[1],dim(acc_LOBOV[[bin_pheno_available[i]]])[1]-1),]
    a<-cbind(a,Uploaded="BBMRI-nl biobanks")
    a<-rbind(a,uploaded_data=data.frame(AUC=AUCs[bin_pheno_available[i],], Percentages=NaN, N=NaN, Uploaded="uploaded dataset"))
    a$Uploaded<-as.factor(a$Uploaded)
    
    acc<-data.frame(biobank=rownames(a),a,outcome=rep(bin_pheno_available[i], length=dim(a)[1]))
  }
  acc_bin$biobank<-as.factor(acc_bin$biobank)
  
  p <- ggplot2::ggplot(data= acc_bin, ggplot2::aes(x=outcome, y=AUC, color = Uploaded)) +
    #ggplot2::geom_jitter(ggplot2::aes(text=paste("Biobank: ", biobank)), width=0.25, alpha=1) +
    ggplot2::geom_jitter(width=0.25, alpha=1) +
    ggplot2::geom_boxplot(outlier.size = 0,outlier.shape = NA) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=15, face= "bold"),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12))+
    ggplot2::labs(title="AUC of the current dataset compared to the LOBOV",x="Clinical variables", y="AUC")
  p
  
  # Need to modify the plotly object and make outlier points have opacity equal to 0
  fig <- plotly::ggplotly(p)
  
  fig$x$data <- lapply(fig$x$data, FUN = function(x){
    
    if (x$type == "box") {
      x$marker = list(opacity = 0)
    }
    return(x)
  })
  return(fig)
}



#' scatterplot_predictions
#' 
#' Function to visualize a scatter-plot comparing two variables
#'
#' @param x numeric vector
#' @param p second numeric vector
#' @param title string vector with the title
#' @param xname string vector with the name of the variable on the x axis
#' @param yname string vector with the name of the variable on the y axis
#' @return plotly image with the scatterplot
#' 
#' @examples
#' library(plotly)
#' #load the dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' phenotypes <- synthetic_phenotypic_dataset
#' 
#' #Pre-process the metabolic features
#' prepped_met<-QCprep(as.matrix(metabolic_measures), MiMIR::PARAM_metaboAge)
#' #Apply the metaboAge
#' metaboAge<-apply.fit(prepped_met, FIT=PARAM_metaboAge$FIT_COEF)
#' 
#' age<-data.frame(phenotypes$age)
#' rownames(age)<-rownames(phenotypes)
#' scatterplot_predictions(age, metaboAge, title="Chronological Age vs MetaboAge")
#' 
#' @export
#'
scatterplot_predictions <- function(x, p, title, xname="x", yname="predicted x") {
  x<-data.frame(ID=row.names(x), original=x)
  colnames(x)[2]<-"original"
  p<-data.frame(ID=row.names(p), predicted=p)
  colnames(p)[2]<-"predicted"
  df<-merge(x,p,by.x="ID", by.y = "ID")
    
    axis_font <- list(
      family = "Arial",
      size = 14
    )
    title_font <- list(
      family = "Arial",
      size = 18,
      margin=10
    )
    
    line.fmt = list(dash="solid", width = 1.5, color=NULL)
    fit <- stats::lm(predicted ~ original, data = df)
    cordf<-stats::cor(df$original, df$predicted, method = c("pearson"))
    
    cPlot <- plotly::plot_ly(df, x= ~original, y=~predicted, mode="markers", alpha=0.7) %>% 
      plotly::add_markers(y= ~predicted) %>% 
      plotly::add_lines(x=~original, y=stats::predict(fit), name="Regression line", 
                line = list( width= 3, color="#E41A1C")) %>% 
      plotly::add_segments(x = min(df$original), xend = max(df$original), 
                 y = min(df$predicted), yend = max(df$predicted), 
                 name="x=y", line = list(dash = "dash", width= 1, color="black"))%>%
      plotly::layout(title = list(text=paste('<b>',title,'<b>'),font=title_font, y = 0.98),
             xaxis = list(title = paste('<b>',xname,'<b>'),titlefont = axis_font),
             yaxis = list (title = paste('<b>',yname,'<b>'),titlefont = axis_font),
             # #margin = list(l = 50, r = 140, b = 30, t = 30, pad = 4),
             annotations = list(text = paste0("R= ", as.character(round(cordf, digits= 3)),", R^2= ", as.character(round((cordf)^2, digits= 3)),
                                  "\n med. error=", as.character(round(stats::median(abs(df$original - df$predicted)), digits=3))),
                                xref='paper',
                                yref='paper',
                                x = 1.1, y = 0.25,
                                showarrow=FALSE,
                                bordercolor=c("black"),
                                borderwidth=2
                                )
             )
    return(cPlot)
}



#' multi_hist
#' 
#' #' Function to plot the histograms for all the variables in dat
#'
#' @param dat data.frame or matrix with the variables to plot
#' @param color colors selected for all the variables
#' @param scaled logical to z-scale the variables
#' @return plotly image with the histograms for all the variables in dat
#' @export
#'
#' @examples
#' \donttest{
#' library(plotly)
#' library(MiMIR)
#' 
#' #load the dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' 
#' multi_hist(metabolic_measures[,MiMIR::metabolites_subsets$MET14], scaled=T)
#' }
#' 
#' @export
#'
#'
multi_hist<-function(dat, color=MiMIR::c21, scaled=FALSE){
  
  axis_font <- list(
    family = "Arial",
    size = 80
  )
  tick_font <- list(
    family = "Arial",
    size = 6
  )
  title_font <- list(
    family = "Arial",
    size = 22,
    margin=10
  )
  hist_images<-lapply(colnames(dat), function(i){
    plot<-plotly::plot_ly(x =~dat[,i],opacity = 1, type = "histogram", marker = list(color = color[i]), name = i)%>% 
      plotly::layout(
        xaxis = list(tickfont = tick_font),
        yaxis = list (tickfont = tick_font),
        showlegend = F,
        annotations = list(
          text = list(paste0('<b>',i,'<b>')),
        xref = "paper",
        yref = "paper",
        yanchor = c("bottom"),
        xanchor = c("center"),
        align = "center",
        x = c(0.5),
        y = c(0.933),
        showarrow=F,
        textangle=c(0),
        font=list(size=10))
        )
  })
  return(plotly::subplot(hist_images, nrows = 5, shareX = F,shareY = F, titleX = F, titleY = F, which_layout = 1))
}

#' hist_plots
#' 
#' #' Function to plot the histograms for all the variables in dat
#'
#' @param dat data.frame or matrix with the variables to plot
#' @param x_name string with the names of the selected variables in dat
#' @param color colors selected for all the variables
#' @param scaled logical to z-scale the variables
#' @param datatype a character vector indicating what data type is beeing plotted
#' @param main title of the plot
#' @return plotly image with the histograms of the selected variables
#'
#' @examples
#' require(MiMIR)
#' require(plotly)
#' require(matrixStats)
#' #load the metabolites dataset
#' m <- synthetic_metabolic_dataset
#' 
#' #Apply a surrogate models and plot the ROC curve
#' surrogates<-calculate_surrogate_scores(m, PARAM_surrogates=MiMIR::PARAM_surrogates, roc=FALSE)
#' #Plot the histogram of the surrogate sex values scaled 
#' hist_plots(surrogates$surrogates, x_name="s_sex", scaled=TRUE)
#'
#' @export
#'
hist_plots<-function(dat, x_name, color=MiMIR::c21, scaled=FALSE, datatype="metabolic score", main="Predictors Distributions"){
  dat<-as.matrix(dat[,x_name])
  colnames(dat)<-x_name
  if(scaled){
    if(length(x_name)>1){
      sd<-matrixStats::colSds(dat,na.rm = T)
      names(sd)<-colnames(dat)
      dat<-apply.scale(dat, MEAN = colMeans(dat,na.rm = T), SD=sd)
    }else{
      dat <-as.matrix((dat-mean(dat,na.rm=T))/sd(dat,na.rm=T))
      colnames(dat)<-x_name
    }
  }
  axis_font <- list(
    family = "Arial",
    size = 18
  )
  title_font <- list(
    family = "Arial",
    size = 20,
    margin = 10
  )
  plot<-plotly::plot_ly(x = ~dat[,x_name[1]],opacity = 0.45, type = "histogram", histnorm = "probability", marker = list(color =color[x_name[1]]), name = x_name[1])
  x_title<-datatype
   if(length(x_name)>1){
    for(i in c(2:length(x_name))){
      plot <- plot %>% 
        plotly::add_histogram(x = dat[,x_name[i]],  histnorm = "probability", marker = list(color = color[x_name[i]]), name = x_name[i])
    }
     x_title<-datatype
  }
  if(scaled){
    x_title<-paste(datatype,"scaled")
  }
  plot<-plot %>% 
    plotly::layout( title =list(text=paste("<b>",main,"<b>"),font=title_font,y = 0.98),
                         xaxis = list(
                           title = paste("<b>",x_title,"<b>"),
                           titlefont = axis_font), 
                         yaxis = list(
                           title = "<b>Frequency<b>",
                           titlefont = axis_font),
                         barmode = "overlay")
  print(plot)
}
  

#' multi_hist
#' 
#' Function to plot the ~60 metabolites used for the metabolomics-based scores and compare them to to their distributions in BBMRI-nl
#'
#' @param dat data.frame or matrix with the metabolites
#' @param x_name string with the name of the selected variable
#' @param color colors selected for all the variables
#' @param scaled logical to z-scale the variables
#' @param datatype a character vector indicating what data type is being plotted
#' @param main title of the plot
#' @return plotly image with the histogram of the selected variable compared to the distributions in BBMRI-nl
#'
#'
#' @examples
#' library(plotly)
#' library(MiMIR)
#' 
#' #load the metabolites dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' 
#' BBMRI_hist_plot(metabolic_measures, x_name="alb", scaled=TRUE)
#' 
#' @details 
#' This function plots the distribution of a metabolic feature in the uploaded dataset, compared to their distributions in BBMRI-nl.
#' The selection of features available is done following the metabolic scores features.
#' 
#' @references 
#' The selection of metabolic features available is the one selected by the papers:
#' Deelen,J. et al. (2019) A metabolic profile of all-cause mortality risk identified in an observational study of 44,168 individuals. Nature Communications, 10, 1-8, <doi:10.1038/s41467-019-11311-9>
#' Ahola-Olli,A.V. et al. (2019) Circulating metabolites and the risk of type 2 diabetes: a prospective study of 11,896 young adults from four Finnish cohorts. Diabetologia, 62, 2298-2309, <doi:10.1007/s00125-019-05001-w>
#' Wurtz,P. et al. (2015) Metabolite profiling and cardiovascular event risk: a prospective study of 3 population-based cohorts. Circulation, 131, 774-785, <doi:10.1161/CIRCULATIONAHA.114.013116>
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541-547, <doi:10.1161/CIRCGEN.119.002610>
#' 
#' @export
#'
BBMRI_hist_plot<-function(dat, x_name, color=MiMIR::c21, scaled=FALSE, datatype="metabolite", main="Comparison with the metabolites measures in BBMRI"){
  dat<-as.matrix(dat[,x_name])
  colnames(dat)<-x_name
  if(scaled){
    dat <-as.matrix((dat-mean(dat,na.rm=T))/stats::sd(dat,na.rm=T))
    colnames(dat)<-x_name
    BBMRI<- BBMRI_hist_scaled[[x_name]]
  }else{
    BBMRI<- BBMRI_hist[[x_name]]
  }
  axis_font <- list(
    family = "Arial",
    size = 18
  )
  title_font <- list(
    family = "Arial",
    size = 20,
    margin = 10
  )
  
  res<-graphics::hist(stats::na.omit(dat[,x_name[1]]), breaks = BBMRI$breaks, plot = F)
  
  plot<-plotly::plot_ly(x = res$breaks[-1], y=res$density, opacity = 0.60, type = 'bar', 
                        marker = list(color = color[x_name[1]]), name = x_name[1])
  plot <- plot %>% 
    plotly:: add_trace(y = ~BBMRI$density,opacity = 0.50, type = 'bar', 
                       marker = list(color = "grey"), name = paste("BBMRI",x_name[1]))
  
  x_title<-datatype
  if(scaled){
    x_title<-paste(datatype,"scaled")
  }
  plot<-plot %>% 
    plotly::layout( title =list(text=paste("<b>",main,"<b>"),font=title_font,y = 0.98),
                    xaxis = list(
                      title = paste("<b>",x_title,"<b>"),
                      titlefont = axis_font), 
                    yaxis = list(
                      title = "<b>Frequency<b>",
                      titlefont = axis_font),
                    barmode = "overlay")
  print(plot)
}


#' hist_plots_mortality
#' 
#' #' Function to plot the histogram of the mortality score separated for different age ranges as a plotly image
#'
#' @param mort_score data.frame containing the mortality score
#' @param phenotypes data.frame containing age
#' @return plotly image with the histogram of the mortality score separated in 3 age ranges
#' 
#' @examples
#' library(MiMIR)
#' library(plotly)
#' #' #load the dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' phenotypes <- synthetic_phenotypic_dataset
#' 
#' #Compute the mortality score
#' mortScore<-comp.mort_score(metabolic_measures,quiet=TRUE)
#' #Plot the mortality score histogram at different ages
#' hist_plots_mortality(mortScore, phenotypes)
#'
#' @export
#'
hist_plots_mortality<-function(mort_score,phenotypes){
  age<-data.frame(age=phenotypes[rownames(mort_score),"age"])
  age$age_ranges<-base::cut(age$age, breaks= c(0,45,65, Inf), labels = c("age<45","45>=age<65","age>=65"))
  
  mort_score$age_ranges<-age$age_ranges
  colnames(mort_score)[2]<-"age_ranges"
  mort_score$age_ranges<-factor(mort_score$age_ranges, levels=c("age<45","45>=age<65","age>=65"))
  
  mu<-foreach(i=levels(mort_score$age_ranges), .combine = "rbind") %do% {
    m<-base::mean(mort_score[ which(mort_score$age_ranges==i),"mortScore"],na.rm=T)
    m<-data.frame(age_ranges=i,mean=m)
    }
  
  axis_font <- list(
    family = "Arial",
    size = 14
  )
  title_font <- list(
    family = "Arial",
    size = 18,
    margin=10
  )
  fig <- plotly::plot_ly(alpha = 0.4) 
    if(length(which(mort_score$age_ranges=="age<45"))>0){
      fig<-fig %>% 
        plotly::add_histogram(x = ~mort_score$mortScore[which(mort_score$age_ranges=="age<45")], name="age<45", marker = list(color = "#B3DE69"), opacity=0.7) %>% 
        plotly::add_segments(x = mu[1,2], xend = mu[1,2], y = 0, yend = 100, name="mean(age<45)")
    }
  if(length(which(mort_score$age_ranges=="45>=age<65"))>0){
    fig<-fig %>% 
      plotly::add_histogram(x = ~mort_score$mortScore[which(mort_score$age_ranges=="45>=age<65")], name="45>=age<65", marker = list(color = "#80B1D3"), opacity=0.7) %>% 
      plotly::add_segments(x = mu[2,2], xend = mu[2,2], y = 0, yend = 100, name="mean(45>=age<65)")
  }
  if(length(which(mort_score$age_ranges=="age>=65"))>0){
    fig<-fig %>% 
      plotly::add_histogram(x = ~mort_score$mortScore[which(mort_score$age_ranges=="age>=65")], name="age>=65", marker = list(color = "#FB8072"), opacity=0.7) %>% 
      plotly::add_segments(x = mu[3,2], xend = mu[3,2], y = 0, yend = 100, name="mean(age>=65)")
  }  
  
    fig<-fig %>% 
      plotly::layout(title =list(text=paste("<b> Mortality divided in age ranges <b>"),font=title_font,y = 0.98),
           xaxis = list(
             title = paste("<b> Counts <b>"),
             titlefont = axis_font), 
           yaxis = list(
             title = "<b> Mortality score <b>",
             titlefont = axis_font),barmode = "overlay")
  fig
}


#' ttest_scores
#' 
#' #' Function that creates a boxplot with a continuous variable split using the binary variable
#'
#' @param dat The data.frame containing the 2 variables
#' @param pred character indicating the y variable
#' @param pheno character indicating the binary variable
#' @return plotly boxplot with the continuous variable split using the binary variable
#' @export
#' 
#' 
#' @examples
#' library(MiMIR)
#' library(plotly)
#' 
#' #load the dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' phenotypes <- synthetic_phenotypic_dataset
#' 
#' #Compute the mortality score
#' mortScore<-comp.mort_score(metabolic_measures,quiet=TRUE)

#' dat<-data.frame(predictor=mortScore, pheno=phenotypes$sex)
#' colnames(dat)<-c("predictor","pheno")
#' ttest_scores(dat = dat, pred= "mortScore", pheno="sex")
#'
#' @export
#'
ttest_scores<-function(dat, pred, pheno){
  
  dat[which(is.na(dat),arr.ind=T)]<-NaN
  dat$pheno<-factor(dat$pheno)

  tt<-stats::t.test(dat[which(dat$pheno==1),"predictor"],dat[which(dat$pheno==0),"predictor"], var.equal = TRUE )
  
  pl<-plotly::plot_ly(dat,
                            x = ~ pheno,
                            y = ~ predictor,
                            color = ~pheno,
                            type = "box",
                            colors = c("#377EB8","#E41A1C", "grey")
  ) %>% 
    plotly::layout(
                   title = list(text=paste0("<b>",pred," distribution split for ",pheno,", pvalue=",formatC(tt$p.value, format = "e", digits = 3),"<b>"),y = 0.98),
                   xaxis = list(title = paste("<b>",pheno,"<b>"),
                                zeroline = FALSE),
                   yaxis = list(title = paste("<b>",pred,"<b>"),
                                zeroline = FALSE)
    )
  
  return(pl)
}


#' kapmeier_scores
#' 
#' #' Function that creates a Kaplan Meier comparing first and last tertile of a metabolic score
#'
#' @param predictors The data.frame containing the predictors
#' @param pheno The data.frame containing the phenotypes
#' @param score a character string indicating which predictor to use
#' @param Eventname a character string with the name of the event to print on the plot
#' @return plotly with a Kaplan Meier comparing first and last tertile of a metabolic score
#' @export
#' 
#' @examples
#' require(MiMIR)
#' require(plotly)
#' require(survminer)
#' require(ggfortify)
#' require(ggplot2)
#' 
#' #load the dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' phenotypes <- synthetic_phenotypic_dataset
#'
#' #Compute the mortality score
#' mortScore<-comp.mort_score(metabolic_measures,quiet=TRUE)
#' 
#' #Plot a Kaplan Meier
#' kapmeier_scores(predictors=mortScore, pheno=phenotypes, score="mortScore")
#'
#' @export
#'
kapmeier_scores<-function(predictors, pheno, score, Eventname="Event"){
  dat<-data.frame(pheno[,c("age","Event","EventAge")], predictors[,score])
  colnames(dat)[4]<-"score"
  dat<-dat[which(!is.na(dat$age)),]
  dat<-dat[which(!is.na(dat$EventAge)),]
  dat<-dat[which(!is.na(dat$Event)),]
  dat<-dat[which(!is.na(dat[,"score"])),]
  if(length(which(dat$Event==-1))!=0){
    dat <-dat[-which(dat$Event==-1),]
  }
  
  dat$tertile <- with(dat, factor(
    findInterval(score, c(-Inf,
                         quantile(score, probs=c(0.3333333, 0.6666666)), Inf) ), 
    labels=c("Tertile 1","Tertile 2","Tertile 3")
  ))
  
  dat<-dat[-which(dat$tertile=="Tertile 2"),]
  
  km_trt_fit <- survminer::surv_fit(Surv(time = dat$EventAge-dat$age , event=dat$Event) ~ 
                                      tertile, data=dat)
  
  pl<-ggplot2::autoplot(km_trt_fit, xlab = "Time (in years)", ylab = "Survival",
                        main = paste0("<b>Kaplan Meier of ",score," associated with ", Eventname,"<b>",
                                     "\np=",formatC(survminer::surv_pvalue(km_trt_fit)[,2], format = "e", digits = 3),
                                     ", N=", dim(dat)[1],", Nevent=",sum(dat$Event)))
  
  
  pl<-plotly::ggplotly(pl) 
  
  return(pl)
}

#' cor_assoc
#' 
#' Function to calulate the correlation between 2 matrices
#'
#' @param dat1 matrix 1
#' @param dat2 matrix 2
#' @param feat1 vector of strings with the names of the selected variables in dat
#' @param feat2 vector if strings with the names of the selected variables in dat2
#' @param method indicates which methods of the correlation to use
#' @param quiet logical to suppress the messages in the console
#' @return correlations of the selected variables in the 2 martrices
#' 
#' @examples
#' library(stats)
#' 
#' #load the dataset
#' m <- as.matrix(synthetic_metabolic_dataset)
#' 
#' #Compute the pearson correlation of all the variables in the data.frame metabolic_measures
#' cors<-cor_assoc(m, m, MiMIR::metabolites_subsets$MET63,MiMIR::metabolites_subsets$MET63)
#' 
#' @seealso 
#' plot_corply
#' 
#' @export
#' 
cor_assoc <- function(dat1,dat2,feat1,feat2,method="pearson",quiet=FALSE){
  res <- lapply(feat1, function(F1){
    t(sapply(feat2,function(F2){
      unlist(stats::cor.test(x=dat1[,F1],y=dat2[,F2],method=method)[c("estimate","p.value","statistic")])
    }))
  })
  names(res) <- feat1
  return(res)
}


#' plot_corply
#' 
#' Function creating plottig the correlation between 2 datasets, dat1 x dat2 on basis of (partial) correlations
#'
#' @param res associations obtained with cor.assoc
#' @param main title of the plot
#' @param zlim max association to plot
#' @param reorder.x logical indicating if the function should reorder the x axis based on clustering
#' @param reorder.y logical indicating if the function should reorder the y axis based on clustering
#' @param resort_on_p logical indicating if the function should reorder x and y axis based on the pvalues of the associations
#' @param abs logical indicating if the function should  reorder based the absolute values
#' @param cor.abs logical indicating if the function should reorder the plot base on the absolute values 
#' @param reorder_dend Tlogical indicating if the function should reorder the plot based on dendrogram
#' @return heatmap with the results of cor.assoc
#' @export
#'
#' @examples
#' library(stats)
#' 
#' #load the dataset
#' m <- as.matrix(synthetic_metabolic_dataset)
#' 
#' #Compute the pearson correlation of all the variables in the data.frame metabolic_measures
#' cors<-cor_assoc(m, m, MiMIR::metabolites_subsets$MET63,MiMIR::metabolites_subsets$MET63)
#' #Plot the correlations
#' plot_corply(cors, main="Correlations metabolites")
#' 
#' @seealso 
#' cor_assoc
#' 
#' @export
#' 
plot_corply <- function(res,main=NULL,zlim=NULL,reorder.x=FALSE,reorder.y=reorder.x,resort_on_p=FALSE,
                           abs=FALSE,cor.abs=FALSE,reorder_dend=FALSE){

  if(reorder.x|reorder.y){
    if (resort_on_p){
      IND <- resort.on.p(res,abs=abs)
    }else{
      IND <- resort.on.s(res,abs=abs)
    }
  }
  s <- get.s(res)
  p <- get.p(res)
  if(reorder.x){
    s <- s[IND$y,]
    p <- p[IND$y,]
  }
  if(reorder.y){
    s <- s[,IND$x]
    p <- p[,IND$x]
  }
  if(is.null(main)){
    main <- "Phenotypes vs. Features"
  }
  # Plot:
  if(cor.abs){s<-abs(s)}
  if(reorder_dend){
    heatmaply::heatmaply_cor(s,main=paste("<b>",main,"<b>"), margins=c(0,0,40,0))
  }else{
    heatmaply::heatmaply_cor(s,Rowv=F, Colv=F, main=paste("<b>",main,"<b>"),margins=c(0,0,40,0))
  }
  
}


#' plattCalibration
#' 
#' Function that calculates the Platt Calibrations
#'
#' @param r.calib observed binary phenotype
#' @param p.calib predicted probabilities
#' @param nbins number of bins to create the plots
#' @param pl logical indicating if the function should plot the Reliability diagram and histogram of the calibrations
#' @return list with samples, responses, calibrations, ECE, MCE and calibration plots if save==T
#'
#' @examples
#' library(stats)
#' library(plotly)
#' 
#' #load the dataset
#' met <- synthetic_metabolic_dataset
#' phen <- synthetic_phenotypic_dataset
#' 
#' #Calculating the binarized surrogates
#' b_phen<-binarize_all_pheno(phen)
#' #Apply a surrogate models and plot the ROC curve
#' surr<-calculate_surrogate_scores(met, phen,MiMIR::PARAM_surrogates, bin_names=colnames(b_phen))
#' #Calibration of the surrogate sex
#' real_data<-as.numeric(b_phen$sex)
#' pred_data<-surr$surrogates[,"s_sex"]
#' plattCalibration(r.calib=real_data, p.calib=pred_data, nbins = 10, pl=TRUE)
#' 
#' @details 
#' Many popular machine learning algorithms produce inaccurate predicted probabilities, especially when applied on a dataset different than the training set.
#' Platt (1999) proposed an adjustment, in which the original probabilities are used as a predictor in a single-variable logistic regression to produce more accurate adjusted predicted probabilities.
#' The function will also help the evaluation of the calibration, by plotting: reliability diagrams and distributions of the calibrated and non-calibrated probabilities.
#' The reliability diagrams plots the mean predicted value within a certain range of posterior probabilities, against the fraction of accurately predicted values.
#' Finally, we also report accuracy measures for the calibrations: the ECE, MCE and the Log-Loss of the probabilities before and after calibration.
#' 
#' @references
#' This is a function originally created for the package in eRic, under the name prCalibrate and modified ad hoc for our purposes
#' (\href{https://rdrr.io/github/etlundquist/eRic/man/prCalibrate.html}{Github})
#' 
#' J. C. Platt, 'Probabilistic Outputs for Support Vector Machines and Comparisons to Regularized Likelihood Methods', in Advances in Large Margin Classifiers, 1999, pp. 61-74.
#' 
#' 
#' @export
#' 
plattCalibration<- function (r.calib, p.calib, nbins = 10, pl=FALSE) {

  pred <- p.calib
  resp <- r.calib
  sample <- "calibration"
  
  pred[pred == 0] <- 1e-08
  pred[pred == 1] <- 1 - 1e-08
  cmodel <- stats::glm(y ~ x, data.frame(y = r.calib, x = p.calib), 
                family = "binomial")
  calibrated <- stats::predict(cmodel, data.frame(y = resp, x = pred), 
                        type = "response")
  
  raw.bins <- cut(pred, nbins, include.lowest = TRUE)
  raw.xval <- tapply(pred, raw.bins, mean)
  raw.yval <- tapply(resp, raw.bins, mean)
  raw.cali <- data.frame(method = rep("Original", nbins), 
                         x = raw.xval, y = raw.yval)
  raw.logl <- suppressWarnings((-1/length(resp)) * sum(resp * log(pred) + 
                                        (1 - resp) * (log(1 - pred)), na.rm = TRUE))
  cal.bins <- cut(calibrated, nbins, include.lowest = TRUE)
  cal.xval <- tapply(calibrated, cal.bins, mean)
  cal.yval <- tapply(resp, cal.bins, mean)
  cal.cali <- data.frame(method = rep("Calibrated", nbins), 
                         x = cal.xval, y = cal.yval)
  cal.logl <- (-1/length(resp)) * sum(resp * log(calibrated) + (1 - resp) * (log(1 - calibrated)), na.rm = TRUE)
 
  ## Calibration Errors ##
  ece_orig <- getECE(resp, pred, nbins)
  ece_calib <- getECE(resp, calibrated, nbins)
  ECE<-data.frame(ece_orig, ece_calib)
  mce_orig <- getMCE(resp, pred, nbins)
  mce_calib <- getMCE(resp, calibrated, nbins)
  MCE<-data.frame(mce_orig, mce_calib)
  
  if(pl){
    ## Plots ##
    #Layouts
    axis_font <- list(
      family = "Arial",
      size = 14
    )
    title_font <- list(
      family = "Arial",
      size = 18,
      margin=10
    )
    
    cdata<-cbind(raw.cali[,c(2,3)], cal.cali[,c(2,3)])
    colnames(cdata)<-c("x.original","y.original","x.calibrated","y.calibrated")
    
    cPlot <- plotly::plot_ly(cdata, x= ~x.original, y= ~y.original, name="Original", 
                     type = 'scatter', mode = 'lines+markers', line = list(width = 2), legendgroup = ~x.original) %>% 
      plotly::add_trace(x= ~x.calibrated, y= ~y.calibrated, name = 'Calibrated', 
                mode = 'lines+markers', color= I("red"), line = list(width = 2)) %>% 
      plotly::add_segments(x = 0, xend = 1, y = 0, yend = 1, name="Perfect", line = list(dash = "dash", width= 1, color="black")) %>%
      plotly::layout(title = list(text='<b>Reliability Diagram<b>',font=title_font, y = 0.98),
             xaxis = list(title = '<b>Confidence per bin<b>',titlefont = axis_font),
             yaxis = list (title = '<b>Fraction of positives<b>',titlefont = axis_font),
             margin = list(l = 50, r = 140, b = 30, t = 30, pad = 4),
             annotations = list(
               text = list(paste("ECE Original=",round(ece_orig, digits = 3),
                                 "\nMCE Original=",round(mce_orig, digits = 3)),
                           paste("ECE Calibrated=",round(ece_calib, digits = 3),
                                 "\nMCE Calibrated=",round(mce_calib, digits = 3))),
               xref='paper',
               yref='paper',
               x = c(1.2, 1.2), y = c(0, -0.2),
               showarrow=FALSE,
               bordercolor=c('blue','red'),
               borderwidth=2))
    
    hdata<-data.frame(pred, calibrated)
    hPlot<-plotly::plot_ly(hdata, x = ~pred, opacity = 0.45, type = "histogram", 
                   marker = list(color = "blue"), name = "Original",legendgroup = ~pred) %>% 
      plotly::add_histogram(x = calibrated, marker = list(color = "red"), 
                    name = "Calibrated") %>% 
      plotly::layout( title =list(text=paste("<b> Distributions <b>"),font=title_font,y = 0.98),
              xaxis = list(
                title = paste("<b> Predicted Probability <b>"),
                titlefont = axis_font), 
              yaxis = list(
                title = "<b> Counts <b>",
                titlefont = axis_font),barmode = "overlay")
    
    return(list(sample = sample,responses = resp, raw.probs = pred, cal.probs = calibrated, 
                raw.logloss = raw.logl, cal.logloss = cal.logl, 
                cal.model = cmodel, cal.Plot = cPlot, prob.hist = hPlot, ECE=ECE, MCE=MCE))
  }else{
    return(list(sample = sample, responses = resp, raw.probs = pred, cal.probs = calibrated, 
                 raw.logloss = raw.logl, cal.logloss = cal.logl, 
                cal.model = cmodel, ECE=ECE, MCE=MCE))
  }
}


#' MetaboWAS
#' 
#' Function to calculate a Metabolome Wide Association study
#'
#' @param met numeric data.frame with the metabolomics features
#' @param pheno data.frame containing the phenotype of interest
#' @param test_variable string vector with the name of the phenotype of interest
#' @param covariates string vector with the name of the variables to be added as a covariate
#' @param img logical indicating if the function should plot a Manhattan plot
#' @param adj_method multiple testing correction method
#' @return res= the results of the MetaboWAS, manhplot= the Manhattan plot made with plotly, N_hits= the number of significant hits
#' 
#'
#' @examples
#' require(MiMIR)
#' require(plotly)
#' require(ggplot2)
#' 
#' #' #load the dataset
#' metabolic_measures <- synthetic_metabolic_dataset
#' phenotypes <- synthetic_phenotypic_dataset
#' 
#' #Computing a MetaboWAS for age corrected by sex
#' MetaboWAS(met=metabolic_measures, pheno=phenotypes, test_variable="age", covariates= "sex")
#' 
#' @details 
#' This is a function to compute linear associations individually for each variable in the first data.frame
#' with the test variable and corrected for the selected covariates. 
#' This function to computes linear regression modelindividually for each variable in the first data.frame
#' with the test variable and adjusted for potential confounders. 
#' False Discovery Rate (FDR) is applied to account for multiple testing correction. 
#' The user has the faculty to select the test variable and the potential covariates within the pool of variables in the phenotypic file input.
#' The results of the associations are reported in a Manhattan plot 
#' 
#' The p-value of the association is then corrected using Benjamini Hochberg.
#' Finally we use plotly to plot a Manhattan Plot, which reports on the x-axis the list of metabolites reported in the Nightingale Health,
#' divided in groups, and on the y-axis the -log (adjusted p-value).
#' 
#' @references
#' This method is also described and used in:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, <doi:10.1016/j.ebiom.2021.103764>
#' 
#' @export
#'
MetaboWAS<-function(met, pheno, test_variable, covariates, img=TRUE, adj_method="BH"){
  res<-do_metabowas(phen=pheno, dat= met, test_variable = test_variable, covariates=covariates, adj_method=adj_method)
  N_hits<-length(which(res$pval.adj<=0.05))
  
  if(img){
    metabo_names_translator<-MiMIR::metabo_names_translator[which(!is.na(MiMIR::metabo_names_translator$BBMRI_names)),]
    res<-res[MiMIR::metabo_names_translator$BBMRI_names,]
    res$groups<-MiMIR::metabo_names_translator$group
    res$ord<-1:dim(res)[1]
    res$met_name<-rownames(res)
    
    axis_set <- res %>% 
      group_by(groups) %>% 
      summarise(center = mean(ord))
    
    ylim <- res %>% 
      filter(pval.adj == min(pval.adj,na.rm = T)) %>% 
      mutate(ylim = abs(floor(log10(pval.adj))) + 2)
    ylim<-ylim$ylim
    
    if(!is.null(covariates)){
      if(length(covariates)>1){
        title<-paste("<b>MetaboWAS:",test_variable,"corrected for",paste(covariates, collapse=","),"<b>",
                     "\nNumber of significant hits:",N_hits)
      }else{
        title<-paste("MetaboWAS:",test_variable,"corrected for",covariates,
                     "\nNumber of significant hits:",N_hits)
      }
      
    }else{
      title<-paste("<b>MetaboWAS:",test_variable,"<b>",
                   "\nNumber of significant hits:",N_hits)
    }
    
    manhplot <- ggplot2:: ggplot(res, ggplot2::aes(label = met_name, label2=groups, label3=pval.adj)) +
      ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
      ggplot2::geom_point(ggplot2::aes(x = ord, y = -log10(pval.adj),
                     color = as.factor(groups), size = -log10(pval.adj)), alpha = 0.75) +
      ggplot2::ggtitle(title)+
      ggplot2::scale_x_continuous(label = axis_set$groups, breaks = axis_set$center) +
      ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      ggplot2::scale_size_continuous(range = c(0.5,3)) +
      ggplot2::labs(x = "metabolites", 
           y = "-log(p adjusted)") + 
      ggplot2::theme_minimal() +
      ggplot2::theme( 
        legend.position = "none",
        panel.border = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.x = ggplot2:: element_text(angle = 60, size = 8, vjust = 0.5)
      )
    
    manhplot<-plotly::ggplotly(manhplot, 
                       tooltip=c("met_name", "groups", "pval.adj"))
  }
  
  list(res=res, manhplot=manhplot, N_hits=N_hits)
}










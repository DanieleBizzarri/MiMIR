#' subset_metabolites_overlap
#' 
#' Helper function that subsets the NH-metabolomics matrix features to the selection in metabolites needed for the metabolic score
#'
#' @param x numeric data-frame with Nightingale-metabolomics
#' @param metabos vector of strings containing the names of the metabolic features to be selected
#' @param quiet logical to suppress the messages in the console
#' @return matrix with the selected Nightingale-metabolomics features 
#'
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' #Select the metabolic features
#' mat <- subset_metabolites_overlap(x=metabolic_measures,metabos=PARAM_metaboAge$MET)
#' }
#' 
#' @references 
#' This function is constructed to be able to apply the metaboAge as described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541–547, doi:10.1161/CIRCGEN.119.002610
#' 
#' @seealso 
#' QCprep, apply.fit, subset_samples_miss, subset_samples_zero, subset_samples_sd, impute.miss, apply.scale, and report.dim
#'  
#' @keywords internal
#' 
#' 
subset_metabolites_overlap<-function(x,metabos,quiet=FALSE){
  x <- x[,intersect(colnames(x),metabos),drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header="Selecting metabolites"))
  }
  return(invisible(x))
}

#' subset_samples_miss
#' 
#' Helper function that subsets the NH-metabolomics matrix to the samples with less than Nmax missing values
#'
#' @param x numeric data-frame with Nightingale-metabolomics
#' @param Nmax integer indicating  the max number of missing values allowed per sample (N suggested= 1)
#' @param quiet logical to suppress the messages in the console
#' @return matrix with the samples with limited amount of missing values in the Nightingale-metabolomics dataset
#'
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' #Select the samples with only 1 missing value
#' mat <- subset_samples_miss(x=metabolic_measures, Nmax=1)
#' 
#' }
#' 
#' @references 
#' This function is constructed to be able to apply the metaboAge as described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541–547, doi:10.1161/CIRCGEN.119.002610
#' 
#' @seealso 
#' QCprep, apply.fit, subset_metabolites_overlap, subset_samples_zero, subset_samples_sd, impute.miss, apply.scale, and report.dim
#'  
#' @keywords internal
#' 
subset_samples_miss<-function(x,Nmax=1,quiet=FALSE){
  MISS <- colSums(is.na(t(x)))
  x <- x[which(MISS<=Nmax),,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on missing values [Nmax>=",Nmax,"]")))
  }
  return(invisible(x))
}

#' subset_samples_miss
#' 
#' Helper function that subsets the NH-metabolomics matrix to the samples with less than Nmax zeros

#' @param x numeric data-frame with Nightingale-metabolomics
#' @param Nmax integer indicating  the max number of missing values allowed per sample (N suggested= 1)
#' @param quiet logical to suppress the messages in the console
#' @return matrix with the samples with limited amount of zeros in the Nightingale-metabolomics dataset
#'
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' #Select samples with only 1 zero
#' mat <- subset_samples_zero(x=metabolic_measures, Nmax=1)
#' 
#' }
#' 
#' @references 
#' This function is constructed to be able to apply the metaboAge as described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541–547, doi:10.1161/CIRCGEN.119.002610
#' 
#' @seealso 
#' QCprep, apply.fit, subset_metabolites_overlap, subset_samples_miss, subset_samples_sd, impute.miss, apply.scale, and report.dim
#' 
#' @keywords internal
#' 
subset_samples_zero<-function(x,Nmax=1,quiet=FALSE){
  ZERO <- colSums(t(x==0),na.rm=TRUE)
  x <- x[which(ZERO<=Nmax),,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on zero values [Nmax>=",Nmax,"]")))
  }
  return(invisible(x))
}

#' subset_samples_sd
#' 
#' Helper function that subsets the NH-metabolomics matrix to the samples with limited numbers of outliers
#'
#' @param x numeric data-frame with Nightingale-metabolomics
#' @param MEAN numeric vector indicating the mean of the metabolites in x
#' @param SD numeric vector indicating the standard deviations of the metabolites in x
#' @param quiet logical to suppress the messages in the console
#' @return matrix with the samples with limited amount of outliers in the Nightingale-metabolomics dataset
#' 
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' #Select the samples with low outliers
#' mat <- subset_samples_sd(x=metabolic_measures, Nmax=1)
#' }
#' 
#' @references 
#' This function is constructed to be able to apply the metaboAge as described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541–547, doi:10.1161/CIRCGEN.119.002610
#' 
#' @seealso 
#' QCprep, apply.fit, subset_metabolites_overlap, subset_samples_miss, subset_samples_zero, impute.miss, apply.scale, and report.dim
#' 
#' @keywords internal 
#'
subset_samples_sd<-function(x,MEAN,SD,quiet=FALSE){
  MEAN <- MEAN[colnames(x)]
  SD <- SD[colnames(x)]
  Dummi <- log(x)
  Dummi[which(Dummi==-Inf)] <- NA
  # Exclude persons being an outlier:
  outl_samp <- rownames(Dummi)[unique(which(((Dummi > t(replicate(nrow(Dummi),MEAN)) + 5*t(replicate(nrow(Dummi),SD))) | (Dummi < t(replicate(nrow(Dummi),MEAN)) - 5*t(replicate(nrow(Dummi),SD)))),arr.ind=TRUE)[,"row"])]
  if(!is.null(outl_samp)){
    sample_names <- setdiff(rownames(Dummi),outl_samp)
    x <- x[sample_names,,drop=FALSE]
  }
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on 5SD")))
  }
  return(invisible(x))
}

#' impute.miss
#' 
#' Helper function that subsets the NH-metabolomics matrix to the samples with less than Nmax zeros
#' 
#' Function created that subsets the NH-metabolomics matrix samples to the ones for which the metabolites
#' included in MetaboAge for which the log of the metabolic concentrations are not more than 5SD away from their mean
#'
#' @param x numeric data-frame with Nightingale-metabolomics
#' @return matrix of the Nightingale-metabolomics dataset with missing values imputed to zero
#' 
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' #Imputing missing values
#' mat <- impute.miss(metabolic_measures)
#' }
#' 
#' @references 
#' This function is constructed to be able to apply the metaboAge as described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541–547, doi:10.1161/CIRCGEN.119.002610
#' 
#' @seealso 
#' QCprep, apply.fit, subset_metabolites_overlap, subset_samples_miss, subset_samples_zero, subset_samples_sd, apply.scale, and report.dim
#' 
#' @keywords internal
#'
impute.miss<-function(x){
  ## This is an boiler-plate solution :)
  x[which(is.na(x))] <- 0
  return(x)
}


#' apply.scale
#' 
#' Helper function created to scale the NH-metabolomics matrix samples
#'
#' @param dat numeric data-frame with Nightingale-metabolomics
#' @param MEAN numeric vector indicating the mean of the metabolites present in dat
#' @param SD numeric vector indicating the standard deviations of the metabolites present in dat
#' @param quiet Tlogical to suppress the messages in the console
#' @return The matrix z-scaling the Nightingale-metabolomics dataset using the given Means and SDs
#' @export
#'
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' #Apply the scaling to the metabolic features
#' mat <- apply.scale(metabolic_measures, MEAN=PARAM_metaboAge$MEAN, SD=PARAM_metaboAge$SD)
#' }
#' 
#' @references 
#' This function is constructed to be able to apply the metaboAge as described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541–547, doi:10.1161/CIRCGEN.119.002610
#' 
#' @seealso 
#' QCprep, apply.fit, subset_metabolites_overlap, subset_samples_miss, subset_samples_zero, subset_samples_sd, impute.miss, and report.dim
#' 
#' @keywords internal
#'
apply.scale <- function(dat,MEAN,SD,quiet=TRUE){
  if(!quiet){
    cat("| Apply scaling ... ")
  }
  COLNAMES <- colnames(dat)
  #names(SD)<-COLNAMES
  dat <- sapply(COLNAMES,function(x) (dat[,x]-MEAN[x])/SD[x])
  if(!quiet){
    cat("Done! \n")
  }
  return(dat)
}

#' report.dim
#' 
#' Helper function to report on the console the dimension of the NH metabolomics matrix
#'
#' @param x numeric data-frame with Nightingale-metabolomics
#' @param header string describing the sub-sampling of the NH-metabolomics matrix
#' @param trailing number of digits to show
#' @return The report of the NH-metabolomics matrix dimension
#'
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' #Apply the scaling to the metabolic features
#' cat(report.dim(x, header=paste0("Pruning samples on 5SD")))
#' }
#' 
#' @references 
#' This function is constructed to be able to apply the metaboAge as described in:
#' van den Akker Erik B. et al. (2020) Metabolic Age Based on the BBMRI-NL 1H-NMR Metabolomics Repository as Biomarker of Age-related Disease. Circulation: Genomic and Precision Medicine, 13, 541–547, doi:10.1161/CIRCGEN.119.002610
#' 
#' @seealso 
#' QCprep, apply.fit, subset_metabolites_overlap, subset_samples_miss, subset_samples_zero, subset_samples_sd, impute.miss, and apply.scale
#' 
#' @keywords internal
#'
report.dim<-function(x,header,trailing="0"){
  return(paste0(sprintf(paste0("%-",trailing,"s"),paste0("| ",header,":\n")),sprintf("%4s",ncol(x))," metabolites x ",sprintf("%4s",nrow(x))," samples \n"))
}

#' subset_samples_sd_surrogates
#' 
#' Helper function that subsets the NH-metabolomics matrix to the samples with limited numbers of outliers
#'
#' @param x numeric data-frame with Nightingale-metabolomics
#' @param MEAN numeric vector indicating the mean of the metabolites in x
#' @param SD numeric vector indicating the standard deviations of the metabolites in x
#' @param N numeric vector indicating the amount of standard deviations away from the mean after which we consider an outlier (N suggested=5)
#' @param quiet logical to suppress the messages in the console
#' @return matrix with the samples with limited amount of outliers in the Nightingale-metabolomics dataset
#' 
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' #Select the samples with low outliers
#' mat <- subset_samples_sd_surrogates(x=metabolic_measures, Nmax=1)
#' }
#' 
#' @details 
#' Bizzarri et al. built multivariate models,using 56 metabolic features quantified by Nightingale, to predict the 19 binary characteristics of an individual. 
#' The binary variables are: sex, diabetes status, metabolic syndrome status, lipid medication usage, blood pressure lowering medication,
#' current smoking, alcohol consumption, high age, middle age, low age, high hsCRP, high triglycerides, high ldl cholesterol,
#' high total cholesterol, low hdl cholesterol, low eGFR, low white blood cells, low hemoglobin levels.
#' 
#' @references 
#' This function was made to vidualize the binarized variables calculated following the rules indicated in the article:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, doi: 10.1016/j.ebiom.2021.103764
#' 
#' @seealso 
#' QCprep_surrogates, calculate_surrogate_scores, apply.fit_surro
#' 
#' @keywords internal
#'
subset_samples_sd_surrogates<-function(x,MEAN,SD, N=5, quiet=FALSE){
  MEAN <- MEAN[colnames(x)]
  SD <- SD[colnames(x)]
  Dummi <- x
  # Exclude persons being an outlier:
  outl_samp <- rownames(Dummi)[unique(which(((Dummi > t(replicate(nrow(Dummi),MEAN)) + N*t(replicate(nrow(Dummi),SD))) | (Dummi < t(replicate(nrow(Dummi),MEAN)) - 5*t(replicate(nrow(Dummi),SD)))),arr.ind=TRUE)[,"row"])]
  if(!is.null(outl_samp)){
    sample_names <- setdiff(rownames(Dummi),outl_samp)
    x <- x[sample_names,,drop=FALSE]
  }
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on", N, "SD")))
  }
  return(invisible(x))
}


#' apply.fit_surro
#' 
#' Function that apply on of the surrogates models to the NH-metabolomics concentrations
#'
#' @param mat numeric data-frame with Nightingale-metabolomics
#' @param FIT The betas of the logistic regressions composing the surrogates by Bizzarri et al.
#' @param post logical to obtain posterior probabilities
#' @return numeric data.frame with the metabolomics-based surrogates by Bizzarri et al.
#' 
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' # Do the pre-processing steps to the metabolic measures
#' metabolic_measures<-QCprep_surrogates(as.matrix(metabolic_measures), Nmax_miss=1,Nmax_zero=1)
#' 
#' #load the phenotypic dataset
#' phenotypes <- read.csv("phenotypes_file_path",header = TRUE, row.names = 1)
#' #Calculating the binarized surrogates
#' bin_pheno<-binarize_all_pheno(phenotypes)
#' 
#' #Apply the surrogate models
#' surrogates<-foreach::foreach(i=MiMIR::phenotypes_names$out_surro, .combine="cbind") %do% {
#' pred<-apply.fit_surro(as.matrix(metabo_measures), 
#' PARAM_surrogates$models_betas[i,])}
#' 
#' }
#' 
#' @details 
#' Bizzarri et al. built multivariate models,using 56 metabolic features quantified by Nightingale, to predict the 19 binary characteristics of an individual. 
#' The binary variables are: sex, diabetes status, metabolic syndrome status, lipid medication usage, blood pressure lowering medication,
#' current smoking, alcohol consumption, high age, middle age, low age, high hsCRP, high triglycerides, high ldl cholesterol,
#' high total cholesterol, low hdl cholesterol, low eGFR, low white blood cells, low hemoglobin levels.
#' 
#' @references 
#' This function was made to vidualize the binarized variables calculated following the rules indicated in the article:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, doi: 10.1016/j.ebiom.2021.103764
#' 
#' @seealso 
#' QCprep_surrogates, calculate_surrogate_scores, subset_samples_sd_surrogates, predictions_surrogates
#' 
#' @keywords internal
#'
apply.fit_surro<-function(mat, FIT, post=TRUE){
  # Resort:
  BETA <- FIT[colnames(mat)]
  INTC <- FIT[1]
  if(post){
    # Predict surrogate:
    p<-as.numeric(1/(1+exp(-1 * (INTC+ rowSums(mat %*% BETA)))))
  }else{
    p <- as.numeric(as.vector(mat %*% BETA) + as.vector(INTC),stringsAsFactors=FALSE)
  }
  
  names(p)<-rownames(mat)
  return(p)
}


#' predictions_surrogates
#' 
#' Helper function that apply a surrogate model and plot a ROC curve the accuracy
#'
#' @param FIT numeric vector with betas of the logistic regressions composing the surrogates by Bizzarri et al.
#' @param data numeric data-frame with Nightingale-metabolomics and the binarized phenotype to predict
#' @param title_img string with title of the image
#' @param plot logical to obtain the ROC curve
#' @return If plot==T The surrogate predictions and the roc curve. If plot==F only the surrogate predictions
#'
#' @examples
#' \dontrun{
#' library(MiMIR)
#' 
#' #load the Nightignale metabolomics dataset
#' metabolic_measures <- read.csv("Nightingale_file_path",header = TRUE, row.names = 1)
#' # Do the pre-processing steps to the metabolic measures
#' metabolic_measures<-QCprep_surrogates(as.matrix(metabolic_measures), Nmax_miss=1,Nmax_zero=1)
#' 
#' #load the phenotypic dataset
#' phenotypes <- read.csv("phenotypes_file_path",header = TRUE, row.names = 1)
#' #Calculating the binarized surrogates
#' bin_pheno<-binarize_all_pheno(phenotypes)
#' 
#' #Apply a surrogate models and plot the ROC curve
#' data<-data.frame(out=factor(phenotypes_names$bin_names[,1]), metabo_measures)
#' colnames(data)[1]<-"out"
#' pred<-predictions_surrogates(PARAM_surrogates$models_betas["s_sex",], data=data, title_img="s_sex")
#' 
#' }
#' 
#' @details 
#' Bizzarri et al. built multivariate models,using 56 metabolic features quantified by Nightingale, to predict the 19 binary characteristics of an individual. 
#' The binary variables are: sex, diabetes status, metabolic syndrome status, lipid medication usage, blood pressure lowering medication,
#' current smoking, alcohol consumption, high age, middle age, low age, high hsCRP, high triglycerides, high ldl cholesterol,
#' high total cholesterol, low hdl cholesterol, low eGFR, low white blood cells, low hemoglobin levels.
#' 
#' @references 
#' This function was made to vidualize the binarized variables calculated following the rules indicated in the article:
#' Bizzarri,D. et al. (2022) 1H-NMR metabolomics-based surrogates to impute common clinical risk factors and endpoints. EBioMedicine, 75, 103764, doi: 10.1016/j.ebiom.2021.103764
#' 
#' @seealso 
#' QCprep_surrogates, calculate_surrogate_scores, subset_samples_sd_surrogates, apply.fit_surro
#' 
#' @keywords internal
#'
predictions_surrogates<- function(FIT, data, title_img=FALSE, plot=TRUE){
  
  pred<-apply.fit_surro(as.matrix(data[,-1]),FIT)
  
  all(rownames(pred) == rownames(data))
  
  if(plot){
    #Evaluation of the predictions
    pred_eval<-pROC::roc(data[,"out"], pred, plot=T, col="#377eb8", quiet = T,
                         lwd=4, print.auc=TRUE, main = title_img, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
    
    res<-list(predictions=pred, predictions_eval=pred_eval)
  }else{
    res<-list(predictions=pred)
  }
  return(res)
}


#' model_coeff_heat
#' 
#' Function to plot the scaled coefficients of the metabolic scores
#'
#' @param mort_betas dataframe withthe coefficients of the mortality score
#' @param metaboAge_betas dataframe with the coefficients of the metaboAge
#' @param surrogates_betas dataframe with the coefficients of the surrogates
#' @param Ahola_Olli_betas dataframe with the coefficients of the T2D score
#' @param CVD_score_betas dataframe with the coefficients of the CVD score
#' @param COVID_score_betas ataframe with the coefficients of the COVID_score
#' @return heatmapply with the scaled coefficients of the metabolic scores
#' @keywords internal
#'
model_coeff_heat<-function(mort_betas, metaboAge_betas, surrogates_betas, Ahola_Olli_betas, CVD_score_betas, COVID_score_betas){
  models_betas<-matrix(0, 24, 62)
  rownames(models_betas)<-c("mortScore","MetaboAge", rownames(surrogates_betas),"T2D-score","CVD_score", "COVID_score")
  colnames(models_betas)<-MiMIR::metabolites_subsets$MET62
  models_betas["mortScore",mort_betas$Abbreviation]<-mort_betas$Beta_value
  models_betas["MetaboAge", names(metaboAge_betas)[-1]]<-metaboAge_betas[-1]
  #surrogates
  for(x in rownames(surrogates_betas)){
    models_betas[x,colnames(surrogates_betas)[-1]]<-surrogates_betas[x,-1]
  }
  met<-Ahola_Olli_betas$Abbreviation[which(Ahola_Olli_betas$Abbreviation %in% MiMIR::metabolites_subsets$MET62)]
  models_betas["T2D-score",met]<-Ahola_Olli_betas$Beta_value[which(Ahola_Olli_betas$Abbreviation %in% MiMIR::metabolites_subsets$MET62)]
  met<-CVD_score_betas$Abbreviation[which(CVD_score_betas$Abbreviation %in% MiMIR::metabolites_subsets$MET62)]
  models_betas["CVD_score",met]<-CVD_score_betas$Beta_value[which(CVD_score_betas$Abbreviation %in% MiMIR::metabolites_subsets$MET62)]
  met<-COVID_score_betas$Abbreviation[which(COVID_score_betas$Abbreviation %in% MiMIR::metabolites_subsets$MET62)]
  models_betas["COVID_score",met]<-COVID_score_betas$Beta_value[which(COVID_score_betas$Abbreviation %in% MiMIR::metabolites_subsets$MET62)]
  
  
  mat <-models_betas
  mat <- sapply(colnames(mat), function(colname) {
    sapply(rownames(mat), function(rowname) {
      paste("The original coefficient is", round(models_betas[rowname, colname],digits=3))
    })
  })
  
  heatmaply::heatmaply(models_betas,main=paste("<b> Models Coefficients scaled on the rows<b>"),
                       scale="row",custom_hovertext = mat,
                       scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", mid="white", high = "red", limits=c(-0.1146100, 0.14128)),
                       margins=c(0,0,40,0))
}



#' resort.on.s
#' 
#' helper function for clustering cor.assoc results reordering based on the associations
#'
#' @param res Results of cor.assoc
#' @param abs TRUE/FALSE. If TRUE it will cluster the absolute values
#' @return Returns the clustered order of the associations
#' @keywords internal
#'
resort.on.s <- function(res,abs=FALSE){
  s <- get.s(res)
  if(abs){
    s <- abs(s)
  }
  order.x <- stats::hclust(stats::as.dist(1-(stats::cor(s))))$order
  if(is.sym(res)){
    order.y <- order.x
  } else {
    order.y <- stats::hclust(stats::as.dist(1-(stats::cor(t(s)))))$order
  }
  return(list(x=order.x,y=order.y))
}

#' resort.on.p
#'
#' helper function for clustering cor.assoc results reordering based on the pvalues
#'
#' @param res Results of cor.assoc
#' @param abs TRUE/FALSE. If TRUE it will cluster the absolute values
#' @return Returns the clustered order of the associations based on the pvalues
#' 
#' @keywords internal
#'
resort.on.p <- function(res,abs=FALSE){
  p <- get.p(res)
  if(abs){
    p <- abs(p)
  }
  order.x <- stats::hclust(stats::as.dist(1-(stats::cor(log(p+1e-323)))))$order
  if(is.sym(res)){
    order.y <- order.x
  } else {
    order.y <- stats::hclust(stats::as.dist(1-(stats::cor(t(log(p+1e-323))))))$order
  }
  return(list(x=order.x,y=order.y))
}

#' get.s
#' 
#' helper function for extracting statistics from cor.assoc results
#'
#' @param res Results of cor.assoc
#' @return the matrix of associations
#' 
#' @keywords internal
#'
get.s <- function(res){
  return(
    sapply(res,function(x){
      if(is.data.frame(x)){
        # Data.frame coming from feat.expr
        x[,"logFC"]
      }else{
        # Matrix coming from cor.assoc
        if(is.matrix(x)){
          x[,"estimate.cor"]
        } else {
          sapply(x,function(y){
            # List coming from feat.assoc?
            if("Estimate" %in% colnames(y)){
              y["FEAT","Estimate"]
            } else {
              NA
            }
          })
        }
      }
    })
  )
}

#' get.p
#' 
#' helper function for extracting pvalues from cor.assoc results
#'
#' @param res Results of cor.assoc
#' @return the matrix of the pvalues of the associations
#' 
#' @keywords internal
#'
get.p <- function(res){
  return(
    sapply(res,function(x){
      if(is.data.frame(x)){
        # Data.frame coming from feat.expr
        x[,"PValue"]
      }else{
        # Matrix coming from cor.assoc
        if(is.matrix(x)){
          x[,"p.value"]
        }else{
          # List coming from feat.assoc
          sapply(x,function(y){
            if("Pr(>|t|)" %in% colnames(y)){
              y["FEAT","Pr(>|t|)"]
            } else {
              NA
            }
          })
        }
      }
    })
  )
}

#' is.sym
#'
#'Defining a helper function to check whether a supplied matrix is symmetric:
#' helper function to check whether a supplied matrix is symmetric
#'
#' @param res Results of cor.assoc
#' @return TRUE/FALSE
#' 
#' @keywords internal
#'
is.sym <- function(res){
  s <- get.s(res)
  if(nrow(s)==ncol(s)){
    if(all(colnames(s)==rownames(s))){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}


#' calibration_surro
#' 
#' helper function that calculates the Platt Calibrations  for all the surrogates
#'
#' @param bin_phenotypes data.frame with binary phenotypes resulting form binarize_all_pheno
#' @param surrogates data.frame with surrogates resulting from calculate_surrogate_scores
#' @param bin_names string with the names of the binarize clinical variables
#' @param bin_pheno_available string with the names of the binarize clinical variables available in the dataset
#' @param pl TRUE/FALSE. If TRUE creates the calibration plots
#' @param nbins number of bins for the plots
#' @return list with the calibrated surrogates
#' 
#' @seealso plattCalibration
#' 
#' @keywords internal
#'
calibration_surro<-function(bin_phenotypes, surrogates, bin_names, bin_pheno_available, pl=FALSE, nbins=10){
  bin_surro<-paste0("s_",bin_names)
  calib<-lapply(1:length(bin_surro), function(i){
    orig<-as.numeric(bin_phenotypes[,bin_names[i]])
    names(orig)<-rownames(bin_phenotypes)
    pred<-as.numeric(surrogates[,bin_surro[i]])
    names(pred)<-rownames(surrogates)
    ind<-which(!is.na(orig))
    if(length(ind)!=0){
      calibration<-plattCalibration(orig,pred, nbins = nbins, pl=pl)
      return(calibration)
    }
  })
  names(calib)<-bin_names
  calib_df<-calib_data_frame(calib, bin_phenotypes, bin_pheno_available)
  
  calib<-append(calib,list(calib_df))
  names(calib)<-c(bin_names,"calib_df")
  
  return(calib)
}


#' calib_data_frame
#' 
#' helper function that creates a data.frame with the Platt Calibrations
#'
#'@param calibrations list result of calibration_surro
#' @param bin_phenotypes data.frame with binary phenotypes resulting form binarize_all_pheno
#' @param bin_pheno_available string with the names of the binarize clinical variables available in the dataset
#' @return data.frame with the calibrated surrogates
#' 
#' @keywords internal
#'
calib_data_frame<-function(calibrations, bin_phenotypes, bin_pheno_available){
  cal<-matrix(NA, nrow = dim(bin_phenotypes)[1], ncol = length(bin_pheno_available))
  rownames(cal)<-rownames(bin_phenotypes)
  colnames(cal)<-bin_pheno_available
  for(i in 1:length(bin_pheno_available)){
    #cal[,i]<-calibrations[[bin_pheno_available[i]]]$cal.probs
    cal[names(calibrations[[bin_pheno_available[i]]]$cal.probs),i]<-calibrations[[bin_pheno_available[i]]]$cal.probs
  }
  colnames(cal)<-bin_pheno_available
  return(as.data.frame(cal))
}


#' Function that plots the Platt Calibrations using plotly
#'
#' @param r binary real data
#' @param p predicted probabilities
#' @param p.orig the uncalibrated posterior probabilities
#' @param name character string indicating the name of the calibrated variable
#' @param nbins number of bins to create the plots
#' @param annot_x integer indicating the x axis points in which the ECE and MCE values will be plotted
#' @param annot_y integer indicating the y axis points in which the ECE and MCE values will be plotted
#' @return list with Reliability diagram and histogram with calibrations and original predictions
#' 
#' @keywords internal
#'
plattCalib_evaluation<-function(r, p, p.orig, name, nbins = 10, annot_x=c(1,1),annot_y=c(0.1,0.3)){
  
  i<-which(is.na(r))
  if(length(i!=0)){
    r<-r[-i]
    p<-p[-i]
    p.orig<-p.orig[-i]
  }
  
  set.seed(222)
  train_ind<- caret::createDataPartition(r, p=0.8, list=FALSE)
  r.calib<-as.numeric(r[train_ind])
  p.calib<-p[train_ind]
  resp<-as.numeric(r[-train_ind])
  pred<-p[-train_ind]
  pred.orig<-p.orig[-train_ind]
  
  # add/subtract epsilon to predictions at (0,1) to avoid infinite log-loss
  
  pred[pred == 0] <- 1e-8
  pred[pred == 1] <- 1 - 1e-8
  
  # fit the calibration model and calculate calibrated probabilities
  
  cmodel <- stats::glm(y ~ x, data.frame(y = r.calib, x = p.calib), family = 'binomial')
  calibrated <- stats::predict(cmodel, data.frame(y = resp, x = pred), type = 'response')
  
  # calculate visualization/return measures for original probabilities (on either calibration or validation data) 
  
  raw.bins <- cut(pred.orig, nbins, include.lowest = TRUE)
  raw.xval <- tapply(pred.orig, raw.bins, mean)
  raw.yval <- tapply(resp, raw.bins, mean) 
  raw.cali <- data.frame(method = rep('Original', nbins), x = raw.xval, y = raw.yval)
  raw.logl <- (-1/length(resp)) * sum(resp*log(pred.orig) + (1-resp)*(log(1-pred.orig)), na.rm = TRUE)
  
  # calculate needed measures using transformed probabilities
  
  cal.bins <- cut(calibrated, nbins, include.lowest = TRUE)
  cal.xval <- tapply(calibrated, cal.bins, mean)
  cal.yval <- tapply(resp, cal.bins, mean)
  cal.cali <- data.frame(method = rep('Calibrated', nbins), x = cal.xval, y = cal.yval)
  cal.logl <- (-1/length(resp)) * sum(resp*log(calibrated) + (1-resp)*(log(1-calibrated)), na.rm = TRUE)
  
  ## Calibration Errors ##
  ece_orig <- getECE(resp, pred.orig, nbins)
  ece_calib <- getECE(resp, calibrated, nbins)
  ECE<-data.frame(ece_orig, ece_calib)
  mce_orig <- getMCE(resp, pred.orig, nbins)
  mce_calib <- getMCE(resp, calibrated, nbins)
  MCE<-data.frame(mce_orig, mce_calib)
  
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
                           type = 'scatter', mode = 'lines+markers', line = list(width = 2)) %>% 
    plotly::add_trace(x= ~x.calibrated, y= ~y.calibrated, name = 'Calibrated', 
                      mode = 'lines+markers', color= I("red"), line = list(width = 2)) %>% 
    plotly::add_segments(x = 0, xend = 1, y = 0, yend = 1, name="Perfect", line = list(dash = "dash", width= 1, color="black")) %>%
    plotly::layout(title = list(text=paste('<b>Reliability Diagram<b>',name,'<b>'),font=title_font, y = 0.98),
                   xaxis = list(title = '<b>Confidence per bin<b>',titlefont = axis_font),
                   yaxis = list (title = '<b>Fraction of positives<b>',titlefont = axis_font),
                   margin = list(l = 50, r = 140, b = 30, t = 30, pad = 4),
                   annotations = list(
                     text = list(paste("ECE Orig=",round(ece_orig, digits = 3),
                                       "\nMCE Orig=",round(mce_orig, digits = 3),
                                       "\nLogLoss Orig=",round(raw.logl, digits = 3)),
                                 paste("ECE Calib=",round(ece_calib, digits = 3),
                                       "\nMCE Calib=",round(mce_calib, digits = 3),
                                       "\nLogLoss Calib=",round(cal.logl, digits = 3))),
                     xref='paper',
                     yref='paper',
                     x =annot_x, y = annot_y,
                     showarrow=FALSE,
                     bordercolor=c('blue','red'),
                     borderwidth=2))
  
  hdata<-data.frame(pred.orig, calibrated)
  hPlot<-plotly::plot_ly(hdata, x = ~pred.orig, opacity = 0.45, type = "histogram", nbinsx = nbins,
                         marker = list(color = "blue"), name = "Original") %>% 
    plotly::add_histogram(x = calibrated, marker = list(color = "red"), 
                          name = "Calibrated") %>% 
    plotly::layout( title =list(text=paste("<b> Distributions <b>",name,'<b>'),font=title_font,y = 0.98),
                    xaxis = list(
                      title = paste("<b> Predicted Probability <b>"),
                      titlefont = axis_font), 
                    yaxis = list(
                      title = "<b> Counts <b>",
                      titlefont = axis_font),barmode = "overlay")
  
  return(list(train_ind= train_ind, responses = resp, raw.probs = pred, cal.probs = calibrated, 
              raw.logloss = raw.logl, cal.logloss = cal.logl, 
              cal.model = cmodel, cal.Plot = cPlot, prob.hist = hPlot))
}

#' getECE
#'
#' helper function to calculate the ECE of calibrations
#'
#' @param actual observed binary phenotype
#' @param predicted predicted values
#' @param n_bins the number of bins
#' @return ECE value
#' @keywords internal
#'
getECE<-function (actual, predicted, n_bins = 10){
  predicted <- predicted
  labels <- actual
  idx <- order(predicted)
  pred_actual <- (cbind(predicted[idx], labels[idx]))
  N <- nrow(pred_actual)
  rest <- N%%n_bins
  S <- 0
  W <- c()
  B <- min(N, n_bins)
  groups <- list()
  for (i in 1:B) {
    if (i <= rest) {
      group_pred <- (pred_actual[(((i - 1) * ceiling(N/n_bins) + 
                                     1):(i * ceiling(N/n_bins))), 1])
      group_actual <- (pred_actual[(((i - 1) * ceiling(N/n_bins) + 
                                       1):(i * ceiling(N/n_bins))), 2])
    }
    else {
      group_pred <- (pred_actual[((rest + (i - 1) * floor(N/n_bins) + 
                                     1):(rest + i * floor(N/n_bins))), 1])
      group_actual <- (pred_actual[((rest + (i - 1) * 
                                       floor(N/n_bins) + 1):(rest + i * floor(N/n_bins))), 
                                   2])
    }
    n_ <- length(group_pred)
    expected <- mean(group_pred, na.rm=TRUE)
    observed <- mean(group_actual, na.rm=TRUE)
    S[i] <- abs(observed - expected)
    W[i] <- n_/N
    groups[[i]] <- group_pred
  }
  mean_prediction <- lapply(groups, mean, na.rm=TRUE)
  min_group <- lapply(groups, min, na.rm=TRUE)
  max_group <- lapply(groups, max, na.rm=TRUE)
  res <- t(S) %*% W
  return(as.numeric(res))
}

#' helper function to calculate the MCE of the calibrations
#'
#' @param actual real values of the variables
#' @param predicted predicted values by one of the surrogates
#' @param n_bins the number of bins
#' @return MCE value
#' @keywords internal
#'
getMCE<-function (actual, predicted, n_bins = 10) 
{
  predicted <- predicted
  labels <- actual
  idx <- order(predicted)
  pred_actual <- (cbind(predicted[idx], actual[idx]))
  N <- nrow(pred_actual)
  rest <- N%%n_bins
  B <- min(N, n_bins)
  S <- 0
  W <- c()
  for (i in 1:B) {
    if (i <= rest) {
      group_pred <- (pred_actual[(((i - 1) * ceiling(N/n_bins) + 
                                     1):(i * ceiling(N/n_bins))), 1])
      group_actual <- (pred_actual[(((i - 1) * ceiling(N/n_bins) + 
                                       1):(i * ceiling(N/n_bins))), 2])
    }
    else {
      group_pred <- (pred_actual[((rest + (i - 1) * floor(N/n_bins) + 
                                     1):(rest + i * floor(N/n_bins))), 1])
      group_actual <- (pred_actual[((rest + (i - 1) * 
                                       floor(N/n_bins) + 1):(rest + i * floor(N/n_bins))), 
                                   2])
    }
    n <- length(group_pred)
    expected <- mean(group_pred, na.rm=TRUE)
    observed <- mean(group_actual, na.rm=TRUE)
    S[i] <- abs(observed - expected)
    W[i] <- n/N
  }
  res <- max(S * W)
  return(res)
}



#' plotly_NA_message
#' 
#' helper function to create a plotly indicating a problem with a plotly image
#'
#' @param main Message to plot
#' @return plotly image
#' 
#' @keywords internal
#'
plotly_NA_message<-function(main="Phenotype not available!"){
  title_font <- list(
    family = "Arial",
    size = 20,
    margin=10
  )
  
  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  
  plotly::plot_ly(x=0,y=0, opacity=0, type = "histogram")%>% 
    plotly::layout(title = list(text=paste('<b>',main,'<b>'),font=title_font, y = 0.5),
                   xaxis = ax, yaxis = ax)
}


#' Helper function to compute MetaboWASs
#' 
#' This helper function is called when doing Metabolites wide association analysis. It reports the results of linear regression models to study the association of a test variable to each metabolites individually and corrected for the covariates indicated.
#'
#' @param phen phenotypes data.frame
#' @param dat metabolites data.frame
#' @param test_variable the variable to be investigated 
#' @param covariates the covariates that you want to add
#' @param adj_method correction method.
#' @param quiet if FALSE it will plot the amount of people avaialble
#' @return results= the results of the MetaboWAS (estimate, tstatistics, pvalue, BH corrected pvalue)
#' 
#' @keywords internal
#'
do_metabowas<-function(phen, dat, test_variable="age", covariates=c("sex"), adj_method="BH", quiet=T){
  if(!is.null(covariates)){
    vars <- phen[, c(test_variable, covariates)]
  }else{
    vars <- as.data.frame(phen[,test_variable])
    rownames(vars)<-rownames(phen)
    colnames(vars)<-test_variable
  }
  vars <- stats:: na.omit(vars)
  if(!quiet){
    print(paste("The number of samples is:",dim(vars)[1]))
  }
  dat <- dat[match(rownames(vars), rownames(dat)),]
  allmisscols <- sapply(dat, function(x)all(is.na(x)))
  dat<-dat[,-which(allmisscols==T)]
  data<-merge(vars,dat,by=0)
  rownames(data)<-data$Row.names
  data<-data[,-1]
  fit<-foreach::foreach(i=colnames(dat), .combine="rbind") %do% {
    res <-data.frame(stats::coefficients(summary(stats::lm(stats::as.formula(paste0(test_variable,"~",paste(covariates,collapse="+"),"+",i)), data))))
    colnames(res) <- c("estimate","stdErr","tstat","pval")
    res<-res[,c("estimate","tstat","pval")]
    return(res[i,])
  }
  fit$pval.adj<-stats::p.adjust(fit$pval, method = adj_method)
  fit<-fit[order(fit$pval.adj, decreasing = F),]
  return(fit)
}


#' NA_message
#' 
#' helper function to create a plot indicating a problem with a plotly image
#'
#' @param main Message to plot
#' @return plot
#' 
#' @keywords internal
#'
NA_message<-function(main="Metabolites are missing, please check your upload!"){
  graphics::plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  graphics::text(x = 0.5, y = 0.5, main,cex = 1.6, col = "black", font=2)
}



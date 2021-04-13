#####################
### Libraries #######
#####################
require(shiny)
require(shinydashboard)
require(DT)
require(shinycssloaders)
require(glmnet)
require(glmnet)
require(pheatmap)
require(foreach)
require(RColorBrewer)
require(pROC)
require(plotly)
require(heatmaply)
require(matrixStats)
require(plyr)
###################################
### DEFINITIONS / CONSTANTS #######
###################################
# Metabolite markers may include all or some of:
METABO <- c(
  "xxl_vldl_p", "xxl_vldl_l", "xxl_vldl_pl", "xxl_vldl_c", "xxl_vldl_ce", "xxl_vldl_fc", "xxl_vldl_tg", 
  "xl_vldl_p",  "xl_vldl_l",  "xl_vldl_pl",  "xl_vldl_c",  "xl_vldl_ce",  "xl_vldl_fc",  "xl_vldl_tg", 
  "l_vldl_p",   "l_vldl_l",   "l_vldl_pl",   "l_vldl_c",   "l_vldl_ce",   "l_vldl_fc",   "l_vldl_tg", 
  "m_vldl_p",   "m_vldl_l",   "m_vldl_pl",   "m_vldl_c",   "m_vldl_ce",   "m_vldl_fc",   "m_vldl_tg", 
  "s_vldl_p",   "s_vldl_l",   "s_vldl_pl",   "s_vldl_c",   "s_vldl_ce",   "s_vldl_fc",   "s_vldl_tg", 
  "xs_vldl_p",  "xs_vldl_l",  "xs_vldl_pl",  "xs_vldl_c",  "xs_vldl_ce",  "xs_vldl_fc",  "xs_vldl_tg", 
  "idl_p",      "idl_l",      "idl_pl",      "idl_c",      "idl_ce",      "idl_fc",      "idl_tg", 
  "l_ldl_p",    "l_ldl_l",    "l_ldl_pl",    "l_ldl_c",    "l_ldl_ce",    "l_ldl_fc",    "l_ldl_tg", 
  "m_ldl_p",    "m_ldl_l",    "m_ldl_pl",    "m_ldl_c",    "m_ldl_ce",    "m_ldl_fc",    "m_ldl_tg", 
  "s_ldl_p",    "s_ldl_l",    "s_ldl_pl",    "s_ldl_c",    "s_ldl_ce",    "s_ldl_fc",    "s_ldl_tg", 
  "xl_hdl_p",   "xl_hdl_l",   "xl_hdl_pl",   "xl_hdl_c",   "xl_hdl_ce",   "xl_hdl_fc",   "xl_hdl_tg", 
  "l_hdl_p",    "l_hdl_l",    "l_hdl_pl",    "l_hdl_c",    "l_hdl_ce",    "l_hdl_fc",    "l_hdl_tg", 
  "m_hdl_p",    "m_hdl_l",    "m_hdl_pl",    "m_hdl_c",    "m_hdl_ce",    "m_hdl_fc",    "m_hdl_tg", 
  "s_hdl_p",    "s_hdl_l",    "s_hdl_pl",    "s_hdl_c",    "s_hdl_ce",    "s_hdl_fc",    "s_hdl_tg", 
  "xxl_vldl_pl_percentage", "xxl_vldl_c_percentage", "xxl_vldl_ce_percentage", "xxl_vldl_fc_percentage", "xxl_vldl_tg_percentage", 
  "xl_vldl_pl_percentage",  "xl_vldl_c_percentage",  "xl_vldl_ce_percentage",  "xl_vldl_fc_percentage",  "xl_vldl_tg_percentage", 
  "l_vldl_pl_percentage",   "l_vldl_c_percentage",   "l_vldl_ce_percentage",   "l_vldl_fc_percentage",   "l_vldl_tg_percentage", 
  "m_vldl_pl_percentage",   "m_vldl_c_percentage",   "m_vldl_ce_percentage",   "m_vldl_fc_percentage",   "m_vldl_tg_percentage", 
  "s_vldl_pl_percentage",   "s_vldl_c_percentage",   "s_vldl_ce_percentage",   "s_vldl_fc_percentage",   "s_vldl_tg_percentage", 
  "xs_vldl_pl_percentage",  "xs_vldl_c_percentage",  "xs_vldl_ce_percentage",  "xs_vldl_fc_percentage",  "xs_vldl_tg_percentage", 
  "idl_pl_percentage",      "idl_c_percentage",      "idl_ce_percentage",      "idl_fc_percentage",      "idl_tg_percentage", 
  "l_ldl_pl_percentage",    "l_ldl_c_percentage",    "l_ldl_ce_percentage",    "l_ldl_fc_percentage",    "l_ldl_tg_percentage", 
  "m_ldl_pl_percentage",    "m_ldl_c_percentage",    "m_ldl_ce_percentage",    "m_ldl_fc_percentage",    "m_ldl_tg_percentage", 
  "s_ldl_pl_percentage",    "s_ldl_c_percentage",    "s_ldl_ce_percentage",    "s_ldl_fc_percentage",    "s_ldl_tg_percentage", 
  "xl_hdl_pl_percentage",   "xl_hdl_c_percentage",   "xl_hdl_ce_percentage",   "xl_hdl_fc_percentage",   "xl_hdl_tg_percentage", 
  "l_hdl_pl_percentage",    "l_hdl_c_percentage",    "l_hdl_ce_percentage",    "l_hdl_fc_percentage",    "l_hdl_tg_percentage", 
  "m_hdl_pl_percentage",    "m_hdl_c_percentage",    "m_hdl_ce_percentage",    "m_hdl_fc_percentage",    "m_hdl_tg_percentage", 
  "s_hdl_pl_percentage",    "s_hdl_c_percentage",    "s_hdl_ce_percentage",    "s_hdl_fc_percentage",    "s_hdl_tg_percentage", 
  "vldl_d", "ldl_d", "hdl_d", "serum_c", "vldl_c", "remnant_c", "ldl_c", "hdl_c", "hdl2_c", "hdl3_c", "estc", "freec", 
  "serum_tg", "vldl_tg", "ldl_tg", "hdl_tg", "totpg", "tg_pg", "pc", "sm", "totcho", "apoa1", "apob", "apob_apoa1", "totfa", "unsatdeg", 
  "dha", "la", "faw3", "faw6", "pufa", "mufa", "sfa", "dha_fa", "la_fa", "faw3_fa", "faw6_fa", "pufa_fa", "mufa_fa", 
  "sfa_fa", "glc", "lac", "cit", "ala", "gln", "his", "ile", "leu", "val", "phe", "tyr", "ace", "acace", "bohbut", "crea", "alb", "gp",
  "pyr", "glycerol", "glycine")

# Metabolite markers often deemed to be more or less independent and not numerically derived from eachother:
MET63 <- tolower(c("Ala","Gln","His","Phe","Tyr","Ile","Leu","Val","Glc","Lac","Pyr","Cit","Ace",
                   "AcAce","bOHBut","Crea","Alb","Gp","XXL_VLDL_L","XL_VLDL_L","L_VLDL_L","M_VLDL_L",
                   "S_VLDL_L","XS_VLDL_L","IDL_L","L_LDL_L","M_LDL_L","S_LDL_L","XL_HDL_L","L_HDL_L",
                   "M_HDL_L","S_HDL_L","IDL_C","Serum_C","VLDL_C","LDL_C","HDL_C","HDL2_C","HDL3_C",
                   "VLDL_D","LDL_D","HDL_D","Serum_TG","TotPG","PC","SM","TotCho","ApoA1","ApoB",
                   "TotFA","DHA","LA","FAw3","FAw6","PUFA","MUFA","SFA","FAw3_FA","FAw6_FA","PUFA_FA",
                   "MUFA_FA","SFA_FA","UnsatDeg"))

# The final set of metabolites used by thesurrogates
MET56 <- tolower(c("Ala","Gln","His","Phe","Tyr","Ile","Leu","Val","Glc","Lac","Cit","Ace",
                   "AcAce","Crea","Alb","Gp","M_VLDL_L",
                   "S_VLDL_L","XS_VLDL_L","IDL_L","L_LDL_L","M_LDL_L","S_LDL_L",
                   "M_HDL_L","S_HDL_L","IDL_C","Serum_C","VLDL_C","LDL_C","HDL_C","HDL2_C","HDL3_C",
                   "VLDL_D","LDL_D","HDL_D","Serum_TG","TotPG","PC","SM","TotCho","ApoA1","ApoB",
                   "TotFA","DHA","LA","FAw3","FAw6","PUFA","MUFA","SFA","FAw3_FA","FAw6_FA","PUFA_FA",
                   "MUFA_FA","SFA_FA","UnsatDeg"))
# The final set of metabolites used by the models (metaboAge,surro and mortality score)
MET57 <- tolower(c("Ala","Gln","His","Phe","Tyr","Ile","Leu","Val","Glc","Lac","Cit","Ace",
                   "AcAce","Crea","Alb","Gp","M_VLDL_L",
                   "S_VLDL_L","XS_VLDL_L","IDL_L","L_LDL_L","M_LDL_L","S_LDL_L",
                   "M_HDL_L","S_HDL_L","IDL_C","Serum_C","VLDL_C","LDL_C","HDL_C","HDL2_C","HDL3_C",
                   "VLDL_D","LDL_D","HDL_D","Serum_TG","TotPG","PC","SM","TotCho","ApoA1","ApoB",
                   "TotFA","DHA","LA","FAw3","FAw6","PUFA","MUFA","SFA","FAw3_FA","FAw6_FA","PUFA_FA",
                   "MUFA_FA","SFA_FA","UnsatDeg","XXL_VLDL_L"))

MET14<-c("pufa_fa","gp","glc","s_hdl_l","xxl_vldl_l","alb","phe","acace","ile","vldl_d",
               "leu","val","his","lac")

# Mortality score betas
mort_betas <- data.frame(
  Abbreviation=c("pufa_fa","gp","glc","s_hdl_l","xxl_vldl_l","alb","phe","acace","ile","vldl_d",
                 "leu","val","his","lac"),
  Metabolite=c("Ratio of polyunsaturated fatty acids to total fatty acids (%)",
               "Glycoprotein acetyls",
               "Glucos",
               "Total lipids in small HD",
               "Total lipids in chylomicrons and extremely large VLDL",
               "Albumin",
               "Phenylalanine",
               "Acetoacetate",
               "Isoleucine",
               "Mean diameter for VLDL particles",
               "Leucine",
               "Valine",
               "Histidine",
               "Lactate"),
  Beta_value=c(-0.251264056,0.280179378,0.148392648,-0.141912413,-0.223685945,-0.113290967,
               0.124227231,0.078704926,0.205706859,-0.160137595,-0.195741112,-0.143698197,
               -0.069223326,0.062262328),
  stringsAsFactors = FALSE)



#Surrogates variables
out_list<-c("sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking","metabolic_syndrome", "alcohol_consumption",
            "age", "age","age","BMI", "ln_hscrp",
            "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb")

bin_names<-c("sex","diabetes","lipidmed", "blood_pressure_lowering_med", "current_smoking","metabolic_syndrome","alcohol_consumption",
             "high_age","middle_age","low_age","obesity", "high_hscrp",
             "high_triglycerides","high_ldl_chol","low_hdlchol","high_totchol","low_eGFR","low_wbc","low_hgb")
out_surro<-paste0("s_",out_list)
bin_surro<-paste0("s_",bin_names)

pheno_names<-c("sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking",
               "metabolic_syndrome", "alcohol_consumption",
               "age","BMI", "ln_hscrp","waist_circumference","weight","height",
               "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb")

#MetaboAge PARAMETERS
#load("data/PARAM__2018-06-18_02-16-17.457.RData")
#Surrogates PARAMETERS
#PARAM_surrogates<-readRDS("data/PARAM_surrogates_2021_03_26.RData")

#Set colors for the arrays of images
n <- 21
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
c21=colors[c(1:22)[-4]]
names(c21)<-c("mortScore", "MetaboAge", bin_surro)
###############################
## Mortality Score functions ##
###############################
## Defining a function to prepare metabolite data for Joris' mortality score:
prep.mort_data <- function(dat,featID=c("pufa_fa","gp","glc","s_hdl_l","xxl_vldl_l","alb","phe","acace","ile","vldl_d",
                                        "leu","val","his","lac"),quiet=FALSE){
  if(!quiet){
    cat("|| Preparing data ... \n")
  }
  ## 1. Check for zeroes:
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

## Defining a function to compute Joris' mortality score on metabolite data:
comp.mort_score <- function(dat,betas=mort_betas,quiet=FALSE){
  ## 1. Prepare data:
  if(!quiet){
    cat("=== Computing mortality score === \n")
  }
  prepped_dat <- prep.mort_data(dat,featID=betas$Abbreviation,quiet=quiet)
  ## 2. Compute:
  if(!quiet){
    cat("| Computing score .. \n")
  }
  mortScore <- as.data.frame(as.matrix(prepped_dat[,betas$Abbreviation]) %*% betas$Beta_value)
  colnames(mortScore)<-"mortScore"
  rownames(mortScore)<-rownames(dat)
  # if("LLnr" %in% colnames(dat)){
  #   names(mortScore) <- dat$LLnr
  # }
  if(!quiet){
    cat("Done!\n\n")
  }
  return(mortScore)
}

#########################
## MetaboAge functions ##
#########################
QCprep<-function(mat,PARAM,quiet=FALSE,Nmax_zero=1,Nmax_miss=1){
  ## mat is the input NH matrix; it may contain a mixture of flags and metabolites in the columns.
  ## PARAM is a list holding the parameters of the pipeline
  # 0. Start:
  if(!quiet){
    cat(report.dim(mat,header="Start"))
  }
  # 1. Subset required metabolites:
  mat <- subset.metabolites.overlap(mat,metabos=PARAM$MET,quiet=quiet)
  # 2. Subset samples on missingness:
  mat <- subset.samples.miss(mat,Nmax=Nmax_miss,quiet=quiet)
  # 3. Subset samples on zeros:
  mat <- subset.samples.zero(mat,Nmax=Nmax_zero,quiet=quiet)
  # 4. Subset samples on SD:
  mat <- subset.samples.sd(as.matrix(mat),MEAN=PARAM$logMEAN,SD=PARAM$logSD,quiet=quiet)
  # 5. Perform scaling:
  mat <- apply.scale(mat,MEAN=PARAM$MEAN,SD=PARAM$SD)
  # if(!quiet){
  #   #cat("| Performing scaling ... ")
  #   #cat(" DONE!\n")
  # }
  # 6. Perform imputation:
  mat <- impute.miss(mat)
  # if(!quiet){
  #   #cat("| Imputation ... ")
  #   #cat(" DONE!\n")
  # }
  return(mat)
}

apply.fit<-function(mat,FIT){
  # Resort:
  BETA <- FIT[colnames(mat)]
  INTC <- FIT[1]
  # Predict age:
  AGE <- data.frame(MetaboAge=as.vector(mat %*% BETA) + as.vector(INTC),stringsAsFactors=FALSE)
  rownames(AGE)<-rownames(mat)
  return(AGE)
}

subset.metabolites.overlap<-function(x,metabos,quiet=FALSE){
  x <- x[,intersect(colnames(x),metabos),drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header="Selecting metabolites"))
  }
  return(invisible(x))
}

subset.samples.miss<-function(x,Nmax=1,quiet=FALSE){
  MISS <- colSums(is.na(t(x)))
  x <- x[which(MISS<=Nmax),,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on missing values [Nmax>=",Nmax,"]")))
  }
  return(invisible(x))
}

subset.samples.zero<-function(x,Nmax=1,quiet=FALSE){
  ZERO <- colSums(t(x==0),na.rm=TRUE)
  x <- x[which(ZERO<=Nmax),,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on zero values [Nmax>=",Nmax,"]")))
  }
  return(invisible(x))
}

subset.samples.sd<-function(x,MEAN,SD,quiet=FALSE){
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

impute.miss<-function(x){
  ## This is an boiler-plate solution :)
  x[which(is.na(x))] <- 0
  return(x)
}


## Defining a function to apply the scaling parameters:
apply.scale <- function(dat,MEAN,SD,quiet=FALSE){
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

######################
## Surrogate models ##
######################
## Defining a function to calculate BMI LDL_chol and eGFR
BMI_LDL_eGFR<-function(phenotypes, metabo_measures){
  #cat("\n\n| Calculating BMI, LDL-cholesterol and eGFR using CKD-EPI\n\n")
  #BMIcalculation
  if(!c("BMI") %in% colnames(phenotypes)){
    if(all(c("weight","height") %in% colnames(phenotypes))){
      phenotypes<-cbind(phenotypes, data.frame(BMI=phenotypes$weight/(phenotypes$height/100)^2))
    }
  }
  
  
  #ldl-chol calculation
  if(!c("BMI") %in% colnames(phenotypes)){
    if(all(c("hdlchol","triglycerides", "totchol") %in% colnames(phenotypes))){
    phenotypes<-cbind(phenotypes, data.frame(
      ldl_chol= (phenotypes$totchol-(phenotypes$hdlchol + (0.45 * phenotypes$triglycerides)))))
    phenotypes[which(phenotypes$triglycerides>=4.52),"ldl_chol"]<-NA
    }
  }
  
  #eGFR calculation
  if(!c("BMI") %in% colnames(phenotypes)){
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

#Function to binarize all the phenotypes
binarize_all_pheno<-function(data){
  available_pheno<-intersect(colnames(data),c("sex","diabetes","lipidmed","current_smoking", "blood_pressure_lowering_med", "alcohol_consumption", "metabolic_syndrome"))
  binarized_phenotypes<-cbind(data[,available_pheno],
                              matrix(data=NA,nrow=dim(data)[1],ncol=12))
  colnames(binarized_phenotypes)<-c(available_pheno,
                                    "high_age", "middle_age","low_age","high_hscrp", "high_triglycerides","low_hdlchol","high_ldl_chol","high_totchol",
                                    "low_eGFR","obesity", "low_wbc","low_hgb")
  #hscrp
  if(c("hscrp") %in% colnames(data)){
    binarized_phenotypes[which(data$hscrp<2),"high_hscrp"]<-0
    binarized_phenotypes[which(data$hscrp>=2),"high_hscrp"]<-1
    binarized_phenotypes$high_hscrp<-as.factor(binarized_phenotypes$high_hscrp)
  }else if(c("ln_hscrp") %in% colnames(data)){
    hscrp<-exp(data[,"ln_hscrp"])
    binarized_phenotypes[which(hscrp<2),"high_hscrp"]<-0
    binarized_phenotypes[which(hscrp>=2),"high_hscrp"]<-1
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
  #lipidmed
  if(length(levels(binarized_phenotypes$lipidmed))>1){
    binarized_phenotypes[which(binarized_phenotypes$lipidmed==2),"lipidmed"]<-NaN
    binarized_phenotypes$lipidmed<-factor(binarized_phenotypes$lipidmed)
  }
  
  rownames(binarized_phenotypes)<-rownames(data)
  return(binarized_phenotypes)
}

QCprep_metabotypes<-function(mat,PARAM_surrogates,quiet=FALSE, MET=MET56, 
                             Nmax_miss=1,Nmax_zero=1){
  ## mat is the input NH matrix; it may contain a mixture of flags and metabolites in the columns.
  ## PARAM is a list holding the parameters of the pipeline
  # 0. Start:
  if(!quiet){
    cat(report.dim(mat,header="Start"))
  }
  # 1. Subset required metabolites:
  mat <- subset.metabolites.overlap(mat,metabos=MET56,quiet=quiet)
  # 2. Subset samples on missingness:
  mat <- subset.samples.miss(mat,Nmax=Nmax_miss,quiet=quiet)
  # 3. Subset samples on zeros:
  mat <- subset.samples.zero(mat,Nmax=Nmax_zero,quiet=quiet)
  # 4. Subset samples on SD:
  mat <- subset.samples.sd_no_log(as.matrix(mat),MEAN=as.matrix(PARAM_surrogates$mean),
                                  SD=as.matrix(PARAM_surrogates$sd),quiet=quiet)
  sample_names<-rownames(mat)
  # 5. Perform scaling:
  #if(!quiet){
    #cat("| Performing scaling ... ")
    mat <- sapply(PARAM_surrogates$MET,function(x) (mat[,x]-as.numeric(PARAM_surrogates$mean[x]))/as.numeric(PARAM_surrogates$sd[x]))
    rownames(mat)<-sample_names
    #cat(" DONE!\n")
  #}
  # 6. Perform imputation:
  #if(!quiet){
    #cat("| Imputation ... ")
    mat <- impute.miss(mat)
    #cat(" DONE!\n")
  #}
  return(mat)
}

calculate_surrogate_scores <- function(met, pheno, PARAM_surrogates, bin_names, bin_surro=paste0("s_",bin_names),Nmax_miss=1,Nmax_zero=1, sex_sel="all",log_file=NULL,roc=FALSE,quiet=FALSE){
  #QC of the metabolites
  metabo_measures<-QCprep_metabotypes(as.matrix(met[,MET63]), PARAM_surrogates, quiet=quiet,Nmax_miss=Nmax_miss,Nmax_zero=Nmax_zero)
  
  if(roc){
    sex_sel<-"all"
    
    phenotypes<-pheno[which(rownames(pheno)%in%rownames(metabo_measures)),]
    bin_pheno<-binarize_all_pheno(phenotypes)
    
    #BBMRI surrogates Calculation
    all(rownames(bin_pheno)==rownames(metabo_measures))
    
    surrogates<-foreach(i=bin_names, .combine="cbind") %do% {
      pred<-predict(PARAM_surrogates$models_betas[[i]]$finalModel, newx=as.matrix(metabo_measures[,rownames(PARAM_surrogates$models_betas[[i]]$finalModel$beta)]), 
                    s=PARAM_surrogates$models_betas[[i]]$best_tune$lambda, type= "response")
    }
    ind<-sapply(bin_names, function(x) !all(is.na(bin_phenotypes[x])))
    par(mfrow=c(5,4))
    surro_images<-foreach(i=bin_names[ind], .combine="cbind") %do% {
      a<-data.frame(out=factor(bin_pheno[,i]), metabo_measures)
      colnames(a)[1]<-"out"
      pred<-predictions_glmnet(PARAM_surrogates$models_betas[[i]], data=a,title_img=i)
      return(pred$predictions_eval)
    }
  }else{
    surrogates<-foreach(i=bin_names, .combine="cbind") %do% {
      pred<-predict(PARAM_surrogates$models_betas[[i]]$finalModel, newx=as.matrix(metabo_measures[,rownames(PARAM_surrogates$models_betas[[i]]$finalModel$beta)]), 
                    s=PARAM_surrogates$models_betas[[i]]$best_tune$lambda, type= "response")
    }
  }
  colnames(surrogates)<-bin_surro
  rownames(surrogates)<-rownames(metabo_measures)
  if(roc){
    return(list(dat=list(phenotypes=phenotypes,metabo_measures=metabo_measures,bin_phenotypes=bin_phenotypes),surrogates=surrogates, roc_curves=surro_images))
  }else{
    return(list(metabo_measures=met,surrogates=surrogates))
  }
}

subset.samples.sd_no_log<-function(x,MEAN,SD,quiet=FALSE){
  MEAN <- MEAN[colnames(x)]
  SD <- SD[colnames(x)]
  Dummi <- x
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

report.dim<-function(x,header,trailing="0"){
  return(paste0(sprintf(paste0("%-",trailing,"s"),paste0("| ",header,":\n")),sprintf("%4s",ncol(x))," metabolites x ",sprintf("%4s",nrow(x))," samples \n"))
}


predictions_glmnet<- function(en, data, folderIMG, log_file=NULL, img_name, title_img=FALSE, save=FALSE, peval=TRUE){
  if (is.factor(data[,"out"])){
    pred<-as.numeric(predict(en$finalModel, as.matrix(data[,rownames(en$finalModel$beta)]), s=en$best_tune$lambda, type="response"))
  }else{
    pred<-as.numeric(predict(en$finalModel, as.matrix(data[,rownames(en$finalModel$beta)]), s=en$best_tune$lambda))
  }
  all(rownames(pred) == rownames(data))
  if(peval){
    #Evaluation of the predictions
    if (is.factor(data[,"out"])){
      pred_eval<-eval_prediction_binary(data[,"out"], pred, log_file= log_file,
                                        title=title_img,
                                        img_name=paste0(img_name), save=save)
      pROC::auc(pred_eval)
    }else{
      pred_eval<-eval_prediction_continuous(data[,"out"], pred, 
                                            log_file= log_file, save=save)
      pred_eval
    }
    res<-list(predictions=pred, predictions_eval=pred_eval)
  }else{
    res<-list(predictions=pred)
  }
  return(res)
}
####################
## Plot functions ##
####################
### Defining a function for plotting a heatmap indcating missing & zero values:
plot.na.heatmap  <- function(dat){
  Dummi <- matrix(NA,ncol=ncol(dat),nrow=nrow(dat))
  Dummi[which(!is.na(dat==0))] <- 1
  #Dummi[which(dat!=0)] <- 1
  graphics::layout(mat=matrix(c(1,2,3,4),ncol=2,byrow=TRUE),widths=c(4,1),heights=c(1,4))
  par(xaxs = "i")  # THIS IS EVIL MAGIC!
  par(yaxs = "i")
  par(mar=c(0,3,1,0))
  XCOUNT <- nrow(Dummi)-colSums(Dummi,na.rm=TRUE)
  XPERC  <- XCOUNT/nrow(Dummi)*100
  YCOUNT <- ncol(Dummi)-rowSums(Dummi,na.rm=TRUE)
  YPERC  <- YCOUNT/ncol(Dummi)*100
  if(max(XPERC)<10){
    ylim <- c(0,10)
  } else {
    ylim <- c(0,max(pretty(XPERC)))
  }
  par(mar=c(0.5,13,2.5,0))
  barplot(XPERC,axes=TRUE,main="Missingness Plot",col="lightblue",border="lightblue",
          las=2,cex.main=3,ylim=ylim,cex.axis=1.2,font.axis=1.5)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1))
  legend("topright",fill=c("white","grey30"),legend=c("missing","value"))
  par(mar=c(2,13,0,0))
  image(t(Dummi[nrow(Dummi):1,]),axes=F, col=c("grey30"),xlim=c(0,1),ylim=c(0,1))
  mtext(text=rownames(dat), side=2, line=0.3, at=c(0.985,seq(0.945,0.055, length=dim(dat)[1]-2),0.03),
        las=1, cex=0.8)
  mtext(text="samples",side=1,font=1.5,cex=1.5, line = 0.6)
  par(mar=c(2,0,0,0))
  if(max(YPERC)<10){
    xlim <- c(0,10)
  } else {
    xlim <- c(0,max(pretty(YPERC)))
  }
  barplot(rev(YPERC),horiz=TRUE,axes=TRUE,col="lightblue",border="lightblue",cex.axis=1,font.axis=2)
}


#Evaluating predictions
eval_prediction_binary<-function(x,pred, log_file=FALSE, title, img_name=FALSE, save, color="#377eb8"){
  xx<-roc(x, as.numeric(pred), plot=TRUE, col=color, quiet = TRUE,
          lwd=4, print.auc=TRUE, main = title, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
  # if(save){
  #   dev.copy(png,img_name)
  #   dev.off()
  #   cat(paste("\nThe AUC of the predictions is:", pROC::auc(xx)), 
  #       file=log_file, append=TRUE)
  # }
  return(xx)
}


roc_surro_subplots<-function(surrogates, bin_phenotypes, log_file=FALSE, img_name=FALSE, save=F){
  
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
      ROC_curve<-roc(bin_phenotypes[,i], as.numeric(surrogates[,surro]), plot=F, col=c21[surro], quiet = TRUE,
              lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
      roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
    
      pl <- plot_ly(roc_df, x= ~Specificity, y= ~Sensitivity, 
                    type = 'scatter', mode = 'lines', line = list(color= c21[surro],width = 4)) %>%
        add_segments(x = 1, xend =0, y = 0, yend = 1, line = list(dash = "dash", width= 1, color="black")) %>% 
        layout(
          #title = list(text=paste0('<b>',i,'<b> AUC=<b>',round(as.numeric(ROC_curve$auc),digits = 3)),font=title_font, y = 0.98),
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
    #subplot(surro_images, nrows = 5, shareX = F,shareY = F, titleX = F, titleY = F, which_layout = 1)
    return(subplot(surro_images, nrows = 5, shareX = F,shareY = F, titleX = F, titleY = F, which_layout = 1))
}

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
  # ROC_curve<-eval_prediction_binary(bin_phenotypes[,selected[1]], surrogates[,paste0("s_",selected[1])],
  #                                   log_file=FALSE, color = c21[surro],
  #                                   title=selected[1], save=F)
  ROC_curve<-roc(bin_phenotypes[,selected[1]], as.numeric(surrogates[,paste0("s_",selected[1])]), plot=F, col=c21[surro], quiet = TRUE,
                 lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
  
  roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
  pl <- plot_ly(roc_df, x= ~Specificity, y= ~Sensitivity, 
                type = 'scatter', mode = 'lines', line = list(color= c21[surro],width = 4),
                name=paste0(surro, ":\nAUC=",round(as.numeric(ROC_curve$auc),digits = 3))) %>%
    add_segments(x = 1, xend =0, y = 0, yend = 1, line = list(dash = "dash", width= 1, color="black"),
                 name="random:\nAUC=0.5")
  
  if(length(selected)>1){
    surro_images<-foreach(i=2:length(selected), .combine="cbind") %do% {
      surro<-paste0("s_",selected[i])
  
      # ROC_curve<-eval_prediction_binary(bin_phenotypes[,selected[i]], surrogates[,paste0("s_",selected[i])],
      #                                   log_file=FALSE, color = c21[surro],
      #                                   title=selected[i], save=F)
      ROC_curve<-roc(bin_phenotypes[,selected[i]], as.numeric(surrogates[,paste0("s_",selected[I])]), plot=F, col=c21[surro], quiet = TRUE,
                     lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
      roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
      
      pl<-pl %>% add_trace(x= roc_df$Specificity, y= roc_df$Sensitivity, 
                           type = 'scatter', mode = 'lines', line = list(color= c21[surro],width = 4),
                           name=paste0(selected[i], ":\nAUC=",round(as.numeric(ROC_curve$auc),digits = 3)))
    }
  }
  
    if(length(sel_not_available)>0){
      if(length(sel_not_available)>1){
        ann<-paste("Not available:\n",noquote(paste(sel_not_available, collapse = ',\n')))
      }else{ann<-paste("Not available:\n",sel_not_available)}
      
      pl<-pl %>% layout(
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
      pl <- pl %>% layout(
        title = list(text=paste0('<b> ROC curves <b>'),font=title_font, y = 0.98),
        xaxis = list(autorange = "reversed", 
                     title = paste("<b> Specificity <b>")),
        yaxis = list (title = paste("<b> Sensitivity <b>"))
      )
    }
    return(pl)
}


ttest_surrogates<-function(surrogates,bin_phenotypes){
  available<-colnames(bin_phenotypes)[!colSums(is.na(bin_phenotypes))==nrow(bin_phenotypes)]
  
  Surrogates<-foreach(i=available, .combine="rbind") %do% {
    surro<-paste0("s_",i)
    comp<-data.frame(ID=rownames(surrogates), 
                     value=bin_phenotypes[rownames(surrogates),i],
                     surrogate=surrogates[,surro], variable=i, 
                     stringsAsFactors = F)
  }
  
  tt_values<-foreach(i=available, .combine='rbind') %do%{
    comp<-Surrogates[which(Surrogates$variable==i),]
    if(length(which(is.na(comp[,"value"])))>0){
      comp<-comp[-which(is.na(comp[,"value"])),]
    }
    comp$value<-factor(comp$value)
    tt<-t.test(comp[which(comp$value==1),"surrogate"],comp[which(comp$value==0),"surrogate"], var.equal = TRUE )
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

  pl_surro<-plot_ly(Surrogates,
              x = ~variable,
              y = ~surrogate,
              color = ~value,
              type = "box",
              colors = c("#377EB8","#E41A1C")) %>% 
    layout(boxmode = "group",
           title = list(text="<b>Surrogates' distributions split for their original values<b>",y = 0.98),
           xaxis = list(title = "<b>Clinical Variables<b>",
                        zeroline = FALSE),
           yaxis = list(title = "<b>Surrogate values<b>",
                        zeroline = FALSE)
           )
  
  return(pl_surro)
}



scatterplot_predictions <- function(x, p, title, xname, yname) {
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
    fit <- lm(predicted ~ original, data = df)
    cordf<-cor(df$original, df$predicted, method = c("pearson"))
    
    cPlot <- plot_ly(df, x= ~original, y=~predicted, mode="markers", alpha=0.7) %>% 
      add_markers(y= ~predicted) %>% 
      add_lines(x=~original, y=predict(fit), name="Regression line", 
                line = list( width= 3, color="#E41A1C")) %>% 
    add_segments(x = min(df$original), xend = max(df$original), 
                 y = min(df$predicted), yend = max(df$predicted), 
                 name="x=y", line = list(dash = "dash", width= 1, color="black"))%>%
      layout(title = list(text=paste('<b>',title,'<b>'),font=title_font, y = 0.98),
             xaxis = list(title = paste('<b>',xname,'<b>'),titlefont = axis_font),
             yaxis = list (title = paste('<b>',yname,'<b>'),titlefont = axis_font),
             # #margin = list(l = 50, r = 140, b = 30, t = 30, pad = 4),
             annotations = list(text = paste0("R= ", as.character(round(cordf, digits= 3)),", R^2= ", as.character(round((cordf)^2, digits= 3)),
                                  "\n med. error=", as.character(round(median(abs(df$original - df$predicted)), digits=3))),
                                xref='paper',
                                yref='paper',
                                x = 1.24, y = 0.25,
                                showarrow=FALSE,
                                bordercolor=c("black"),
                                borderwidth=2
                                )
             )
    return(cPlot)
}



#Function to plot the barplots with density
metabo_barplots<-function(metabo_measures,color){
  par(mfrow = c(3, 4))
  for (i in 1:56) { # Loop over loop.vector
    hist(na.omit(metabo_measures[,i]), # histogram
         border="black",
         prob = TRUE, # show densities instead of frequencies
         breaks = 1000,
         xlab="",
         main = colnames(metabo_measures)[i])
    lines(density(na.omit(metabo_measures[,i])), # density plot
          lwd = 3, # thickness of line
          col = color)
  }
  
}

one_hist<-function(dat, color=c21, scaled=FALSE){
  
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
    plot<-plot_ly(x =~dat[,i],opacity = 1, type = "histogram", marker = list(color = color[i]), name = i)%>% 
      layout(
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
  return(subplot(hist_images, nrows = 5, shareX = F,shareY = F, titleX = F, titleY = F, which_layout = 1))
}

hist_plots<-function(dat,x_name, color=c21, scaled=FALSE, main="Predictors Distributions"){
  # bin<-max(dat[,"x"],na.rm = T)/100
  # plot<-ggplot(data=dat, aes(x=x)) +
  #   geom_histogram(binwidth=bin,position="identity", fill=color[x_name],alpha=1)+
  #   theme(plot.title = element_text(size=16, face= "bold"),
  #         axis.title.x = element_text(size=14, face= "bold"),
  #         axis.title.y = element_text(size=14, face= "bold"),
  #         axis.text.x = element_text(hjust = 1, size = 10),
  #         axis.text.y = element_text(hjust = 1, size = 10),
  #         legend.position = "none")+
  #   labs(title=paste("Distribution of",x_name),
  #        x=x_name)
  dat<-as.matrix(dat[,x_name])
  colnames(dat)<-x_name
  if(scaled){
    if(length(x_name)>1){
      sd<-colSds(dat,na.rm = T)
      names(sd)<-colnames(dat)
      dat<-apply.scale(dat,MEAN = colMeans(dat,na.rm = T), SD=sd)
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
  plot<-plot_ly(x = ~dat[,x_name[1]],opacity = 0.45, type = "histogram", marker = list(color = color[x_name[1]]), name = x_name[1])
  x_title<-"score"
   if(length(x_name)>1){
    for(i in c(2:length(x_name))){
      plot <- plot %>% add_histogram(x = dat[,x_name[i]], marker = list(color = color[x_name[i]]), name = x_name[i])
    }
     x_title<-"scores"
  }
  
  plot<-plot %>% layout( title =list(text=paste("<b>",main,"<b>"),font=title_font,y = 0.98),
                         xaxis = list(
                           title = paste("<b>",x_title,"<b>"),
                           titlefont = axis_font), 
                         yaxis = list(
                           title = "<b>Counts<b>",
                           titlefont = axis_font),
                         barmode = "overlay")
  print(plot)
}
  


hist_plots_mortality<-function(mort_score,phenotypes){
  age<-data.frame(age=phenotypes[rownames(mort_score),"age"])
  age[which(age$age<45),"age_ranges"]<-"age<45"
  age[intersect(which(age$age>=45),which(age$age<65)),"age_ranges"]<-"45>=age<65"
  age[which(age$age>=65),"age_ranges"]<-"age>=65"
  
  mort_score$age_ranges<-age$age_ranges
  colnames(mort_score)[2]<-"age_ranges"
  mort_score$age_ranges<-factor(mort_score$age_ranges, levels=c("age<45","45>=age<65","age>=65"))
  
  mu<-ddply(mort_score,~age_ranges,summarise,mean=mean(mortScore,na.rm=T))
  
  axis_font <- list(
    family = "Arial",
    size = 14
  )
  title_font <- list(
    family = "Arial",
    size = 18,
    margin=10
  )
  fig <- plot_ly(alpha = 0.4) 
    if(length(which(mort_score$age_ranges=="age<45"))>0){
      fig<-fig %>% 
        add_histogram(x = ~mort_score$mortScore[which(mort_score$age_ranges=="age<45")], name="age<45", marker = list(color = "#B3DE69"), opacity=0.7) %>% 
        add_segments(x = mu[1,2], xend = mu[1,2], y = 0, yend = 100, name="mean(age<45)")
    }
  if(length(which(mort_score$age_ranges=="45>=age<65"))>0){
    fig<-fig %>% 
      add_histogram(x = ~mort_score$mortScore[which(mort_score$age_ranges=="45>=age<65")], name="45>=age<65", marker = list(color = "#80B1D3"), opacity=0.7) %>% 
      add_segments(x = mu[2,2], xend = mu[2,2], y = 0, yend = 100, name="mean(45>=age<65)")
  }
  if(length(which(mort_score$age_ranges=="age>=65"))>0){
    fig<-fig %>% 
      add_histogram(x = ~mort_score$mortScore[which(mort_score$age_ranges=="age>=65")], name="age>=65", marker = list(color = "#FB8072"), opacity=0.7) %>% 
      add_segments(x = mu[3,2], xend = mu[3,2], y = 0, yend = 100, name="mean(age>=65)")
  }  
  
    fig<-fig %>% layout(title =list(text=paste("<b> Mortality divided in age ranges <b>"),font=title_font,y = 0.98),
           xaxis = list(
             title = paste("<b> Counts <b>"),
             titlefont = axis_font), 
           yaxis = list(
             title = "<b> Mortality score <b>",
             titlefont = axis_font),barmode = "overlay")
  fig
}


cor.assoc <- function(dat,dat2,featID,phenID,covID=NULL,method="pearson",quiet=FALSE){
  res <- lapply(featID, function(FEAT){
    t(sapply(phenID,function(PHEN){
      unlist(cor.test(x=dat[,FEAT],y=dat2[,PHEN],method=method)[c("estimate","p.value","statistic")])
    }))
  })
  names(res) <- featID
  return(res)
}


## Defining a SubFunction for plotting feature x phenotype associations on basis of (partial) correlations:
plot.cor.assoc <- function(res,main=NULL,zlim=NULL,reorder.x=FALSE,reorder.y=reorder.x,resort_on_p=FALSE,
                           sep.x=NULL,abs=FALSE,cor.abs=FALSE,text=TRUE,value=FALSE,cex.lab=1, cex.num=8, 
                           myColor=NULL,myBreaks=NULL, annotations=NULL, gaps=NULL,  mycolors=NULL, fontsize=8){
  #suppressPackageStartupMessages(require(WGCNA))
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
  if(is.null(zlim)){
    MAX_s <- max(abs(c(min(s),max(s))))
    zlim <- c(min(pretty(c(-MAX_s,MAX_s))),max(pretty(c(-MAX_s,MAX_s))))
  }
  if(cor.abs){s<-abs(s)}
  if(text){
    if(value){
      if(dim(s)[1]>35 || dim(s)[2]>35){
        textMatrix <- paste(signif(s, 2), "(",signif(p, 2), ")", sep = "")
      }else{
        textMatrix <- paste(signif(s, 2), "\n(",signif(p, 2), ")", sep = "")
      }
    } else {
      textMatrix <- rep("",(nrow(p)*ncol(p))) #it repeats the "" for all the possible values
      IND1 <- which(p<=(0.05/(nrow(p)*ncol(p)))) #it compares the p-values with the Bonferroni corrected threshold
      IND2 <- which(p.adjust(p,method="fdr")<=0.05) #it compares the p-values adjusted with the fdr with thresh=0.05
      #compose the textMatrix with * if they are selected with the bonferroni thresh, 
      #. if they are selected using the fdr but not with bonferroni
      #X if they were not selected on the previous steps
      textMatrix[IND1] <- "**"
      textMatrix[setdiff(IND2,IND1)] <- ".."
      textMatrix[which(is.na(p))] <- "X"
    }
    dim(textMatrix) <- dim(p)
    if(is.sym(res)){
      textMatrix[cbind(1:nrow(p),1:ncol(p))] <- "X"
    }
  } else {
    textMatrix <- NULL
  }
  if(text){
    pheatmap(s, scale = "none", main=main, display_numbers = textMatrix,
                 cluster_rows = F, cluster_cols = F, margin = c(5,5), fontsize_number = cex.num, cex=1.2,
                 show_rownames=1, show_colnames=1, colors = colorRampPalette(c("blue","white", "red"))(20))
  }else{
    pheatmap(s, scale = "none", main=main,
             cluster_rows = F, cluster_cols = F, margin = c(5,5), fontsize_number = cex.num, cex=1.2,
             fontsize=fontsize,
             show_rownames=1, show_colnames=1, colors = colorRampPalette(c("blue","white", "red"))(20))
  }
}



plot.corply <- function(res,main=NULL,zlim=NULL,reorder.x=FALSE,reorder.y=reorder.x,resort_on_p=FALSE,
                           sep.x=NULL,abs=FALSE,cor.abs=FALSE,reorder_dend=F){
  #suppressPackageStartupMessages(require(WGCNA))
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
    heatmaply_cor(s,main=main, margins=c(0,0,40,0))
  }else{
    heatmaply_cor(s,Rowv=F, Colv=F, main=main,margins=c(0,0,40,0))
  }
  
}

plot.vis.mis <- function(dat){
  vis_miss(dat)
}

## Defining a helper function for clustering results:
resort.on.s <- function(res,abs=FALSE){
  s <- get.s(res)
  if(abs){
    s <- abs(s)
  }
  order.x <- hclust(as.dist(1-(cor(s))))$order
  if(is.sym(res)){
    order.y <- order.x
  } else {
    order.y <- hclust(as.dist(1-(cor(t(s)))))$order
  }
  return(list(x=order.x,y=order.y))
}

resort.on.p <- function(res,abs=FALSE){
  p <- get.p(res)
  if(abs){
    p <- abs(p)
  }
  order.x <- hclust(as.dist(1-(cor(log(p+1e-323)))))$order
  if(is.sym(res)){
    order.y <- order.x
  } else {
    order.y <- hclust(as.dist(1-(cor(t(log(p+1e-323))))))$order
  }
  return(list(x=order.x,y=order.y))
}

## Defining a helper function for extracting statistics from results:
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

## Defining a helper function for extracting pvalues from results:
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


## Defining a helper function to check whether a supplied matrix is symmetric:
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

# a<-CalibratR::calibrate(as.numeric(LLS_pheno$sex), as.numeric(surro[,"s_sex"]), model_idx = c(1, 2, 3, 4, 5),
#           evaluate_no_CV_error = TRUE, evaluate_CV_error = TRUE, folds = 5)
# library(eRic)
# platt<-prCalibrate(r.calib=as.numeric(LLS_pheno$sex), p.calib= as.numeric(surro[,"s_sex"]), nbins = 50)
# 
# r.calib<-as.numeric(LLS_pheno$sex)
# p.calib<-as.numeric(surro[,"s_sex"])
# cmodel <- glm(y ~ x, data.frame(y = r.calib, x = p.calib), 
#               family = "binomial")
# calibrated <- predict(cmodel, data.frame(y = r.calib, x = p.calib), 
#                       type = "response")
# all(platt$cal.probs == calibrated)
# 
# ece1 <- getECE(r.calib, p.calib, 50)
# ece2 <- getECE(r.calib, calibrated, 50)

getECE<-function (actual, predicted, n_bins = 10) 
{
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
    expected <- mean(group_pred)
    observed <- mean(group_actual)
    S[i] <- abs(observed - expected)
    W[i] <- n_/N
    groups[[i]] <- group_pred
  }
  mean_prediction <- lapply(groups, mean)
  min_group <- lapply(groups, min)
  max_group <- lapply(groups, max)
  res <- t(S) %*% W
  return(as.numeric(res))
}

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
    expected <- mean(group_pred)
    observed <- mean(group_actual)
    S[i] <- abs(observed - expected)
    W[i] <- n/N
  }
  res <- max(S * W)
  return(res)
}

plattCalibration<-function (r.calib, p.calib, r.valid = NULL, p.valid = NULL, 
          nbins = 10, plot=FALSE) {
  
  #cat("\n\n| Performing calibration!!\n\n")
  if (!require(ggplot2) || !require(gridExtra)) {
    stop("packages [ggplot2, gridExtra] are required for this function")
  }
  # if (!is.numeric(r.calib) || (length(unique(r.calib)) != 2) ||
  #     (min(r.calib) != 0) || (max(r.calib) != 1)) {
  #   print("bubi")
  #   stop("[r.calib] must be a binary numeric vector with no missing values")
  # }
  # if (!is.null(r.valid) && (!is.numeric(r.valid) || (length(unique(r.valid)) != 
  #                                                    2) || (min(r.valid) != 0) || (max(r.valid) != 1))) {
  #   stop("[r.valid] must be a binary numeric vector with no missing values if included")
  # }
  # if (!is.numeric(p.calib) || length(p.calib[is.na(p.calib)]) > 
  #     0 || min(p.calib) < 0 || max(p.calib) > 1) {
  #   stop("[p.calib] must be a numeric vector with values between [0,1]")
  # }
  # if (!is.null(p.valid) && (!is.numeric(p.valid) || length(p.valid[is.na(p.valid)]) > 
  #                           0 || min(p.valid) < 0 || max(p.valid) > 1)) {
  #   stop("[p.valid] must be a numeric vector with values between [0,1] if included")
  # }
  # if (!is.numeric(nbins) || length(nbins) != 1 || nbins != 
  #     round(nbins) || nbins <= 0) {
  #   stop("[nbins] must be an integer value greater than zero")
  # }
  if (is.null(r.valid) || is.null(p.valid)) {
    pred <- p.calib
    resp <- r.calib
    sample <- "calibration"
  }else{
    pred <- p.valid
    resp <- r.valid
    sample <- "validation"
  }
  pred[pred == 0] <- 1e-08
  pred[pred == 1] <- 1 - 1e-08
  cmodel <- glm(y ~ x, data.frame(y = r.calib, x = p.calib), 
                family = "binomial")
  calibrated <- predict(cmodel, data.frame(y = resp, x = pred), 
                        type = "response")
  raw.bins <- cut(pred, nbins, include.lowest = TRUE)
  raw.xval <- tapply(pred, raw.bins, mean)
  raw.yval <- tapply(resp, raw.bins, mean)
  raw.cali <- data.frame(method = rep("Original", nbins), 
                         x = raw.xval, y = raw.yval)
  raw.logl <- (-1/length(resp)) * sum(resp * log(pred) + (1 - 
                                                            resp) * (log(1 - pred)), na.rm = TRUE)
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
  
  if(plot){
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
    
    cPlot <- plot_ly(cdata, x= ~x.original, y= ~y.original, name="Original", 
                     type = 'scatter', mode = 'lines+markers', line = list(width = 2), legendgroup = ~x.original) %>% 
      add_trace(x= ~x.calibrated, y= ~y.calibrated, name = 'Calibrated', 
                mode = 'lines+markers', color= I("red"), line = list(width = 2)) %>% 
      add_segments(x = 0, xend = 1, y = 0, yend = 1, name="Perfect", line = list(dash = "dash", width= 1, color="black")) %>%
      layout(title = list(text='<b>Reliability Diagram<b>',font=title_font, y = 0.98),
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
    hPlot<-plot_ly(hdata, x = ~pred, opacity = 0.45, type = "histogram", 
                   marker = list(color = "blue"), name = "Original",legendgroup = ~pred) %>% 
      add_histogram(x = calibrated, marker = list(color = "red"), 
                    name = "Calibrated") %>% 
      layout( title =list(text=paste("<b> Distributions <b>"),font=title_font,y = 0.98),
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

plattCalib_plot<-function(calibration,name, nbins = 10, annot_x=c(1.27,1.27),annot_y=c(0.6,0.4)){
  resp<-calibration$responses
  pred<-calibration$raw.probs
  calibrated<-calibration$cal.probs
  cmodel<-calibration$cal.model
  
  raw.bins <- cut(pred, nbins, include.lowest = TRUE)
  raw.xval <- tapply(pred, raw.bins, mean)
  raw.yval <- tapply(resp, raw.bins, mean)
  raw.cali <- data.frame(method = rep("Original", nbins), 
                         x = raw.xval, y = raw.yval)
  raw.logl <- (-1/length(resp)) * sum(resp * log(pred) + (1 - 
                                                            resp) * (log(1 - pred)), na.rm = TRUE)
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
  
  cPlot <- plot_ly(cdata, x= ~x.original, y= ~y.original, name="Original", 
                   type = 'scatter', mode = 'lines+markers', line = list(width = 2)) %>% 
    add_trace(x= ~x.calibrated, y= ~y.calibrated, name = 'Calibrated', 
              mode = 'lines+markers', color= I("red"), line = list(width = 2)) %>% 
    add_segments(x = 0, xend = 1, y = 0, yend = 1, name="Perfect", line = list(dash = "dash", width= 1, color="black")) %>%
    layout(title = list(text=paste('<b>Reliability Diagram<b>',name,'<b>'),font=title_font, y = 0.98),
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
             x =annot_x, y = annot_y,
             showarrow=FALSE,
             bordercolor=c('blue','red'),
             borderwidth=2))
  
  hdata<-data.frame(pred, calibrated)
  hPlot<-plot_ly(hdata, x = ~pred, opacity = 0.45, type = "histogram", nbinsx = nbins,
                 marker = list(color = "blue"), name = "Original") %>% 
    add_histogram(x = calibrated, marker = list(color = "red"), 
                  name = "Calibrated") %>% 
    layout( title =list(text=paste("<b> Distributions <b>",name,'<b>'),font=title_font,y = 0.98),
            xaxis = list(
              title = paste("<b> Predicted Probability <b>"),
              titlefont = axis_font), 
            yaxis = list(
              title = "<b> Counts <b>",
              titlefont = axis_font),barmode = "overlay")

  return(list(cal.Plot = cPlot, prob.hist = hPlot))
}

calibration_surro<-function(bin_phenotypes, surrogates, bin_names, bin_surro, bin_pheno_available, plot=F){
  calib<-lapply(1:length(bin_surro), function(i){
    orig<-as.numeric(bin_phenotypes[,bin_names[i]])-1
    names(orig)<-rownames(bin_phenotypes)
    pred<-as.numeric(surrogates[,bin_surro[i]])
    names(pred)<-rownames(surrogates)
    ind<-which(!is.na(orig))
    if(length(ind)!=0){
      calibration<-plattCalibration(orig[ind],pred[ind], 10, plot=plot)
      return(calibration)
    }
  })
  
  calib_df<-calib_data_frame(calib, bin_phenotypes, bin_pheno_available)

  calib<-append(calib,list(calib_df))
  names(calib)<-bin_names
  return(calib)
}

calib_data_frame<-function(calibrations, bin_phenotypes, bin_pheno_available){
    cal<-matrix(NA, nrow = dim(bin_phenotypes)[1], ncol = length(bin_pheno_available))
    rownames(cal)<-rownames(bin_phenotypes)
    colnames(cal)<-bin_surro[bin_pheno_available]
    for(i in 1:length(bin_pheno_available)){
      cal[names(calibrations[[bin_pheno_available[i]]]$cal.probs),i]<-calibrations[[bin_pheno_available[i]]]$cal.probs
    }
    colnames(cal)<-bin_pheno_available
  return(as.data.frame(cal))
}

pheno_NA<-function(){
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
  
   plot_ly(x=0,y=0, opacity=0, type = "histogram")%>% 
    layout(title = list(text='<b>Phenotype not available!<b>',font=title_font, y = 0.98),
                        xaxis = ax, yaxis = ax)
}

phenos_NA<-function(){
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
  
  plot_ly(x=0,y=0, opacity=0, type = "histogram")%>% 
    layout(title = list(text='<b>Phenotypes not available, please upload them to get this image!<b>',font=title_font, y = 0.98),
           xaxis = ax, yaxis = ax)
}

met_NA<-function(){
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
  
  plot_ly(x=0,y=0, opacity=0, type = "histogram")%>% 
    layout(title = list(text='<b>Metabolites are missing, please check your upload!<b>',font=title_font, y = 0.98),
           xaxis = ax, yaxis = ax)
}



met_NA_image<-function(){
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, "Metabolites are missing, please check your upload!",cex = 1.6, col = "black", font=2)
}
pheno_NA_image<-function(){
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, "Phenotypes not available, please upload them to get this image!",cex = 1.6, col = "black", font=2)
}

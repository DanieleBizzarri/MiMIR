###################################
### DEFINITIONS / CONSTANTS #######
###################################
# Metabolite markers names from Nightingale health
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


BBMRI_METABO <- c("serum_c", "X","remnant_c", "vldl_c", "X","ldl_c", "hdl_c", "hdl2_c", "hdl3_c", "serum_tg", "vldl_tg", "ldl_tg", "hdl_tg",
                  "X", "X", "X", "X", "estc", "X", "X","X", "freec", "X", "X","X", "X", "X", "X","X","X", "X", "X","X","vldl_d", "ldl_d", "hdl_d",
                  "totpg", "tg_pg", "totcho", "pc", "sm", "apob", "apoa1", "apob_apoa1", "totfa", "unsatdeg","faw3", "faw6", 
                  "pufa", "mufa", "sfa", "la", "dha", "faw3_fa", "faw6_fa",  "pufa_fa", "mufa_fa", "sfa_fa", "la_fa","dha_fa", 
                  "X", "X", "ala", "gln", "glycine", "his", "X",  "ile", "leu", "val", "phe", "tyr", "glc", "lac", "pyr",
                  "cit", "glycerol", "bohbut", "ace", "acace", "X", "crea", "alb", "gp",
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
                  "s_hdl_pl_percentage",    "s_hdl_c_percentage",    "s_hdl_ce_percentage",    "s_hdl_fc_percentage",    "s_hdl_tg_percentage")

# Metabolite markers often deemed to be more or less independent and not numerically derived from eachother:
MET63 <- tolower(c("Ala","Gln","His","Phe","Tyr","Ile","Leu","Val","Glc","Lac","Pyr","Cit","Ace",
                   "AcAce","bOHBut","Crea","Alb","Gp","XXL_VLDL_L","XL_VLDL_L","L_VLDL_L","M_VLDL_L",
                   "S_VLDL_L","XS_VLDL_L","IDL_L","L_LDL_L","M_LDL_L","S_LDL_L","XL_HDL_L","L_HDL_L",
                   "M_HDL_L","S_HDL_L","IDL_C","Serum_C","VLDL_C","LDL_C","HDL_C","HDL2_C","HDL3_C",
                   "VLDL_D","LDL_D","HDL_D","Serum_TG","TotPG","PC","SM","TotCho","ApoA1","ApoB",
                   "TotFA","DHA","LA","FAw3","FAw6","PUFA","MUFA","SFA","FAw3_FA","FAw6_FA","PUFA_FA",
                   "MUFA_FA","SFA_FA","UnsatDeg"))

# The final set of metabolites used by the surrogates
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

MET62 <- tolower(c("Ala","Gln","His","Phe","Tyr","Ile","Leu","Val","Glc","Lac","Cit","Ace",
                   "AcAce","Crea","Alb","Gp","M_VLDL_L",
                   "S_VLDL_L","XS_VLDL_L","IDL_L","L_LDL_L","M_LDL_L","S_LDL_L",
                   "M_HDL_L","S_HDL_L","IDL_C","Serum_C","VLDL_C","LDL_C","HDL_C","HDL2_C","HDL3_C",
                   "VLDL_D","LDL_D","HDL_D","Serum_TG","TotPG","PC","SM","TotCho","ApoA1","ApoB",
                   "TotFA","DHA","LA","FAw3","FAw6","PUFA","MUFA","SFA","FAw3_FA","FAw6_FA","PUFA_FA",
                   "MUFA_FA","SFA_FA","UnsatDeg","XXL_VLDL_L", "l_vldl_ce_percentage", "l_hdl_fc", 
                   "apob_apoa1", "faw6_faw3", "glycine"))

#The final set of metabolites used by the mortality score
MET14<-c("pufa_fa","gp","glc","s_hdl_l","xxl_vldl_l","alb","phe","acace","ile","vldl_d",
               "leu","val","his","lac")

MET_COVID<-c("gp","dha","crea","mufa", "apob_apoa1","tyr","ile","sfa_fa","glc","lac","faw6_faw3",
             "phe", "serum_c", "faw6_fa","ala","pufa","glycine","his","pufa_fa","val","leu",
             "alb","faw3","ldl_c","serum_tg")

MET_T2D<-c("phe", "l_vldl_ce_percentage", "l_hdl_fc")

MET_CVD<-c("phe","mufa_fa","faw6","dha")

# Mortality score betas
mort_betas <- data.frame(
  Abbreviation=c("pufa_fa","gp","glc","s_hdl_l","xxl_vldl_l",
                 "alb","phe","acace","ile","vldl_d",
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

# Ahola Olli T2D score betas
Ahola_Olli_betas<-data.frame(
  Abbreviation=c("phe", "l_vldl_ce_percentage", "l_hdl_fc", "sex", "age", "BMI", "glucose"),
  Metabolite=c("phenylanine","Cholesteryl esters to total lipids ratio in large VLDL", 
               "free cholesterol in large HDL", "sex", "age", "BMI", "glucose"),
  Beta_value=c(0.32,-0.474,-0.321, -0.364, 0.0286, 0.462, 0.640),
  stringsAsFactors = FALSE)

# CVD score betas
CVD_score_betas<-data.frame(
  Abbreviation=c("phe","mufa_fa","faw6","dha","sex","systolic_blood_pressure","current_smoking","diabetes",
                 "blood_pressure_lowering_med","lipidmed","totchol","hdlchol"),
  Metabolite=c("phenylanine","mufa ratio to total fatty acids", 
               "omega6 fatty acids"," dha", "sex","systolic_blood_pressure","current_smoking","diabetes",
               "blood_pressure_lowering_med","lipidmed","totchol","hdlchol"),
  Beta_value=c(0.225,-0.145,-0.137,0.1,-0.613,0.145,0.379,0.665,0.103,0.295,0.247,-0.315),
  stringsAsFactors = FALSE)

covid_betas <- data.frame(
  Abbreviation=c("gp","dha","crea","mufa", "apob_apoa1","tyr","ile","sfa_fa","glc","lac","faw6_faw3",
                 "phe", "serum_c", "faw6_fa","ala","pufa","glycine","his","pufa_fa","val","leu",
                 "alb","faw3","ldl_c","serum_tg"),
  Metabolite=c("Glycoprotein acetyls", "docosahexaenoic acid", "creatinine","Monounsaturated fatty acids",
               "Apolipoprotein B/ Apolipoprotein A1","Tyrosine","Isoleucine","Ratio of saturated fatty acids to total fatty acids",
               "Glucose","Lactate","Ratio of omega-6 fatty acids to omega 3 fatty acids","Phenylanine", "Total cholesterol",
               "Ratio of omega-6 fatty acids to total fatty acids","Alanine","Polyunsaturated fatty acids",
               "Glycine","Histidine","Ratio of polyunsaturated fatty acids to total fatty acids",
               "Valine","Leucine","Albumin","Omega-3 fatty acids","Total cholesterol in LDL","Serum total triglycerides"),
  Beta_value=c(0.3713,0.2533, 0.2170, 0.1693, 0.1388, 0.1375, 0.1091, 0.0965, 0.0928,
               0.0772, 0.0642, 0.0289, -0.0116, -0.0493, -0.0498, -0.0578, -0.0648,
               -0.0987, -0.1404, -0.1812, -0.1844, -0.1914, -0.2208, -0.2466, -0.2652),
  stringsAsFactors = FALSE)

#Phenotypic variables names
pheno_names<-c("sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking",
               "metabolic_syndrome", "alcohol_consumption",
               "age","BMI", "hscrp","waist_circumference","weight","height",
               "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb",
               "glucose", "systolic_blood_pressure", "Event", "EventAge")

#Variables names related to each surrogate
out_list<-c("sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking",
            "metabolic_syndrome", "alcohol_consumption",
            "age", "age","age","BMI", "hscrp",
            "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb")
out_surro<-paste0("s_",out_list)

#Binarized variables names related to each surrogate
bin_names<-c("sex","diabetes","lipidmed", "blood_pressure_lowering_med", "current_smoking",
             "metabolic_syndrome","alcohol_consumption",
             "high_age","middle_age","low_age","obesity", "high_hscrp",
             "high_triglycerides","high_ldl_chol","low_hdlchol","high_totchol","low_eGFR","low_wbc","low_hgb")
bin_surro<-paste0("s_",bin_names)

#Set colors for the arrays of images
c21<-c("#7FC97F", "#BEAED4", "#FDC086", "#386CB0", "#F0027F", "#BF5B17", 
       "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
       "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A",
       "#33A02C", "#FB9A99", "#E31A1C")
names(c21)<-c("mortScore", "MetaboAge", bin_surro)

#Global variables
utils::globalVariables(c("bin_phenotypes", "i", "y", "mortScore", 
                         "surro","AUC","outcome", "Uploaded", "biobank",
                         "BBMRI_hist_scaled","BBMRI_hist"))

#Import packages
#' @import foreach
#' @importFrom purrr map_chr
#' @importFrom purrr map_lgl
#' @import dplyr
#' @import shiny
#' @import caret
#' @import ggfortify
#' @import survival
#' @import limma
NULL

######################
## Upload functions ##
######################
#' Function to translate metabolites names alternatives to the BBMRI names
#'
#' @param names The NH-metabolomics names
#' @return The uploaded names and the BBMRI names of the metabolites available in BBMRI.
#' @export

find_BBMRI_names<-function(names){
  new_names <- names %>% purrr::map_chr(function(id) {
    # Look through the alternative_ids
    hits <-
      purrr::map_lgl(
        metabo_names_translator$alternative_names,
        ~ id %in% .
      )
    
    # If one unambiguous hit, return it.
    if (sum(hits) == 1L) {
      return(metabo_names_translator$BBMRI_names[hits])
      # If not found, give a warning and pass through the input.
    } else {
      warning("Biomarker not found: ", id, call. = FALSE)
      return(id)
    } 
  })
  n<-data.frame(uploaded=names,BBMRI_names=new_names)
  #n<-n[which(n$BBMRI_names %in% metabo_names_translator$BBMRI_names),]
  return(n)
}

###############################
## Mortality Score functions ##
###############################
#' Function to prepare the metabolites data for the mortality, T2D and CVD scores.
#'
#' @param dat The NH-metabolomics matrix
#' @param featID String vector with the names of metabolic features used
#' @param plusone TRUE or FALSe if 1.0 should be  added to all metabolites to calculate the score
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @return The NH-metabolomics matrix after checking for zeros and zscale the log metabolites concentrations
#' @export

prep_met_for_scores <- function(dat,featID,plusone=F, quiet=FALSE){
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

#' Function to compute Deelen et al. mortality score on metabolite data.
#'
#' @param dat The NH-metabolomics matrix 
#' @param betas The betas of the linear regression composing the mortality score
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @return The mortality score by Deelen et al.
#' @export
#' 
comp.mort_score <- function(dat,betas=mort_betas,quiet=FALSE){
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

###############################
## Other scores functions ##
###############################
#' Function to compute T2D score by Ahola Olli mortality score on metabolite data.
#'
#' @param dat The NH-metabolomics matrix 
#' @param betas The betas of the linear regression composing the mortality score
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @return The T2D by Ahola Olli
#' @export
#' 
comp.T2D_Ahola_Olli<- function(met,phen,betas=Ahola_Olli_betas,quiet=FALSE){
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

#' Function to compute CVD score by Peter Wurtz on metabolite data.
#'
#' @param dat The NH-metabolomics matrix 
#' @param betas The betas of the linear regression composing the score
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @return The CVD score by Peter Wurtz
#' @export
#' 
comp.CVD_score <- function(met,phen,betas=CVD_score_betas,quiet=FALSE){
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

#' Function to prepare the metabolites data for the COVID score
#'
#' @param dat The NH-metabolomics matrix
#' @param featID String vector with the names of metabolic features used
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @return The NH-metabolomics matrix after checking for zeros and zscale the log metabolites concentrations
#' @export
#' 
prep_data_COVID_score <- function(dat,featID=c("gp","dha","crea","mufa", "apob_apoa1","tyr","ile","sfa_fa","glc","lac","faw6_faw3",
                                               "phe", "serum_c", "faw6_fa","ala","pufa","glycine","his","pufa_fa","val","leu",
                                               "alb","faw3","ldl_c","serum_tg"),quiet=FALSE){
  if(!quiet){
    cat("|| Preparing data ... \n")
  }
  
  if(length(which(colnames(dat)=="faw6_faw3"))==0){
    dat$faw6_faw3<-dat$faw6/dat$faw3
  }
  MEANS=colMedians(as.matrix(dat[,featID]),na.rm=T)
  names(MEANS)<-featID
  SDS=colSds(as.matrix(dat[,featID]),na.rm=T)
  names(SDS)<-featID
  # 1. Subset samples on SD:
  dat <- subset_samples_sd_surrogates(x=as.matrix(dat[,featID]), MEAN=MEANS, SD=SDS,N=4, quiet=quiet)
  
  ## 2. Scale:
  dat[,featID] <- scale(dat[,featID],center=TRUE,scale=TRUE)
  if(!quiet){
    cat("| Perform log transform & scaling to zero mean and unity sd .. ")
  }
  
  if(!quiet){
    cat("Done!\n")
  }
  return(dat)
}

#' Function to compute COVID score by Nightingale Health from the metabolite data.
#'
#' @param dat The NH-metabolomics matrix 
#' @param betas The betas of the linear regression composing the score
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @return The COVID score by Nightingale Health
#' @export
#' 
comp_covid_score <- function(dat,betas=covid_betas ,quiet=FALSE){
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
      covidScore<-as.data.frame(as.matrix(prepped_dat[,betas$Abbreviation]) %*% betas$Beta_value)
      colnames(covidScore)<-"covidScore"
      rownames(covidScore)<-rownames(prepped_dat)
      
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
  
  if(!quiet){
    cat("Done!\n\n")
  }
  return(covidScore)
}


#########################
## MetaboAge functions ##
#########################
#' Function to prepare the NH metabolomics dataset to compute MetaboAge by van den Akker et al.
#'
#' @param mat The NH-metabolomics matrix; it may contain a mixture of flags and metabolites in the columns.
#' @param PARAM is a list holding the parameters to compute the metaboAge
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @param Nmax_zero max number of zeros allowed per sample
#' @param Nmax_miss max number of missing values allowed per sample
#' @return The NH-metabolomics matrix after checking for zeros, missing values, samples>5SD from the BBMRI-mean, imputing the missing values and zscale metabolites concentrations
#' @export
#' 
QCprep<-function(mat,PARAM,quiet=TRUE, Nmax_zero=1, Nmax_miss=1){
  # 0. Start:
  if(!quiet){
    cat(report.dim(mat,header="Start"))
  }
  # 1. Subset required metabolites:
  mat <- subset_metabolites_overlap(mat,metabos=PARAM$MET,quiet=quiet)
  # 2. Subset samples on missingness:
  mat <- subset_samples_miss(mat, Nmax=Nmax_miss, quiet=quiet)
  # 3. Subset samples on zeros:
  mat <- subset_samples_zero(mat, Nmax=Nmax_zero, quiet=quiet)
  # 4. Subset samples on SD:
  mat <- subset_samples_sd(as.matrix(mat), MEAN=PARAM$logMEAN, SD=PARAM$logSD, quiet=quiet)
  # 5. Perform scaling:
  mat <- apply.scale(mat, MEAN=PARAM$MEAN, SD=PARAM$SD)
  if(!quiet){
    cat("| Performing scaling ... ")
    cat(" DONE!\n")
  }
  # 6. Perform imputation:
  mat <- impute.miss(mat)
  if(!quiet){
    cat("| Imputation ... ")
    cat(" DONE!\n")
  }
  return(mat)
}

#' Function that apply the MetaboAge model to the NH-metabolomics concentrations
#'
#' @param mat The NH-metabolomics matrix 
#' @param FIT The betas of the linear regression composing the MetaboAge by van den Akker et al.
#' @return The MetaboAge by van den Akker et al.
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

#' Function that subsets the NH-metabolomics matrix features to the selection in metabos
#'
#' @param x The NH-metabolomics matrix 
#' @param metabos String vector containing the metabolites concentration names to be selected
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @return The NH-metabolomics matrix with the selected metabolites
#' @export
#'
subset_metabolites_overlap<-function(x,metabos,quiet=FALSE){
  x <- x[,intersect(colnames(x),metabos),drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header="Selecting metabolites"))
  }
  return(invisible(x))
}

#' Function that subsets the NH-metabolomics matrix samples to the ones with less than Nmax missing values
#'
#' @param x The NH-metabolomics matrix 
#' @param Nmax integer indicating  the max number of missing values allowed per sample
#' @param quiet TRUE/FALSE. if TRUE if will suppress all messages from the function
#' @return The NH-metabolomics matrix with the selected samples
#' @export
#'
subset_samples_miss<-function(x,Nmax=1,quiet=FALSE){
  MISS <- colSums(is.na(t(x)))
  x <- x[which(MISS<=Nmax),,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on missing values [Nmax>=",Nmax,"]")))
  }
  return(invisible(x))
}

#' Function that subsets the NH-metabolomics matrix samples to the ones with less than Nmax zeros
#'
#' @param x The NH-metabolomics matrix 
#' @param Nmax integer indicating  the max number of zeros allowed per sample
#' @param quiet TRUE/FALSE. if TRUE if will suppress all messages from the function
#' @return The NH-metabolomics matrix with the selected samples
#' @export
#'
subset_samples_zero<-function(x,Nmax=1,quiet=FALSE){
  ZERO <- colSums(t(x==0),na.rm=TRUE)
  x <- x[which(ZERO<=Nmax),,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on zero values [Nmax>=",Nmax,"]")))
  }
  return(invisible(x))
}

#' Function created that subsets the NH-metabolomics matrix samples to the ones for which the metabolites
#' included in MetaboAge for which the log of the metabolic concentrations are not more than 5SD away from their mean
#'
#' @param x The NH-metabolomics matrix 
#' @param MEAN numeric vector indicating the mean of the metabolites in x
#' @param SD numeric vector indicating the standard deviations of the metabolites in x
#' @param quiet TRUE/FALSE. if TRUE if will suppress all messages from the function
#' @return The NH-metabolomics matrix with the selected samples
#' @export
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

#' Function that serves at imputing the missing values with zeros
#'
#' @param x The matrix 
#' @return The matrix with missing values imputed to zero
#' @export
#'
impute.miss<-function(x){
  ## This is an boiler-plate solution :)
  x[which(is.na(x))] <- 0
  return(x)
}


#' Function created to scale the NH-metabolomics matrix samples
#'
#' @param dat The matrix 
#' @param MEAN numeric vector indicating the mean of the metabolites present in dat
#' @param SD numeric vector indicating the standard deviations of the metabolites present in dat
#' @param quiet TRUE/FALSE. if TRUE if will suppress all messages from the function
#' @return The matrix z-scaled using the means and sds given
#' @export
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

#' Function to report on the console the dimension of the NH metabolomics matrix
#'
#' @param x The NH metabolomics matrix 
#' @param header string describing the subsampling of the NH-metabolomics matrix
#' @param trailing number of digits to include
#' @return The report of the NH-metabolomics matrix dimension
#' @export
#'
report.dim<-function(x,header,trailing="0"){
  return(paste0(sprintf(paste0("%-",trailing,"s"),paste0("| ",header,":\n")),sprintf("%4s",ncol(x))," metabolites x ",sprintf("%4s",nrow(x))," samples \n"))
}
######################
## Surrogate models ##
######################
#' Function created to calculate:
#' the BMI using height and weight;
#' the LDL cholesterol using HDL cholesterol, triglycerides, totchol;
#' eGFR using creatinine, sex and age.
#'
#' @param phenotypes data.frame containing height and weight, HDL cholesterol, triglycerides, totchol, sex and age
#' @param metabo_measures NH-metabolomics matrix containing creatinine (crea)
#' @return The phenotypes data.frame with the addition of BMI, LDL cholesterol and eGFR
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

#' Function created to binarize all the phenotypes supplied using the threshold defined 
#' in the surrogates by Bizzarri et al.
#'
#' @param data phenotypes data.frame containing some of the following variables (with the same namings):
#' "sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking",
#' "metabolic_syndrome", "alcohol_consumption", "age","BMI", "ln_hscrp","waist_circumference",
#' "weight","height", "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb" 
#' @return The phenotypic variables binarized following the thresholds in in the surrogates by Bizzarri et al.
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


#' Function created to visualize the binarized variables in a barplot following the 
#' rules described in Bizzarri et al.
#'
#' @param bin_phenotypes the binarized variables
#' @return A plot with a barplot of the binarized variables
#' @export
#'
pheno_barplots<-function(bin_phenotypes){
  bin_phenotypes<-data.frame(ID=rownames(bin_phenotypes), bin_phenotypes)
  molten_pheno <- reshape2::melt(bin_phenotypes, id.vars= "ID", value.name="Value", variable.name="Variable", na.rm=F)
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

#' Function to prepare the metabolomics dataset to compute the surrogates by Bizzarri et al.
#'
#' @param mat The NH-metabolomics matrix; it may contain a mixture of flags and metabolites in the columns.
#' @param PARAM_surrogates is a list holding the parameters to compute the surrogates
#' @param Nmax_zero max number of zeros allowed per sample
#' @param Nmax_miss max number of missing values allowed per sample
#' @param quiet TRUE/FALSE if TRUE if will suppress all messages from the function
#' @return The NH-metabolomics matrix after checking for zeros, missing values, samples>5SD from the BBMRI-mean, 
#' imputing the missing values and zscale metabolites concentrations
#' @export
#' 
QCprep_surrogates<-function(mat,PARAM_surrogates, 
                             Nmax_miss=1, Nmax_zero=1,quiet=FALSE){
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
  mat <- impute.miss(mat)
  if(!quiet){
  cat("| Imputation ... ")
  cat(" DONE!\n")
  }
  return(mat)
}

#' Function created to calculate the surrogate scores by Bizzarri et al. from  the NH-metabolomics matrix
#'
#' @param met The NH-metabolomics matrix 
#' @param pheno The phenotypic matrix including this clinical variables: "sex","diabetes", "lipidmed",  "blood_pressure_lowering_med", "current_smoking",
#' "metabolic_syndrome", "alcohol_consumption", "age","BMI", "ln_hscrp","waist_circumference",
#' "weight","height", "triglycerides", "ldl_chol", "hdlchol", "totchol", "eGFR","wbc","hgb" 
#' @param PARAM_surrogates is a list holding the parameters to compute the surrogates
#' @param bin_names string vector indicating the names of the binary variables obtained
#' @param Nmax_zero max number of zeros allowed per sample
#' @param Nmax_miss max number of missing values allowed per sample
#' @param post TRUE/FALSE. If TRUE it will calculate the posterior probabilities
#' @param roc TRUE/FALSE. If TRUE it will produce an image with all the roc curves 
#' for the surrogate with the avialable phenotypes in the current dataset
#' @param quiet TRUE/FALSE. if TRUE if will suppress all messages from the function
#' @return if pheno is not available: a list with the surrogates and the NH metabolomics matrix after QC.
#' if pheno is available: a list with the surrogates, the roc curves, the phenotypes, the binarized phenotypes and the NH metabolomics matrix after QC,
#' @export
#
calculate_surrogate_scores <- function(met, pheno, PARAM_surrogates, bin_names, Nmax_miss=1,Nmax_zero=1, post=TRUE, roc=FALSE,quiet=FALSE){
  bin_surro=paste0("s_",bin_names)
  #QC of the metabolites
  metabo_measures<-QCprep_surrogates(as.matrix(met[,MET63]), PARAM_surrogates, quiet=quiet,Nmax_miss=Nmax_miss,Nmax_zero=Nmax_zero)
  
  if(roc){
    phenotypes<-pheno[which(rownames(pheno)%in%rownames(metabo_measures)),]
    #Calculating the binarized surrogates
    bin_pheno<-binarize_all_pheno(phenotypes)
    all(rownames(bin_pheno)==rownames(metabo_measures))
    
    #Surrogates Calculation
    surrogates<-foreach::foreach(i=bin_names, .combine="cbind") %do% {
      pred<-apply.fit_surro(as.matrix(metabo_measures),PARAM_surrogates$models_betas[paste0("s_",i),])
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

#' Function created that subsets the NH-metabolomics matrix samples to the ones for which the metabolites
#' included in the surrogates calculations for which the metabolic concentrations are not more than 5SD away from their mean in BBMRI-NL
#'
#' @param x The NH-metabolomics matrix 
#' @param MEAN numeric vector indicating the mean of the metabolites in x
#' @param SD numeric vector indicating the standard deviations of the metabolites in x
#' @param N number of standard deviations away from the mean
#' @param quiet TRUE/FALSE. if TRUE if will suppress all messages from the function
#' @return The NH-metabolomics matrix with the selected samples
#' @export
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

#' Function that apply on of the surrogates models to the NH-metabolomics concentrations
#'
#' @param mat The NH-metabolomics matrix 
#' @param FIT The betas of the logistic regressions composing the surrogates by Bizzarri et al.
#' @return The surrogates by Bizzarri et al.
#' @export
#'
apply.fit_surro<-function(mat,FIT, post=TRUE){
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

#' Function that calculates one of the surrogate and creates figures with their accuracy
#'
#' @param FIT The betas of the logistic regressions composing the surrogates by Bizzarri et al.
#' @param data The NH-metabolomics matrix 
#' @param title_img string with title
#' @param plot TRUE/FALSE. If TRUE it calculates the roc curve
#' @return If plot==T The surrogate predictions and the roc curve. If plot==F only the surrogate predictions
#' @export
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

####################
## Plot functions ##
####################
#' Function for plotting a heatmap indcating missing & zero values
#'
#' @param dat The matrix or data.frame
#' @return Plot with an heatmap and two histogram on the sides indicating zeros and missing values
#' @export
#'
plot_na_heatmap  <- function(dat){
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

#' Function that creates a plotly image of the selected surrogates
#'
#' @param surrogates The data.frame containing the surrogate values by Bizzarri et al.
#' @param bin_phenotypes data.frame with the binarized phenotypes
#' @param x_name string with the names of the selected binary phenotypes for the roc
#' @return plotly image with the ROC curves for the selected variables
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
  
  ROC_curve<-pROC::roc(bin_phenotypes[,selected[1]], as.numeric(surrogates[,paste0("s_",selected[1])]), plot=F, col=c21[surro], quiet = TRUE,
                 lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
  
  roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
  pl <- plotly::plot_ly(roc_df, x= ~Specificity, y= ~Sensitivity, 
                type = 'scatter', mode = 'lines', line = list(color= c21[surro],width = 4),
                name=paste0(surro, ":\nAUC=",round(as.numeric(ROC_curve$auc),digits = 3))) %>%
    plotly::add_segments(x = 1, xend =0, y = 0, yend = 1, line = list(dash = "dash", width= 1, color="black"),
                 name="random:\nAUC=0.5")
  
  if(length(selected)>1){
    surro_images<-foreach::foreach(i=2:length(selected), .combine="cbind") %do% {
      surro<-paste0("s_",selected[i])
  
      ROC_curve<-pROC::roc(bin_phenotypes[,selected[i]], as.numeric(surrogates[,paste0("s_",selected[i])]), plot=F, col=c21[surro], quiet = TRUE,
                     lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
      roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
      
      pl<-pl %>% 
        plotly::add_trace(x= roc_df$Specificity, y= roc_df$Sensitivity, 
                           type = 'scatter', mode = 'lines', line = list(color= c21[surro],width = 4),
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

#' Function that creates a plotly image with the subplots with ROCs of the surrogates of the available variables 
#'
#' @param surrogates data.frame containing the surrogate values by Bizzarri et al.
#' @param bin_phenotypes data.frame with the binarized phenotypes output of binarize_all_pheno
#' @return plotly image with all the ROCs for the available clinical variables
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
    ROC_curve<-pROC::roc(bin_phenotypes[,i], as.numeric(surrogates[,surro]), plot=F, col=c21[surro], quiet = TRUE,
                   lwd=4, print.auc=TRUE, main = i, xaxs="i", yaxs="i", grid=TRUE, asp=NA) #AUC
    roc_df<-data.frame(Specificity=ROC_curve$specificities, Sensitivity=ROC_curve$sensitivities)
    
    pl <- plotly::plot_ly(roc_df, x= ~Specificity, y= ~Sensitivity, 
                  type = 'scatter', mode = 'lines', line = list(color= c21[surro],width = 4)) %>%
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

#' Function that creates a plotly image of the selected surrogates
#'
#' @param surrogates The data.frame containing the surrogate values by Bizzarri et al.
#' @param bin_phenotypes data.frame with the binarized phenotypes
#' @return plotly image with the ROC curves for the selected variables
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
  pl_surro<-plotly::plot_ly(Surrogates,
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
           )
  
  return(pl_surro)
}

#' Function created to visualize the accuracies in the current dataset compared to the
#' accuracies in the Leave One Biobank Out Validation in Bizzarri et al.
#'
#' @param surrogates matrix with the surrogate variables calculated in Bizzarri et all.
#' @param bin_phenotypes magtrix with the binarized variables
#' @param bin_pheno_available character vector with the phenotypes availables
#' @param acc_LOBOV accuracy of LOBOV calculated in Bizzarri et al.
#' @return Boxplot with the accuracies of the LOBOV
#' @export
#'
LOBOV_accuracies<-function(surrogates, bin_phenotypes, bin_pheno_available, acc_LOBOV){
  AUCs<-foreach(i=1:length(bin_pheno_available), .combine="rbind") %do%{
    ROC_curve<-pROC::roc(bin_phenotypes[,bin_pheno_available[i]], as.numeric(surrogates[,paste0("s_",bin_pheno_available[i])]), plot=F, 
                         col=c21[surro], quiet = TRUE,
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

#' Function that creates a scatterplot with the real values and predicted values
#'
#' @param x numeric vector with the real values
#' @param p numeric vector with the predicted values
#' @param title string with the title
#' @param xname string with the name of the variable on the x axis
#' @param yname string with the name of the variable on the y axis
#' @return plotly image with the scatterplot of a real vs predicted variable
#' @export
#'
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
                                x = 1.24, y = 0.25,
                                showarrow=FALSE,
                                bordercolor=c("black"),
                                borderwidth=2
                                )
             )
    return(cPlot)
}


# metabo_barplots<-function(metabo_measures,color){
#   par(mfrow = c(3, 4))
#   for (i in 1:56) { # Loop over loop.vector
#     hist(na.omit(metabo_measures[,i]), # histogram
#          border="black",
#          prob = TRUE, # show densities instead of frequencies
#          breaks = 1000,
#          xlab="",
#          main = colnames(metabo_measures)[i])
#     lines(density(na.omit(metabo_measures[,i])), # density plot
#           lwd = 3, # thickness of line
#           col = color)
#   }
#   
# }


#' Function to plot the histograms for all the variables in dat
#'
#' @param dat data.frame or matrix with the variables to plot
#' @param color colors selected for all the variables
#' @param scaled TRUE/FALSE. If TRUE the variables will be z-scaled
#' @return plotly image with the histograms for all the variables in dat
#' @export
#'
multi_hist<-function(dat, color=c21, scaled=FALSE){
  
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

#' Function to plot the histogram of the selected variables
#'
#' @param dat data.frame or matrix with the variables to plot
#' @param x_name string with the names of the selected variables in dat
#' @param color colors selected for all the variables
#' @param scaled TRUE/FALSE. If TRUE the variables will be z-scaled before plotting them
#' @param datatype a character vector indicating what data type is beeing plotted
#' @param main title of the plot
#' @return plotly image with the histograms of the selected variables
#' @export
#'
hist_plots<-function(dat, x_name, color=c21, scaled=FALSE, datatype="metabolic score", main="Predictors Distributions"){
  dat<-as.matrix(dat[,x_name])
  colnames(dat)<-x_name
  if(scaled){
    if(length(x_name)>1){
      sd<-matrixStats::colSds(dat,na.rm = T)
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
  plot<-plotly::plot_ly(x = ~dat[,x_name[1]],opacity = 0.45, type = "histogram", histnorm = "probability", marker = list(color = color[x_name[1]]), name = x_name[1])
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
  
#' Function to plot the 57 selected metabolites to their distributions in BBMRI-nl
#'
#' @param dat data.frame or matrix with the metabolites
#' @param x_name string with the name of the selected variable
#' @param color colors selected for all the variables
#' @param scaled TRUE/FALSE. If TRUE the variables will be z-scaled before plotting them
#' @param datatype a character vector indicating what data type is beeing plotted
#' @param main title of the plot
#' @return plotly image with the histogram of the selected variable compared to the distributions in BBMRI-nl
#' @export
#'
BBMRI_hist_plot<-function(dat, x_name, color=c21, scaled=FALSE, datatype="metabolite", main="Comparison with the metabolites measures in BBMRI"){
  dat<-as.matrix(dat[,x_name])
  colnames(dat)<-x_name
  if(scaled){
    dat <-as.matrix((dat-mean(dat,na.rm=T))/stats::sd(dat,na.rm=T))
    colnames(dat)<-x_name
    BBMRI<-BBMRI_hist_scaled[[x_name]]
  }else{
    BBMRI<-BBMRI_hist[[x_name]]
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


#' Function to plot the scaled coefficients of the metabolic scores
#'
#' @param mort_betas dataframe withthe coefficients of the mortality score
#' @param metaboAge_betas dataframe with the coefficients of the metaboAge
#' @param surrogates_betas dataframe with the coefficients of the surrogates
#' @param Ahola_Olli_betas dataframe with the coefficients of the T2D score
#' @param CVD_score_betas dataframe with the coefficients of the CVD score
#' @param COVID_score_betas ataframe with the coefficients of the COVID_score
#' @return heatmapply with the scaled coefficients of the metabolic scores
#' @export
#'
model_coeff_heat<-function(mort_betas, metaboAge_betas, surrogates_betas, Ahola_Olli_betas, CVD_score_betas, COVID_score_betas){
  models_betas<-matrix(0, 24, 62)
  rownames(models_betas)<-c("mortScore","MetaboAge", rownames(surrogates_betas),"T2D-score","CVD_score", "COVID_score")
  colnames(models_betas)<-MET62
  models_betas["mortScore",mort_betas$Abbreviation]<-mort_betas$Beta_value
  models_betas["MetaboAge", names(metaboAge_betas)[-1]]<-metaboAge_betas[-1]
  #surrogates
  for(x in rownames(surrogates_betas)){
    models_betas[x,colnames(surrogates_betas)[-1]]<-surrogates_betas[x,-1]
  }
  met<-Ahola_Olli_betas$Abbreviation[which(Ahola_Olli_betas$Abbreviation %in% MET62)]
  models_betas["T2D-score",met]<-Ahola_Olli_betas$Beta_value[which(Ahola_Olli_betas$Abbreviation %in% MET62)]
  met<-CVD_score_betas$Abbreviation[which(CVD_score_betas$Abbreviation %in% MET62)]
  models_betas["CVD_score",met]<-CVD_score_betas$Beta_value[which(CVD_score_betas$Abbreviation %in% MET62)]
  met<-COVID_score_betas$Abbreviation[which(COVID_score_betas$Abbreviation %in% MET62)]
  models_betas["COVID_score",met]<-COVID_score_betas$Beta_value[which(COVID_score_betas$Abbreviation %in% MET62)]
  
  
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

#' Function to plot the histogram of the mortality score separated for different age ranges
#'
#' @param mort_score data.frame of the mortality score
#' @param phenotypes data.frame with age included in it
#' @return plotly image with the histogram of the mortlity score separated in 3 age ranges
#' @export
#'
hist_plots_mortality<-function(mort_score,phenotypes){
  age<-data.frame(age=phenotypes[rownames(mort_score),"age"])
  age[which(age$age<45),"age_ranges"]<-"age<45"
  age[intersect(which(age$age>=45),which(age$age<65)),"age_ranges"]<-"45>=age<65"
  age[which(age$age>=65),"age_ranges"]<-"age>=65"
  
  mort_score$age_ranges<-age$age_ranges
  colnames(mort_score)[2]<-"age_ranges"
  mort_score$age_ranges<-factor(mort_score$age_ranges, levels=c("age<45","45>=age<65","age>=65"))
  
  mu<-plyr::ddply(mort_score,~age_ranges,dplyr::summarise,mean=mean(mortScore,na.rm=T))
  
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

#' Function that creates a boxplot with a continuous variable split using the binary variable
#'
#' @param data The data.frame containing the 2 variables
#' @param pred character indicating the y variable
#' @param pheno character indicating the binary variable
#' @return plotly boxplot with the continuous variable split using the binary variable
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

#' Function that creates a Kaplan Meier for an Event of the predictor divided by its mean 
#'
#' @param predictors The data.frame containing the predictors
#' @param pheno The data.frame containing the phenotypes
#' @param score a character string indicating which predictor to use
#' @return plotly with a Kaplan Meier for an Event of the predictor divided by its mean 
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
    findInterval( score, c(-Inf,
                         quantile(score, probs=c(0.3333333, 0.6666666)), Inf) ), 
    labels=c("Tertile 1","Tertile 2","Tertile 3")
  ))
  
  # M<-mean(predictors[,score],na.rm=T)
# dat$score_divided<-rep(NA,dim(dat)[1])
#   dat$score_divided[which(dat[,score]<=M)]<-paste0(score,"<mean(",score,")")
#   dat$score_divided[which(dat[,score]>M)]<-paste0(score,">mean(",score,")") 
#   
#   dat$score_divided<-as.factor(dat$score_divided)
  
  dat<-dat[-which(dat$tertile=="Tertile 2"),]
  
  km_trt_fit <- survminer::surv_fit(Surv(time = dat$EventAge-dat$age , event=dat$Event) ~ 
                                      tertile, data=dat)
  
  pl<-ggplot2::autoplot(km_trt_fit, xlab = "Time (in years)", ylab = "Survival",
                        main = paste0("Kaplan Meier of ",score," associated with ", Eventname,
                                     "\np=",formatC(survminer::surv_pvalue(km_trt_fit)[,2], format = "e", digits = 3),
                                     ", N=", dim(dat)[1],", Nevent=",sum(dat$Event)))
  
  
  pl<-plotly::ggplotly(pl) 
  
  
  return(pl)
}


#' SubFunction for plotting data1 x data2 associations on basis of (partial) correlations:
#'
#' @param dat1 matrix 1
#' @param dat2 matrix 2
#' @param feat1 string with the names of the selected variables in dat
#' @param feat2 string with the names of the selected variables in dat2
#' @param method indicates which methods of the correlation to use (look at cor.test)
#' @param quiet TRUE/FALSE. if TRUE if will suppress all messages from the function
#' @return correlations of the selected variables in the 2 martrices
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


#' Function for plotting dat1 x dat2 associations on basis of (partial) correlations:
#'
#' @param res associatins obtained with cor.assoc
#' @param main title of the plot
#' @param zlim max association to plot
#' @param reorder.x TRUE/FALSE. If TRUE it will reorder the x axis based on clustering
#' @param reorder.y TRUE/FALSE. If TRUE it will reorder the y axis based on clustering
#' @param resort_on_p TRUE/FALSE. If TRUE it will reorder the x and y axis based on the pvalues of the associations
#' @param abs TRUE/FALSE. If TRUE it will reorder based the absolute values
#' @param cor.abs TRUE/FALSE. If TRUE it will plot the absolute values 
#' @param text TRUE/FALSE. If TRUE it will include text in the plot
#' @param value TRUE/FALSE. If TRUE it will plot the values of the associations
#' @param cex.num integer indicating the size of the text
#' @param fontsize integer indicating the size of the title and axis
#' @return heatmap with the results of cor.assoc
#' @export
#'
plot_cor_assoc <- function(res, main=NULL,zlim=NULL,reorder.x=FALSE,reorder.y=reorder.x,resort_on_p=FALSE,
                           abs=FALSE,cor.abs=FALSE,text=TRUE,value=FALSE, cex.num=8, fontsize=8){

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
      IND2 <- which(stats::p.adjust(p,method="fdr")<=0.05) #it compares the p-values adjusted with the fdr with thresh=0.05
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
    pheatmap::pheatmap(s, scale = "none", main=main, display_numbers = textMatrix,
                 cluster_rows = F, cluster_cols = F, margin = c(5,5), fontsize_number = cex.num, cex=1.2,
                 show_rownames=1, show_colnames=1, colors = grDevices::colorRampPalette(c("blue","white", "red"))(20))
  }else{
    pheatmap::pheatmap(s, scale = "none", main=main,
             cluster_rows = F, cluster_cols = F, margin = c(5,5), fontsize_number = cex.num, cex=1.2,
             fontsize=fontsize,
             show_rownames=1, show_colnames=1, colors = grDevices::colorRampPalette(c("blue","white", "red"))(20))
  }
}


#' Function creating a plotly image with the dat1 x dat2 associations on basis of (partial) correlations:
#'
#' @param res associatins obtained with cor.assoc
#' @param main title of the plot
#' @param zlim max association to plot
#' @param reorder.x TRUE/FALSE. If TRUE it will reorder the x axis based on clustering
#' @param reorder.y TRUE/FALSE. If TRUE it will reorder the y axis based on clustering
#' @param resort_on_p TRUE/FALSE. If TRUE it will reorder the x and y axis based on the pvalues of the associations
#' @param abs TRUE/FALSE. If TRUE it will reorder based the absolute values
#' @param cor.abs TRUE/FALSE. If TRUE it will plot the absolute values 
#' @param reorder_dend TRUE/FALSE. If TRUE it will use the dendrogram calculated by plotly
#' @return heatmap with the results of cor.assoc
#' @export
#'
plot_corply <- function(res,main=NULL,zlim=NULL,reorder.x=FALSE,reorder.y=reorder.x,resort_on_p=FALSE,
                           abs=FALSE,cor.abs=FALSE,reorder_dend=F){

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

#' helper function for clustering cor.assoc results reordering based on the associations
#'
#' @param res Results of cor.assoc
#' @param abs TRUE/FALSE. If TRUE it will cluster the absolute values
#' @return Returns the clustered order of the associations
#' @export
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

#' helper function for clustering cor.assoc results reordering based on the pvalues
#'
#' @param res Results of cor.assoc
#' @param abs TRUE/FALSE. If TRUE it will cluster the absolute values
#' @return Returns the clustered order of the associations based on the pvalues
#' @export
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

#' helper function for extracting statistics from cor.assoc results
#'
#' @param res Results of cor.assoc
#' @return the matrix of associations
#' @export
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

#' helper function for extracting pvalues from cor.assoc results
#'
#' @param res Results of cor.assoc
#' @return the matrix of the pvalues of the associations
#' @export
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


## Defining a helper function to check whether a supplied matrix is symmetric:
#' helper function to check whether a supplied matrix is symmetric
#'
#' @param res Results of cor.assoc
#' @return TRUE/FALSE
#' @export
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

#' Function that calculates the Platt Calibrations
#'
#' @param r.calib binary real data
#' @param p.calib predicted probabilities
#' @param nbins number of bins to create the plots
#' @param pl TRUE/FALSE. If TRUE it will plot the Reliability diagram and histogram of the calibrations
#' @return list with samples, responses, calibrations, ECE, MCE and calibration plots if save==T
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
  raw.logl <- (-1/length(resp)) * sum(resp * log(pred) + 
                                        (1 - resp) * (log(1 - pred)), na.rm = TRUE)
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

#' Function that plots the Platt Calibrations using plotly
#'
#' @param r binary real data
#' @param p predicted probabilities
#' @param name character string indicating the name of the calibrated variable
#' @param nbins number of bins to create the plots
#' @param annot_x integer indicating the x axis points in which the ECE and MCE values will be plotted
#' @param annot_y integer indicating the y axis points in which the ECE and MCE values will be plotted
#' @return list with Reliability diagram and histogram with calibrations and original predictions
#' @export
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
  r.calib<-as.numeric(r[train_ind])-1
  p.calib<-p[train_ind]
  resp<-as.numeric(r[-train_ind])-1
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

#' helper function that calculates the Platt Calibrations  for all the surrogates
#'
#' @param bin_phenotypes data.frame with binary phenotypes resulting form binarize_all_pheno
#' @param surrogates data.frame with surrogates resulting from calculate_surrogate_scores
#' @param bin_names string with the names of the binarize clinical variables
#' @param bin_pheno_available string with the names of the binarize clinical variables available in the dataset
#' @param pl TRUE/FALSE. If TRUE creates the calibration plots
#' @param nbins number of bins for the plots
#' @return list with the calibrated surrogates
#' @export
#'
calibration_surro<-function(bin_phenotypes, surrogates, bin_names, bin_pheno_available, pl=FALSE, nbins=10){
  bin_surro<-paste0("s_",bin_names)
  calib<-lapply(1:length(bin_surro), function(i){
    orig<-as.numeric(bin_phenotypes[,bin_names[i]])-1
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

#' helper function that creates a data.frame with the Platt Calibrations
#'
#'@param calibrations list result of calibration_surro
#' @param bin_phenotypes data.frame with binary phenotypes resulting form binarize_all_pheno
#' @param bin_pheno_available string with the names of the binarize clinical variables available in the dataset
#' @return data.frame with the calibrated surrogates
#' @export
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

#' helper function to calculate the ECE of the calibrations
#'
#' @param actual real values of the variables
#' @param predicted predicted values by one of the surrogates
#' @param n_bins the number of bins
#' @return ECE value
#' @export
#'
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
#' @export
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


#' helper function to create a plotly indicating a problem with a plotly image
#'
#' @param main Message to plot
#' @return plotly image
#' @export
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

#' helper function to create a plot indicating a problem with a plotly image
#'
#' @param main Message to plot
#' @return plot
#' @export
#'
NA_message<-function(main="Metabolites are missing, please check your upload!"){
  graphics::plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  graphics::text(x = 0.5, y = 0.5, main,cex = 1.6, col = "black", font=2)
}





#' Function that produces a MetaboWAS
#'
#' @param met metabolites data.frame
#' @param pheno phenotypes data.frame
#' @param test_variable the variable to be investigated 
#' @param covariates the covariates that you want to add
#' @param img TRUE to plot a Manhattan plot
#' @return res= the results of the MetaboWAS, manhplot= the Manhattan plot, N_hits= the number of significant hits
#' @export
#'
MetaboWAS<-function(met, pheno, test_variable, covariates, img=T){
  res<-do.metabowas(phen=pheno, dat= t(met), test_variable = test_variable, covariates=covariates)
  N_hits<-length(which(res$pval.adj<=0.05))
  
  if(img){
    metabo_names_translator<-metabo_names_translator[which(!is.na(metabo_names_translator$BBMRI_names)),]
    res<-res[metabo_names_translator$BBMRI_names,]
    res$groups<-metabo_names_translator$group
    res$ord<-1:dim(res)[1]
    res$met_name<-rownames(res)
    
    axis_set <- res %>% 
      group_by(groups) %>% 
      summarise(center = mean(ord))
    
    ylim <- res %>% 
      filter(pval.adj == min(pval.adj,na.rm = T)) %>% 
      mutate(ylim = abs(floor(log10(pval.adj))) + 2) %>% 
      pull(ylim)
    if(!is.null(covariates)){
      if(length(covariates)>1){
        title<-paste("MetaboWAS:",test_variable,"corrected for",paste(covariates, collapse=","),
                     "\nNumber of significant hits:",N_hits)
      }else{
        title<-paste("MetaboWAS:",test_variable,"corrected for",covariates,
                     "\nNumber of significant hits:",N_hits)
      }
      
    }else{
      title<-paste("MetaboWAS:",test_variable,
                   "\nNumber of significant hits:",N_hits)
    }
    
    manhplot <- ggplot(res, aes(label = met_name, label2=groups, label3=pval.adj)) +
      geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = "dashed") + 
      geom_point(aes(x = ord, y = -log10(pval.adj),
                     color = as.factor(groups), size = -log10(pval.adj)), alpha = 0.75) +
      ggtitle(title)+
      scale_x_continuous(label = axis_set$groups, breaks = axis_set$center) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "metabolites", 
           y = "-log(p adjusted)") + 
      theme_minimal() +
      theme( 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
      )
    
    manhplot<-ggplotly(manhplot, 
                       tooltip=c("met_name", "groups", "pval.adj"))
  }
  
  list(res=res, manhplot=manhplot, N_hits=N_hits)
}


#' Helper function to do a MetaboWAS
#'
#' @param phen phenotypes data.frame
#' @param dat metabolites data.frame
#' @param test_variable the variable to be investigated 
#' @param covariates the covariates that you want to add
#' @param quiet if FALSE it will plot the amount of people avaialble
#' @param PC the number of Principal compontents that represent the 90% of the variance, in case you want to correct for this number (instead of Bonferroni/BH)
#' @return results= the results of the MetaboWAS (estimate, tstatistics, pvalue, BH corrected pvalue)
#' @export
#'
do.metabowas<-function(phen,dat,test_variable="age",covariates=c("sex"), PC=NULL, quiet=T){
  if(!is.null(covariates)){
    vars <- phen[, c(test_variable, covariates)]
  }else{
    vars <- as.data.frame(phen[,test_variable])
    rownames(vars)<-rownames(phen)
    colnames(vars)<-test_variable
  }
  vars <- na.omit(vars)
  if(!quiet){
    print(paste("The number of samples is:",dim(vars)[1]))
  }
  dat <- dat[, match(rownames(vars), colnames(dat))]
  design <- model.matrix(~ ., vars)
  fit <- lmFit(dat, design)
  fit <- eBayes(fit)
  result <- limma::topTable(fit, coef=2, number=Inf)
  colnames(result) <- c("estimate","AveEpr","tstat","pval","pval.adj")
  result <- result[,c("estimate","tstat","pval","pval.adj")]
  if(!is.null(PC)){
    result$PC.adj<-(result$pval * PC)
  }
  return(result)
}




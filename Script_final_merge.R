#############################################################################################################################
#                                                                                                                           #
#                           Grey and white matter metrics demonstrate distinct and complementary                            #
#                                   prediction of differences in cognitive performance                                      #
#                                     in children: Findings from ABCD (N= 11 876)                                           #
#                                                                                                                           #
#                                            Michel et al, 2023, __Name of the journal__                                    #
#                                                                                                                           #
#---------------------------------------------------------------------------------------------------------------------------#
#                                                                                                                           #
#                                      URL & DOI:                                                                           #
#                                 Supplementary materials:                                                                  #
#                                                                                                                           #
#---------------------------------------------------------------------------------------------------------------------------#
#                                                                                                                           #
#                                         Script for merging the files                                                      #
#                                                                                                                           #
#############################################################################################################################


# -------------- Table of contents --------------------------------------------------------------------------------------------
# -------------- 1. Load libraries (line 35) ----------------------------------------------------------------------------------
# -------------- 2. Load the datasets and select the relevant variables (line 51) ---------------------------------------------
# -------------- 3. Cleaning the datasets (line 111) --------------------------------------------------------------------------
# -------------- 4. Average every regions of interest bilaterally (line 149) --------------------------------------------------
# -------------- 5. Creation of the merge and processed datasets (line 322) ---------------------------------------------------
# -------------- 6. Descriptives (line 340) -----------------------------------------------------------------------------------


# This script reproduces all the analyses from the paper
# To reproduce results exactly, please request raw data from https://nda.nih.gov/nda/access-data-info.html
# Script author: Lea Michel, 2022, using RStudio (version 4.1.0)

-------------------------------------------------------------------------------------------------------------------------------
#Load libraries
-------------------------------------------------------------------------------------------------------------------------------

library(ggplot2)
library(ggpubr) 
library(dplyr)
library(tidyverse)
library(rstatix)
library(tidyr)
library(plyr)
library(udunits2)
library(ggseg)
library(summarytools)
library(corrplot)

#-------------------------------------------------------------------------------------------------------------------------------
#Load the datasets and select the relevant variables
#-------------------------------------------------------------------------------------------------------------------------------

setwd("C:/Users/leamic/Documents/Project 1/Data/Raw data")
data_Youth_NIH_TB_Summary_Scores<-read.delim(file="abcd_tbss01.txt",header=TRUE)
data_Pearson_Scores<-read.delim(file="abcd_ps01.txt",header=TRUE)
data_Little_Man_Task_Summary_Scores<-read.delim(file="lmtp201.txt",header=TRUE)
brain_sMRI<-read.delim(file="abcd_smrip10201.txt",header=TRUE)
brain_DTI<-read.delim(file="abcd_dti_p101.txt",header=TRUE)

#Deletion of the first row
data_Youth_NIH_TB_Summary_Scores<-data_Youth_NIH_TB_Summary_Scores[-1,]
data_Pearson_Scores<-data_Pearson_Scores[-1,]
data_Little_Man_Task_Summary_Scores<-data_Little_Man_Task_Summary_Scores[-1,]
brain_sMRI<-brain_sMRI[-1,]
brain_DTI<-brain_DTI[-1,]

#Selection of the relevant measures and computing the scores
##Names of the tasks:	nihtbx_picvocab_uncorrected, nihtbx_flanker_uncorrected, nihtbx_reading_uncorrected
data_Youth_NIH_TB_Summary_Scores<-data_Youth_NIH_TB_Summary_Scores[,c("subjectkey","interview_age","eventname","sex","nihtbx_picvocab_uncorrected","nihtbx_flanker_uncorrected","nihtbx_reading_uncorrected")]
##For the Rey Auditory Verbal Learning Task: The total number of items correctly recalled across the five learning trials was summed to produce a measure of auditory verbal learning (pea_ravlt_sd_trial_i_tc, pea_ravlt_sd_trial_ii_tc, pea_ravlt_sd_trial_iii_tc, pea_ravlt_sd_trial_iv_tc, pea_ravlt_sd_trial_v_tc)
data_Pearson_Scores<-data_Pearson_Scores[,c("subjectkey","interview_age","eventname","sex","pea_ravlt_sd_trial_i_tc","pea_ravlt_sd_trial_ii_tc","pea_ravlt_sd_trial_iii_tc","pea_ravlt_sd_trial_iv_tc","pea_ravlt_sd_trial_v_tc")]
data_Pearson_Scores<-data_Pearson_Scores %>%
  mutate(pea_ravlt_sum=as.numeric(pea_ravlt_sd_trial_i_tc)+as.numeric(pea_ravlt_sd_trial_ii_tc)+as.numeric(pea_ravlt_sd_trial_iii_tc)+as.numeric(pea_ravlt_sd_trial_iv_tc)+as.numeric(pea_ravlt_sd_trial_v_tc)) %>%
  select(-c("pea_ravlt_sd_trial_i_tc","pea_ravlt_sd_trial_ii_tc","pea_ravlt_sd_trial_iii_tc","pea_ravlt_sd_trial_iv_tc","pea_ravlt_sd_trial_v_tc"))
##For the Little Man Task: The summary score is the percentage of correct answers (lmt_scr_perc_correct)
data_Little_Man_Task_Summary_Scores<-data_Little_Man_Task_Summary_Scores[,c("subjectkey","interview_age","eventname","sex","lmt_scr_perc_correct")]
data_Little_Man_Task_Summary_Scores<-data_Little_Man_Task_Summary_Scores %>%
  mutate(lmt_scr_perc_correct=as.numeric(lmt_scr_perc_correct)*100)
##Datasets for grey matter and for white matter
data_greymatter<-brain_sMRI[,c("subjectkey","interview_age","eventname","sex","smri_thick_cdk_banksstslh","smri_thick_cdk_cdacatelh","smri_thick_cdk_cdmdfrlh","smri_thick_cdk_cuneuslh","smri_thick_cdk_ehinallh","smri_thick_cdk_fusiformlh","smri_thick_cdk_ifpllh","smri_thick_cdk_iftmlh","smri_thick_cdk_ihcatelh","smri_thick_cdk_locclh","smri_thick_cdk_lobfrlh","smri_thick_cdk_linguallh","smri_thick_cdk_mobfrlh","smri_thick_cdk_mdtmlh","smri_thick_cdk_parahpallh","smri_thick_cdk_paracnlh","smri_thick_cdk_parsopclh","smri_thick_cdk_parsobislh","smri_thick_cdk_parstgrislh","smri_thick_cdk_pericclh","smri_thick_cdk_postcnlh","smri_thick_cdk_ptcatelh","smri_thick_cdk_precnlh","smri_thick_cdk_pclh","smri_thick_cdk_rracatelh","smri_thick_cdk_rrmdfrlh","smri_thick_cdk_sufrlh","smri_thick_cdk_supllh","smri_thick_cdk_sutmlh","smri_thick_cdk_smlh","smri_thick_cdk_frpolelh","smri_thick_cdk_tmpolelh","smri_thick_cdk_trvtmlh","smri_thick_cdk_insulalh","smri_thick_cdk_banksstsrh","smri_thick_cdk_cdacaterh","smri_thick_cdk_cdmdfrrh","smri_thick_cdk_cuneusrh","smri_thick_cdk_ehinalrh","smri_thick_cdk_fusiformrh","smri_thick_cdk_ifplrh","smri_thick_cdk_iftmrh","smri_thick_cdk_ihcaterh","smri_thick_cdk_loccrh","smri_thick_cdk_lobfrrh","smri_thick_cdk_lingualrh","smri_thick_cdk_mobfrrh","smri_thick_cdk_mdtmrh","smri_thick_cdk_parahpalrh","smri_thick_cdk_paracnrh","smri_thick_cdk_parsopcrh","smri_thick_cdk_parsobisrh","smri_thick_cdk_parstgrisrh","smri_thick_cdk_periccrh","smri_thick_cdk_postcnrh","smri_thick_cdk_ptcaterh","smri_thick_cdk_precnrh","smri_thick_cdk_pcrh","smri_thick_cdk_rracaterh","smri_thick_cdk_rrmdfrrh","smri_thick_cdk_sufrrh","smri_thick_cdk_suplrh","smri_thick_cdk_sutmrh","smri_thick_cdk_smrh","smri_thick_cdk_frpolerh","smri_thick_cdk_tmpolerh","smri_thick_cdk_trvtmrh","smri_thick_cdk_insularh",
                               "smri_area_cdk_banksstslh","smri_area_cdk_cdacatelh","smri_area_cdk_cdmdfrlh","smri_area_cdk_cuneuslh","smri_area_cdk_ehinallh","smri_area_cdk_fusiformlh","smri_area_cdk_ifpllh","smri_area_cdk_iftmlh","smri_area_cdk_ihcatelh","smri_area_cdk_locclh","smri_area_cdk_lobfrlh","smri_area_cdk_linguallh","smri_area_cdk_mobfrlh","smri_area_cdk_mdtmlh","smri_area_cdk_parahpallh","smri_area_cdk_paracnlh","smri_area_cdk_parsopclh","smri_area_cdk_parsobislh","smri_area_cdk_parstgrislh","smri_area_cdk_pericclh","smri_area_cdk_postcnlh","smri_area_cdk_ptcatelh","smri_area_cdk_precnlh","smri_area_cdk_pclh","smri_area_cdk_rracatelh","smri_area_cdk_rrmdfrlh","smri_area_cdk_sufrlh","smri_area_cdk_supllh","smri_area_cdk_sutmlh","smri_area_cdk_smlh","smri_area_cdk_frpolelh","smri_area_cdk_tmpolelh","smri_area_cdk_trvtmlh","smri_area_cdk_insulalh","smri_area_cdk_banksstsrh","smri_area_cdk_cdacaterh","smri_area_cdk_cdmdfrrh","smri_area_cdk_cuneusrh","smri_area_cdk_ehinalrh","smri_area_cdk_fusiformrh","smri_area_cdk_ifplrh","smri_area_cdk_iftmrh","smri_area_cdk_ihcaterh","smri_area_cdk_loccrh","smri_area_cdk_lobfrrh","smri_area_cdk_lingualrh","smri_area_cdk_mobfrrh","smri_area_cdk_mdtmrh","smri_area_cdk_parahpalrh","smri_area_cdk_paracnrh","smri_area_cdk_parsopcrh","smri_area_cdk_parsobisrh","smri_area_cdk_parstgrisrh","smri_area_cdk_periccrh","smri_area_cdk_postcnrh","smri_area_cdk_ptcaterh","smri_area_cdk_precnrh","smri_area_cdk_pcrh","smri_area_cdk_rracaterh","smri_area_cdk_rrmdfrrh","smri_area_cdk_sufrrh","smri_area_cdk_suplrh","smri_area_cdk_sutmrh","smri_area_cdk_smrh","smri_area_cdk_frpolerh","smri_area_cdk_tmpolerh","smri_area_cdk_trvtmrh","smri_area_cdk_insularh",
                               "smri_vol_cdk_banksstslh","smri_vol_cdk_cdacatelh","smri_vol_cdk_cdmdfrlh","smri_vol_cdk_cuneuslh","smri_vol_cdk_ehinallh","smri_vol_cdk_fusiformlh","smri_vol_cdk_ifpllh","smri_vol_cdk_iftmlh","smri_vol_cdk_ihcatelh","smri_vol_cdk_locclh","smri_vol_cdk_lobfrlh","smri_vol_cdk_linguallh","smri_vol_cdk_mobfrlh","smri_vol_cdk_mdtmlh","smri_vol_cdk_parahpallh","smri_vol_cdk_paracnlh","smri_vol_cdk_parsopclh","smri_vol_cdk_parsobislh","smri_vol_cdk_parstgrislh","smri_vol_cdk_pericclh","smri_vol_cdk_postcnlh","smri_vol_cdk_ptcatelh","smri_vol_cdk_precnlh","smri_vol_cdk_pclh","smri_vol_cdk_rracatelh","smri_vol_cdk_rrmdfrlh","smri_vol_cdk_sufrlh","smri_vol_cdk_supllh","smri_vol_cdk_sutmlh","smri_vol_cdk_smlh","smri_vol_cdk_frpolelh","smri_vol_cdk_tmpolelh","smri_vol_cdk_trvtmlh","smri_vol_cdk_insulalh","smri_vol_cdk_banksstsrh","smri_vol_cdk_cdacaterh","smri_vol_cdk_cdmdfrrh","smri_vol_cdk_cuneusrh","smri_vol_cdk_ehinalrh","smri_vol_cdk_fusiformrh","smri_vol_cdk_ifplrh","smri_vol_cdk_iftmrh","smri_vol_cdk_ihcaterh","smri_vol_cdk_loccrh","smri_vol_cdk_lobfrrh","smri_vol_cdk_lingualrh","smri_vol_cdk_mobfrrh","smri_vol_cdk_mdtmrh","smri_vol_cdk_parahpalrh","smri_vol_cdk_paracnrh","smri_vol_cdk_parsopcrh","smri_vol_cdk_parsobisrh","smri_vol_cdk_parstgrisrh","smri_vol_cdk_periccrh","smri_vol_cdk_postcnrh","smri_vol_cdk_ptcaterh","smri_vol_cdk_precnrh","smri_vol_cdk_pcrh","smri_vol_cdk_rracaterh","smri_vol_cdk_rrmdfrrh","smri_vol_cdk_sufrrh","smri_vol_cdk_suplrh","smri_vol_cdk_sutmrh","smri_vol_cdk_smrh","smri_vol_cdk_frpolerh","smri_vol_cdk_tmpolerh","smri_vol_cdk_trvtmrh","smri_vol_cdk_insularh")]
data_whitematter<-brain_DTI[,c("subjectkey","interview_age","eventname","sex","dmri_dtifa_fiberat_fxrh","dmri_dtifa_fiberat_fxlh","dmri_dtifa_fiberat_cgcrh","dmri_dtifa_fiberat_cgclh","dmri_dtifa_fiberat_cghrh","dmri_dtifa_fiberat_cghlh","dmri_dtifa_fiberat_cstrh","dmri_dtifa_fiberat_cstlh","dmri_dtifa_fiberat_atrrh","dmri_dtifa_fiberat_atrlh","dmri_dtifa_fiberat_uncrh","dmri_dtifa_fiberat_unclh","dmri_dtifa_fiberat_ilfrh","dmri_dtifa_fiberat_ilflh","dmri_dtifa_fiberat_iforh","dmri_dtifa_fiberat_ifolh","dmri_dtifa_fiberat_fmaj","dmri_dtifa_fiberat_fmin","dmri_dtifa_fiberat_cc","dmri_dtifa_fiberat_slfrh","dmri_dtifa_fiberat_slflh","dmri_dtifa_fiberat_tslfrh","dmri_dtifa_fiberat_tslflh","dmri_dtifa_fiberat_pslfrh","dmri_dtifa_fiberat_pslflh","dmri_dtifa_fiberat_scsrh","dmri_dtifa_fiberat_scslh","dmri_dtifa_fiberat_fscsrh","dmri_dtifa_fiberat_fscslh","dmri_dtifa_fiberat_pscsrh","dmri_dtifa_fiberat_pscslh","dmri_dtifa_fiberat_sifcrh","dmri_dtifa_fiberat_sifclh","dmri_dtifa_fiberat_ifsfcrh","dmri_dtifa_fiberat_ifsfclh","dmri_dtifa_fiberat_fxcutrh","dmri_dtifa_fiberat_fxcutlh",
                               "dmri_dtimd_fiberat_fxrh","dmri_dtimd_fiberat_fxlh","dmri_dtimd_fiberat_cgcrh","dmri_dtimd_fiberat_cgclh","dmri_dtimd_fiberat_cghrh","dmri_dtimd_fiberat_cghlh","dmri_dtimd_fiberat_cstrh","dmri_dtimd_fiberat_cstlh","dmri_dtimd_fiberat_atrrh","dmri_dtimd_fiberat_atrlh","dmri_dtimd_fiberat_uncrh","dmri_dtimd_fiberat_unclh","dmri_dtimd_fiberat_ilfrh","dmri_dtimd_fiberat_ilflh","dmri_dtimd_fiberat_iforh","dmri_dtimd_fiberat_ifolh","dmri_dtimd_fiberat_fmaj","dmri_dtimd_fiberat_fmin","dmri_dtimd_fiberat_cc","dmri_dtimd_fiberat_slfrh","dmri_dtimd_fiberat_slflh","dmri_dtimd_fiberat_tslfrh","dmri_dtimd_fiberat_tslflh","dmri_dtimd_fiberat_pslfrh","dmri_dtimd_fiberat_pslflh","dmri_dtimd_fiberat_scsrh","dmri_dtimd_fiberat_scslh","dmri_dtimd_fiberat_fscsrh","dmri_dtimd_fiberat_fscslh","dmri_dtimd_fiberat_pscsrh","dmri_dtimd_fiberat_pscslh","dmri_dtimd_fiberat_sifcrh","dmri_dtimd_fiberat_sifclh","dmri_dtimd_fiberat_ifsfcrh","dmri_dtimd_fiberat_ifsfclh","dmri_dtimd_fiberat_fxcutrh","dmri_dtimd_fiberat_fxcutlh",
                               "dmri_dtivol_fiberat_fxrh","dmri_dtivol_fiberat_fxlh","dmri_dtivol_fiberat_cgcrh","dmri_dtivol_fiberat_cgclh","dmri_dtivol_fiberat_cghrh","dmri_dtivol_fiberat_cghlh","dmri_dtivol_fiberat_cstrh","dmri_dtivol_fiberat_cstlh","dmri_dtivol_fiberat_atrrh","dmri_dtivol_fiberat_atrlh","dmri_dtivol_fiberat_uncrh","dmri_dtivol_fiberat_unclh","dmri_dtivol_fiberat_ilfrh","dmri_dtivol_fiberat_ilflh","dmri_dtivol_fiberat_iforh","dmri_dtivol_fiberat_ifolh","dmri_dtivol_fiberat_fmaj","dmri_dtivol_fiberat_fmin","dmri_dtivol_fiberat_cc","dmri_dtivol_fiberat_slfrh","dmri_dtivol_fiberat_slflh","dmri_dtivol_fiberat_tslfrh","dmri_dtivol_fiberat_tslflh","dmri_dtivol_fiberat_pslfrh","dmri_dtivol_fiberat_pslflh","dmri_dtivol_fiberat_scsrh","dmri_dtivol_fiberat_scslh","dmri_dtivol_fiberat_fscsrh","dmri_dtivol_fiberat_fscslh","dmri_dtivol_fiberat_pscsrh","dmri_dtivol_fiberat_pscslh","dmri_dtivol_fiberat_sifcrh","dmri_dtivol_fiberat_sifclh","dmri_dtivol_fiberat_ifsfcrh","dmri_dtivol_fiberat_ifsfclh","dmri_dtivol_fiberat_fxcutrh","dmri_dtivol_fiberat_fxcutlh")]

#Merge the five cognitive tasks into one dataframe 
data_cognition<-merge(data_Youth_NIH_TB_Summary_Scores,data_Pearson_Scores,by=c("subjectkey","eventname","interview_age","sex"),all=T)
data_cognition<-merge(data_cognition,data_Little_Man_Task_Summary_Scores,by=c("subjectkey","eventname","interview_age","sex"),all=T)


#Unifying labels for easier syntax
names(data_cognition)<- c("subjectkey","eventname","interview_age","sex","Picture_Vocabulary","Flanker","Oral_Reading_Recognition","Rey_Auditory_Verbal_Learning","Little_Man")
names(data_greymatter) <- c("subjectkey","interview_age","eventname","sex","lh_bankssts_ct","lh_caudalanteriorcingulate_ct","lh_caudalmiddlefrontal_ct","lh_cuneus_ct","lh_entorhinal_ct","lh_fusiform_ct","lh_inferiorparietal_ct","lh_inferiortemporal_ct","lh_isthmuscingulate_ct","lh_lateraloccipital_ct","lh_lateralorbitofrontal_ct","lh_lingual_ct","lh_medialorbitofrontal_ct","lh_middletemporal_ct","lh_parahippocampal_ct","lh_paracentral_ct","lh_parsopercularis_ct","lh_parsorbitalis_ct","lh_parstriangularis_ct","lh_pericalcarine_ct","lh_postcentral_ct","lh_posteriorcingulate_ct","lh_precentral_ct","lh_precuneus_ct","lh_rostralanteriorcingulate_ct","lh_rostralmiddlefrontal_ct","lh_superiorfrontal_ct","lh_superiorparietal_ct","lh_superiortemporal_ct","lh_supramarginal_ct","lh_frontalpole_ct","lh_temporalpole_ct","lh_transversetemporal_ct","lh_insula_ct","rh_bankssts_ct","rh_caudalanteriorcingulate_ct","rh_caudalmiddlefrontal_ct","rh_cuneus_ct","rh_entorhinal_ct","rh_fusiform_ct","rh_inferiorparietal_ct","rh_inferiortemporal_ct","rh_isthmuscingulate_ct","rh_lateraloccipital_ct","rh_lateralorbitofrontal_ct","rh_lingual_ct","rh_medialorbitofrontal_ct","rh_middletemporal_ct","rh_parahippocampal_ct","rh_paracentral_ct","rh_parsopercularis_ct","rh_parsorbitalis_ct","rh_parstriangularis_ct","rh_pericalcarine_ct","rh_postcentral_ct","rh_posteriorcingulate_ct","rh_precentral_ct","rh_precuneus_ct","rh_rostralanteriorcingulate_ct","rh_rostralmiddlefrontal_ct","rh_superiorfrontal_ct","rh_superiorparietal_ct","rh_superiortemporal_ct","rh_supramarginal_ct","rh_frontalpole_ct","rh_temporalpole_ct","rh_transversetemporal_ct","rh_insula_ct",
                            "lh_bankssts_sa","lh_caudalanteriorcingulate_sa","lh_caudalmiddlefrontal_sa","lh_cuneus_sa","lh_entorhinal_sa","lh_fusiform_sa","lh_inferiorparietal_sa","lh_inferiortemporal_sa","lh_isthmuscingulate_sa","lh_lateraloccipital_sa","lh_lateralorbitofrontal_sa","lh_lingual_sa","lh_medialorbitofrontal_sa","lh_middletemporal_sa","lh_parahippocampal_sa","lh_paracentral_sa","lh_parsopercularis_sa","lh_parsorbitalis_sa","lh_parstriangularis_sa","lh_pericalcarine_sa","lh_postcentral_sa","lh_posteriorcingulate_sa","lh_precentral_sa","lh_precuneus_sa","lh_rostralanteriorcingulate_sa","lh_rostralmiddlefrontal_sa","lh_superiorfrontal_sa","lh_superiorparietal_sa","lh_superiortemporal_sa","lh_supramarginal_sa","lh_frontalpole_sa","lh_temporalpole_sa","lh_transversetemporal_sa","lh_insula_sa","rh_bankssts_sa","rh_caudalanteriorcingulate_sa","rh_caudalmiddlefrontal_sa","rh_cuneus_sa","rh_entorhinal_sa","rh_fusiform_sa","rh_inferiorparietal_sa","rh_inferiortemporal_sa","rh_isthmuscingulate_sa","rh_lateraloccipital_sa","rh_lateralorbitofrontal_sa","rh_lingual_sa","rh_medialorbitofrontal_sa","rh_middletemporal_sa","rh_parahippocampal_sa","rh_paracentral_sa","rh_parsopercularis_sa","rh_parsorbitalis_sa","rh_parstriangularis_sa","rh_pericalcarine_sa","rh_postcentral_sa","rh_posteriorcingulate_sa","rh_precentral_sa","rh_precuneus_sa","rh_rostralanteriorcingulate_sa","rh_rostralmiddlefrontal_sa","rh_superiorfrontal_sa","rh_superiorparietal_sa","rh_superiortemporal_sa","rh_supramarginal_sa","rh_frontalpole_sa","rh_temporalpole_sa","rh_transversetemporal_sa","rh_insula_sa",
                            "lh_bankssts_gmv","lh_caudalanteriorcingulate_gmv","lh_caudalmiddlefrontal_gmv","lh_cuneus_gmv","lh_entorhinal_gmv","lh_fusiform_gmv","lh_inferiorparietal_gmv","lh_inferiortemporal_gmv","lh_isthmuscingulate_gmv","lh_lateraloccipital_gmv","lh_lateralorbitofrontal_gmv","lh_lingual_gmv","lh_medialorbitofrontal_gmv","lh_middletemporal_gmv","lh_parahippocampal_gmv","lh_paracentral_gmv","lh_parsopercularis_gmv","lh_parsorbitalis_gmv","lh_parstriangularis_gmv","lh_pericalcarine_gmv","lh_postcentral_gmv","lh_posteriorcingulate_gmv","lh_precentral_gmv","lh_precuneus_gmv","lh_rostralanteriorcingulate_gmv","lh_rostralmiddlefrontal_gmv","lh_superiorfrontal_gmv","lh_superiorparietal_gmv","lh_superiortemporal_gmv","lh_supramarginal_gmv","lh_frontalpole_gmv","lh_temporalpole_gmv","lh_transversetemporal_gmv","lh_insula_gmv","rh_bankssts_gmv","rh_caudalanteriorcingulate_gmv","rh_caudalmiddlefrontal_gmv","rh_cuneus_gmv","rh_entorhinal_gmv","rh_fusiform_gmv","rh_inferiorparietal_gmv","rh_inferiortemporal_gmv","rh_isthmuscingulate_gmv","rh_lateraloccipital_gmv","rh_lateralorbitofrontal_gmv","rh_lingual_gmv","rh_medialorbitofrontal_gmv","rh_middletemporal_gmv","rh_parahippocampal_gmv","rh_paracentral_gmv","rh_parsopercularis_gmv","rh_parsorbitalis_gmv","rh_parstriangularis_gmv","rh_pericalcarine_gmv","rh_postcentral_gmv","rh_posteriorcingulate_gmv","rh_precentral_gmv","rh_precuneus_gmv","rh_rostralanteriorcingulate_gmv","rh_rostralmiddlefrontal_gmv","rh_superiorfrontal_gmv","rh_superiorparietal_gmv","rh_superiortemporal_gmv","rh_supramarginal_gmv","rh_frontalpole_gmv","rh_temporalpole_gmv","rh_transversetemporal_gmv","rh_insula_gmv")
names(data_whitematter) <- c("subjectkey","interview_age","eventname","sex","rh_fornix_fa","lh_fornix_fa","rh_cingulatecingulum_fa","lh_cingulatecingulum_fa","rh_parahippocampalcingulum_fa","lh_parahippocampalcingulum_fa","rh_corticospinalpyramidal_fa","lh_corticospinalpyramidal_fa","rh_anteriorthalamicradiations_fa","lh_anteriorthalamicradiations_fa","rh_uncinate_fa","lh_uncinate_fa","rh_inferiorlongitudinalfasiculus_fa","lh_inferiorlongitudinalfasiculus_fa","rh_inferiorfrontooccipitalfasiculus_fa","lh_inferiorfrontooccipitalfasiculus_fa","forcepsmajor_fa","forcepsminor_fa","corpuscallosum_fa","rh_superiorlongitudinalfasiculus_fa","lh_superiorlongitudinalfasiculus_fa","rh_temporalsuperiorlongitudinalfasiculus_fa","lh_temporalsuperiorlongitudinalfasiculus_fa","rh_parietalsuperiorlongitudinalfasiculus_fa","lh_parietalsuperiorlongitudinalfasiculus_fa","rh_superiorcorticostriate_fa","lh_superiorcorticostriate_fa","rh_superiorcorticostriatefrontalcortex_fa","lh_superiorcorticostriatefrontalcortex_fa","rh_superiorcorticostriateparietalcortex_fa","lh_superiorcorticostriateparietalcortex_fa","rh_striatalinferiorfrontalcortex_fa","lh_striatalinferiorfrontalcortex_fa","rh_inferiorfrontalsuperiorfrontalcortex_fa","lh_inferiorfrontalsuperiorfrontalcortex_fa","rh_fornix_exfimbria_fa","lh_fornix_exfimbria_fa",
                             "rh_fornix_md","lh_fornix_md","rh_cingulatecingulum_md","lh_cingulatecingulum_md","rh_parahippocampalcingulum_md","lh_parahippocampalcingulum_md","rh_corticospinalpyramidal_md","lh_corticospinalpyramidal_md","rh_anteriorthalamicradiations_md","lh_anteriorthalamicradiations_md","rh_uncinate_md","lh_uncinate_md","rh_inferiorlongitudinalfasiculus_md","lh_inferiorlongitudinalfasiculus_md","rh_inferiorfrontooccipitalfasiculus_md","lh_inferiorfrontooccipitalfasiculus_md","forcepsmajor_md","forcepsminor_md","corpuscallosum_md","rh_superiorlongitudinalfasiculus_md","lh_superiorlongitudinalfasiculus_md","rh_temporalsuperiorlongitudinalfasiculus_md","lh_temporalsuperiorlongitudinalfasiculus_md","rh_parietalsuperiorlongitudinalfasiculus_md","lh_parietalsuperiorlongitudinalfasiculus_md","rh_superiorcorticostriate_md","lh_superiorcorticostriate_md","rh_superiorcorticostriatefrontalcortex_md","lh_superiorcorticostriatefrontalcortex_md","rh_superiorcorticostriateparietalcortex_md","lh_superiorcorticostriateparietalcortex_md","rh_striatalinferiorfrontalcortex_md","lh_striatalinferiorfrontalcortex_md","rh_inferiorfrontalsuperiorfrontalcortex_md","lh_inferiorfrontalsuperiorfrontalcortex_md","rh_fornix_exfimbria_md","lh_fornix_exfimbria_md",
                             "rh_fornix_wmv","lh_fornix_wmv","rh_cingulatecingulum_wmv","lh_cingulatecingulum_wmv","rh_parahippocampalcingulum_wmv","lh_parahippocampalcingulum_wmv","rh_corticospinalpyramidal_wmv","lh_corticospinalpyramidal_wmv","rh_anteriorthalamicradiations_wmv","lh_anteriorthalamicradiations_wmv","rh_uncinate_wmv","lh_uncinate_wmv","rh_inferiorlongitudinalfasiculus_wmv","lh_inferiorlongitudinalfasiculus_wmv","rh_inferiorfrontooccipitalfasiculus_wmv","lh_inferiorfrontooccipitalfasiculus_wmv","forcepsmajor_wmv","forcepsminor_wmv","corpuscallosum_wmv","rh_superiorlongitudinalfasiculus_wmv","lh_superiorlongitudinalfasiculus_wmv","rh_temporalsuperiorlongitudinalfasiculus_wmv","lh_temporalsuperiorlongitudinalfasiculus_wmv","rh_parietalsuperiorlongitudinalfasiculus_wmv","lh_parietalsuperiorlongitudinalfasiculus_wmv","rh_superiorcorticostriate_wmv","lh_superiorcorticostriate_wmv","rh_superiorcorticostriatefrontalcortex_wmv","lh_superiorcorticostriatefrontalcortex_wmv","rh_superiorcorticostriateparietalcortex_wmv","lh_superiorcorticostriateparietalcortex_wmv","rh_striatalinferiorfrontalcortex_wmv","lh_striatalinferiorfrontalcortex_wmv","rh_inferiorfrontalsuperiorfrontalcortex_wmv","lh_inferiorfrontalsuperiorfrontalcortex_wmv","rh_fornix_exfimbria_wmv","lh_fornix_exfimbria_wmv")

#Filter the datasets for T0
data_cognition <- data_cognition %>%
  filter(eventname == "baseline_year_1_arm_1") 
data_greymatter <- data_greymatter %>%
  filter(eventname == "baseline_year_1_arm_1") 
data_whitematter <- data_whitematter %>%
  filter(eventname == "baseline_year_1_arm_1") 

-------------------------------------------------------------------------------------------------------------------------------
#Cleaning the datasets
-------------------------------------------------------------------------------------------------------------------------------

#Checking the class of the variables
sapply(data_cognition,class)
sapply(data_greymatter,class)
sapply(data_whitematter,class)
data_cognition$sex<-factor(data_cognition$sex)
data_greymatter$sex<-factor(data_greymatter$sex)
data_whitematter$sex<-factor(data_whitematter$sex)
data_cognition$eventname<-factor(data_cognition$eventname)
data_greymatter$eventname<-factor(data_greymatter$eventname)
data_whitematter$eventname<-factor(data_whitematter$eventname)
levels(data_cognition$eventname) <- c("T0")
levels(data_greymatter$eventname) <- c("T0")
levels(data_whitematter$eventname) <- c("T0")
data_cognition[,c(3,5:9)]<-apply(data_cognition[,c(3,5:9)], 2,function(x) as.numeric(as.character(x)))
data_greymatter[,c(2,5:208)]<-apply(data_greymatter[,c(2,5:208)], 2,function(x) as.numeric(as.character(x)))
data_whitematter[,c(2,5:115)]<-apply(data_whitematter[,c(2,5:115)], 2,function(x) as.numeric(as.character(x)))

#Remove outliers
boxplot(data_cognition$Picture_Vocabulary)
boxplot(data_cognition$Flanker)
boxplot(data_cognition$Oral_Reading_Recognition)
boxplot(data_cognition$Rey_Auditory_Verbal_Learning)
boxplot(data_cognition$Little_Man)  
boxplot(data_cognition$Little_Man)$out     # The Little Man score is a percentage, thus value of 0 or above 100 are outliers

data_cognition_clean <- data_cognition %>%       #Before when I had T2 as well there was value over 100 so I was removing them but now the only values that have an issue would be the 0 so not sure if I should remove them or not 
  mutate(Little_Man=ifelse(Little_Man == 0, NA, Little_Man))
boxplot(data_cognition_clean$Little_Man)

#Scaling
data_cognition_clean[,5:9]<-scale(data_cognition_clean[,5:9])
data_greymatter[,5:208]<-scale(data_greymatter[,5:208])
data_whitematter[,5:115]<-scale(data_whitematter[,5:115])

-------------------------------------------------------------------------------------------------------------------------------
#Average every regions of interest bilaterally
-------------------------------------------------------------------------------------------------------------------------------

data_greymatter_lat <- data_greymatter %>%
  mutate(bankssts_ct=(as.numeric(lh_bankssts_ct)+as.numeric(rh_bankssts_ct))/2,
         caudalanteriorcingulate_ct=(as.numeric(lh_caudalanteriorcingulate_ct)+as.numeric(rh_caudalanteriorcingulate_ct))/2,
         caudalmiddlefrontal_ct=(as.numeric(lh_caudalmiddlefrontal_ct)+as.numeric(rh_caudalmiddlefrontal_ct))/2,
         cuneus_ct=(as.numeric(lh_cuneus_ct)+as.numeric(rh_cuneus_ct))/2,
         entorhinal_ct=(as.numeric(lh_entorhinal_ct)+as.numeric(rh_entorhinal_ct))/2,
         fusiform_ct=(as.numeric(lh_fusiform_ct)+as.numeric(rh_fusiform_ct))/2,
         inferiorparietal_ct=(as.numeric(lh_inferiorparietal_ct)+as.numeric(rh_inferiorparietal_ct))/2,
         inferiortemporal_ct=(as.numeric(lh_inferiortemporal_ct)+as.numeric(rh_inferiortemporal_ct))/2,
         isthmuscingulate_ct=(as.numeric(lh_isthmuscingulate_ct)+as.numeric(rh_isthmuscingulate_ct))/2,
         lateraloccipital_ct=(as.numeric(lh_lateraloccipital_ct)+as.numeric(rh_lateraloccipital_ct))/2,
         lateralorbitofrontal_ct=(as.numeric(lh_lateralorbitofrontal_ct)+as.numeric(rh_lateralorbitofrontal_ct))/2,
         lingual_ct=(as.numeric(lh_lingual_ct)+as.numeric(rh_lingual_ct))/2,
         medialorbitofrontal_ct=(as.numeric(lh_medialorbitofrontal_ct)+as.numeric(rh_medialorbitofrontal_ct))/2,
         middletemporal_ct=(as.numeric(lh_middletemporal_ct)+as.numeric(rh_middletemporal_ct))/2,
         parahippocampal_ct=(as.numeric(lh_parahippocampal_ct)+as.numeric(rh_parahippocampal_ct))/2,
         paracentral_ct=(as.numeric(lh_paracentral_ct)+as.numeric(rh_paracentral_ct))/2,
         parsopercularis_ct=(as.numeric(lh_parsopercularis_ct)+as.numeric(rh_parsopercularis_ct))/2,
         parsorbitalis_ct=(as.numeric(lh_parsorbitalis_ct)+as.numeric(rh_parsorbitalis_ct))/2,
         parstriangularis_ct=(as.numeric(lh_parstriangularis_ct)+as.numeric(rh_parstriangularis_ct))/2,
         pericalcarine_ct=(as.numeric(lh_pericalcarine_ct)+as.numeric(rh_pericalcarine_ct))/2,
         postcentral_ct=(as.numeric(lh_postcentral_ct)+as.numeric(rh_postcentral_ct))/2,
         posteriorcingulate_ct=(as.numeric(lh_posteriorcingulate_ct)+as.numeric(rh_posteriorcingulate_ct))/2,
         precentral_ct=(as.numeric(lh_precentral_ct)+as.numeric(rh_precentral_ct))/2,
         precuneus_ct=(as.numeric(lh_precuneus_ct)+as.numeric(rh_precuneus_ct))/2,
         rostralanteriorcingulate_ct=(as.numeric(lh_rostralanteriorcingulate_ct)+as.numeric(rh_rostralanteriorcingulate_ct))/2,
         rostralmiddlefrontal_ct=(as.numeric(lh_rostralmiddlefrontal_ct)+as.numeric(rh_rostralmiddlefrontal_ct))/2,
         superiorfrontal_ct=(as.numeric(lh_superiorfrontal_ct)+as.numeric(rh_superiorfrontal_ct))/2,
         superiorparietal_ct=(as.numeric(lh_superiorparietal_ct)+as.numeric(rh_superiorparietal_ct))/2,
         superiortemporal_ct=(as.numeric(lh_superiortemporal_ct)+as.numeric(rh_superiortemporal_ct))/2,
         supramarginal_ct=(as.numeric(lh_supramarginal_ct)+as.numeric(rh_supramarginal_ct))/2,
         frontalpole_ct=(as.numeric(lh_frontalpole_ct)+as.numeric(rh_frontalpole_ct))/2,
         temporalpole_ct=(as.numeric(lh_temporalpole_ct)+as.numeric(rh_temporalpole_ct))/2,
         transversetemporal_ct=(as.numeric(lh_transversetemporal_ct)+as.numeric(rh_transversetemporal_ct))/2,
         insula_ct=(as.numeric(lh_insula_ct)+as.numeric(rh_insula_ct))/2) %>%
  select(-one_of("lh_bankssts_ct","lh_caudalanteriorcingulate_ct","lh_caudalmiddlefrontal_ct","lh_cuneus_ct","lh_entorhinal_ct","lh_fusiform_ct","lh_inferiorparietal_ct","lh_inferiortemporal_ct","lh_isthmuscingulate_ct","lh_lateraloccipital_ct","lh_lateralorbitofrontal_ct","lh_lingual_ct","lh_medialorbitofrontal_ct","lh_middletemporal_ct","lh_parahippocampal_ct","lh_paracentral_ct","lh_parsopercularis_ct","lh_parsorbitalis_ct","lh_parstriangularis_ct","lh_pericalcarine_ct","lh_postcentral_ct","lh_posteriorcingulate_ct","lh_precentral_ct","lh_precuneus_ct","lh_rostralanteriorcingulate_ct","lh_rostralmiddlefrontal_ct","lh_superiorfrontal_ct","lh_superiorparietal_ct","lh_superiortemporal_ct","lh_supramarginal_ct","lh_frontalpole_ct","lh_temporalpole_ct","lh_transversetemporal_ct","lh_insula_ct"))%>%
  select(-one_of("rh_bankssts_ct","rh_caudalanteriorcingulate_ct","rh_caudalmiddlefrontal_ct","rh_cuneus_ct","rh_entorhinal_ct","rh_fusiform_ct","rh_inferiorparietal_ct","rh_inferiortemporal_ct","rh_isthmuscingulate_ct","rh_lateraloccipital_ct","rh_lateralorbitofrontal_ct","rh_lingual_ct","rh_medialorbitofrontal_ct","rh_middletemporal_ct","rh_parahippocampal_ct","rh_paracentral_ct","rh_parsopercularis_ct","rh_parsorbitalis_ct","rh_parstriangularis_ct","rh_pericalcarine_ct","rh_postcentral_ct","rh_posteriorcingulate_ct","rh_precentral_ct","rh_precuneus_ct","rh_rostralanteriorcingulate_ct","rh_rostralmiddlefrontal_ct","rh_superiorfrontal_ct","rh_superiorparietal_ct","rh_superiortemporal_ct","rh_supramarginal_ct","rh_frontalpole_ct","rh_temporalpole_ct","rh_transversetemporal_ct","rh_insula_ct"))%>%
  mutate(bankssts_sa=(as.numeric(lh_bankssts_sa)+as.numeric(rh_bankssts_sa))/2,
         caudalanteriorcingulate_sa=(as.numeric(lh_caudalanteriorcingulate_sa)+as.numeric(rh_caudalanteriorcingulate_sa))/2,
         caudalmiddlefrontal_sa=(as.numeric(lh_caudalmiddlefrontal_sa)+as.numeric(rh_caudalmiddlefrontal_sa))/2,
         cuneus_sa=(as.numeric(lh_cuneus_sa)+as.numeric(rh_cuneus_sa))/2,
         entorhinal_sa=(as.numeric(lh_entorhinal_sa)+as.numeric(rh_entorhinal_sa))/2,
         fusiform_sa=(as.numeric(lh_fusiform_sa)+as.numeric(rh_fusiform_sa))/2,
         inferiorparietal_sa=(as.numeric(lh_inferiorparietal_sa)+as.numeric(rh_inferiorparietal_sa))/2,
         inferiortemporal_sa=(as.numeric(lh_inferiortemporal_sa)+as.numeric(rh_inferiortemporal_sa))/2,
         isthmuscingulate_sa=(as.numeric(lh_isthmuscingulate_sa)+as.numeric(rh_isthmuscingulate_sa))/2,
         lateraloccipital_sa=(as.numeric(lh_lateraloccipital_sa)+as.numeric(rh_lateraloccipital_sa))/2,
         lateralorbitofrontal_sa=(as.numeric(lh_lateralorbitofrontal_sa)+as.numeric(rh_lateralorbitofrontal_sa))/2,
         lingual_sa=(as.numeric(lh_lingual_sa)+as.numeric(rh_lingual_sa))/2,
         medialorbitofrontal_sa=(as.numeric(lh_medialorbitofrontal_sa)+as.numeric(rh_medialorbitofrontal_sa))/2,
         middletemporal_sa=(as.numeric(lh_middletemporal_sa)+as.numeric(rh_middletemporal_sa))/2,
         parahippocampal_sa=(as.numeric(lh_parahippocampal_sa)+as.numeric(rh_parahippocampal_sa))/2,
         paracentral_sa=(as.numeric(lh_paracentral_sa)+as.numeric(rh_paracentral_sa))/2,
         parsopercularis_sa=(as.numeric(lh_parsopercularis_sa)+as.numeric(rh_parsopercularis_sa))/2,
         parsorbitalis_sa=(as.numeric(lh_parsorbitalis_sa)+as.numeric(rh_parsorbitalis_sa))/2,
         parstriangularis_sa=(as.numeric(lh_parstriangularis_sa)+as.numeric(rh_parstriangularis_sa))/2,
         pericalcarine_sa=(as.numeric(lh_pericalcarine_sa)+as.numeric(rh_pericalcarine_sa))/2,
         postcentral_sa=(as.numeric(lh_postcentral_sa)+as.numeric(rh_postcentral_sa))/2,
         posteriorcingulate_sa=(as.numeric(lh_posteriorcingulate_sa)+as.numeric(rh_posteriorcingulate_sa))/2,
         precentral_sa=(as.numeric(lh_precentral_sa)+as.numeric(rh_precentral_sa))/2,
         precuneus_sa=(as.numeric(lh_precuneus_sa)+as.numeric(rh_precuneus_sa))/2,
         rostralanteriorcingulate_sa=(as.numeric(lh_rostralanteriorcingulate_sa)+as.numeric(rh_rostralanteriorcingulate_sa))/2,
         rostralmiddlefrontal_sa=(as.numeric(lh_rostralmiddlefrontal_sa)+as.numeric(rh_rostralmiddlefrontal_sa))/2,
         superiorfrontal_sa=(as.numeric(lh_superiorfrontal_sa)+as.numeric(rh_superiorfrontal_sa))/2,
         superiorparietal_sa=(as.numeric(lh_superiorparietal_sa)+as.numeric(rh_superiorparietal_sa))/2,
         superiortemporal_sa=(as.numeric(lh_superiortemporal_sa)+as.numeric(rh_superiortemporal_sa))/2,
         supramarginal_sa=(as.numeric(lh_supramarginal_sa)+as.numeric(rh_supramarginal_sa))/2,
         frontalpole_sa=(as.numeric(lh_frontalpole_sa)+as.numeric(rh_frontalpole_sa))/2,
         temporalpole_sa=(as.numeric(lh_temporalpole_sa)+as.numeric(rh_temporalpole_sa))/2,
         transversetemporal_sa=(as.numeric(lh_transversetemporal_sa)+as.numeric(rh_transversetemporal_sa))/2,
         insula_sa=(as.numeric(lh_insula_sa)+as.numeric(rh_insula_sa))/2) %>%
  select(-one_of("lh_bankssts_sa","lh_caudalanteriorcingulate_sa","lh_caudalmiddlefrontal_sa","lh_cuneus_sa","lh_entorhinal_sa","lh_fusiform_sa","lh_inferiorparietal_sa","lh_inferiortemporal_sa","lh_isthmuscingulate_sa","lh_lateraloccipital_sa","lh_lateralorbitofrontal_sa","lh_lingual_sa","lh_medialorbitofrontal_sa","lh_middletemporal_sa","lh_parahippocampal_sa","lh_paracentral_sa","lh_parsopercularis_sa","lh_parsorbitalis_sa","lh_parstriangularis_sa","lh_pericalcarine_sa","lh_postcentral_sa","lh_posteriorcingulate_sa","lh_precentral_sa","lh_precuneus_sa","lh_rostralanteriorcingulate_sa","lh_rostralmiddlefrontal_sa","lh_superiorfrontal_sa","lh_superiorparietal_sa","lh_superiortemporal_sa","lh_supramarginal_sa","lh_frontalpole_sa","lh_temporalpole_sa","lh_transversetemporal_sa","lh_insula_sa"))%>%
  select(-one_of("rh_bankssts_sa","rh_caudalanteriorcingulate_sa","rh_caudalmiddlefrontal_sa","rh_cuneus_sa","rh_entorhinal_sa","rh_fusiform_sa","rh_inferiorparietal_sa","rh_inferiortemporal_sa","rh_isthmuscingulate_sa","rh_lateraloccipital_sa","rh_lateralorbitofrontal_sa","rh_lingual_sa","rh_medialorbitofrontal_sa","rh_middletemporal_sa","rh_parahippocampal_sa","rh_paracentral_sa","rh_parsopercularis_sa","rh_parsorbitalis_sa","rh_parstriangularis_sa","rh_pericalcarine_sa","rh_postcentral_sa","rh_posteriorcingulate_sa","rh_precentral_sa","rh_precuneus_sa","rh_rostralanteriorcingulate_sa","rh_rostralmiddlefrontal_sa","rh_superiorfrontal_sa","rh_superiorparietal_sa","rh_superiortemporal_sa","rh_supramarginal_sa","rh_frontalpole_sa","rh_temporalpole_sa","rh_transversetemporal_sa","rh_insula_sa"))%>%
  mutate(bankssts_gmv=(as.numeric(lh_bankssts_gmv)+as.numeric(rh_bankssts_gmv))/2,
         caudalanteriorcingulate_gmv=(as.numeric(lh_caudalanteriorcingulate_gmv)+as.numeric(rh_caudalanteriorcingulate_gmv))/2,
         caudalmiddlefrontal_gmv=(as.numeric(lh_caudalmiddlefrontal_gmv)+as.numeric(rh_caudalmiddlefrontal_gmv))/2,
         cuneus_gmv=(as.numeric(lh_cuneus_gmv)+as.numeric(rh_cuneus_gmv))/2,
         entorhinal_gmv=(as.numeric(lh_entorhinal_gmv)+as.numeric(rh_entorhinal_gmv))/2,
         fusiform_gmv=(as.numeric(lh_fusiform_gmv)+as.numeric(rh_fusiform_gmv))/2,
         inferiorparietal_gmv=(as.numeric(lh_inferiorparietal_gmv)+as.numeric(rh_inferiorparietal_gmv))/2,
         inferiortemporal_gmv=(as.numeric(lh_inferiortemporal_gmv)+as.numeric(rh_inferiortemporal_gmv))/2,
         isthmuscingulate_gmv=(as.numeric(lh_isthmuscingulate_gmv)+as.numeric(rh_isthmuscingulate_gmv))/2,
         lateraloccipital_gmv=(as.numeric(lh_lateraloccipital_gmv)+as.numeric(rh_lateraloccipital_gmv))/2,
         lateralorbitofrontal_gmv=(as.numeric(lh_lateralorbitofrontal_gmv)+as.numeric(rh_lateralorbitofrontal_gmv))/2,
         lingual_gmv=(as.numeric(lh_lingual_gmv)+as.numeric(rh_lingual_gmv))/2,
         medialorbitofrontal_gmv=(as.numeric(lh_medialorbitofrontal_gmv)+as.numeric(rh_medialorbitofrontal_gmv))/2,
         middletemporal_gmv=(as.numeric(lh_middletemporal_gmv)+as.numeric(rh_middletemporal_gmv))/2,
         parahippocampal_gmv=(as.numeric(lh_parahippocampal_gmv)+as.numeric(rh_parahippocampal_gmv))/2,
         paracentral_gmv=(as.numeric(lh_paracentral_gmv)+as.numeric(rh_paracentral_gmv))/2,
         parsopercularis_gmv=(as.numeric(lh_parsopercularis_gmv)+as.numeric(rh_parsopercularis_gmv))/2,
         parsorbitalis_gmv=(as.numeric(lh_parsorbitalis_gmv)+as.numeric(rh_parsorbitalis_gmv))/2,
         parstriangularis_gmv=(as.numeric(lh_parstriangularis_gmv)+as.numeric(rh_parstriangularis_gmv))/2,
         pericalcarine_gmv=(as.numeric(lh_pericalcarine_gmv)+as.numeric(rh_pericalcarine_gmv))/2,
         postcentral_gmv=(as.numeric(lh_postcentral_gmv)+as.numeric(rh_postcentral_gmv))/2,
         posteriorcingulate_gmv=(as.numeric(lh_posteriorcingulate_gmv)+as.numeric(rh_posteriorcingulate_gmv))/2,
         precentral_gmv=(as.numeric(lh_precentral_gmv)+as.numeric(rh_precentral_gmv))/2,
         precuneus_gmv=(as.numeric(lh_precuneus_gmv)+as.numeric(rh_precuneus_gmv))/2,
         rostralanteriorcingulate_gmv=(as.numeric(lh_rostralanteriorcingulate_gmv)+as.numeric(rh_rostralanteriorcingulate_gmv))/2,
         rostralmiddlefrontal_gmv=(as.numeric(lh_rostralmiddlefrontal_gmv)+as.numeric(rh_rostralmiddlefrontal_gmv))/2,
         superiorfrontal_gmv=(as.numeric(lh_superiorfrontal_gmv)+as.numeric(rh_superiorfrontal_gmv))/2,
         superiorparietal_gmv=(as.numeric(lh_superiorparietal_gmv)+as.numeric(rh_superiorparietal_gmv))/2,
         superiortemporal_gmv=(as.numeric(lh_superiortemporal_gmv)+as.numeric(rh_superiortemporal_gmv))/2,
         supramarginal_gmv=(as.numeric(lh_supramarginal_gmv)+as.numeric(rh_supramarginal_gmv))/2,
         frontalpole_gmv=(as.numeric(lh_frontalpole_gmv)+as.numeric(rh_frontalpole_gmv))/2,
         temporalpole_gmv=(as.numeric(lh_temporalpole_gmv)+as.numeric(rh_temporalpole_gmv))/2,
         transversetemporal_gmv=(as.numeric(lh_transversetemporal_gmv)+as.numeric(rh_transversetemporal_gmv))/2,
         insula_gmv=(as.numeric(lh_insula_gmv)+as.numeric(rh_insula_gmv))/2) %>%
  select(-one_of("lh_bankssts_gmv","lh_caudalanteriorcingulate_gmv","lh_caudalmiddlefrontal_gmv","lh_cuneus_gmv","lh_entorhinal_gmv","lh_fusiform_gmv","lh_inferiorparietal_gmv","lh_inferiortemporal_gmv","lh_isthmuscingulate_gmv","lh_lateraloccipital_gmv","lh_lateralorbitofrontal_gmv","lh_lingual_gmv","lh_medialorbitofrontal_gmv","lh_middletemporal_gmv","lh_parahippocampal_gmv","lh_paracentral_gmv","lh_parsopercularis_gmv","lh_parsorbitalis_gmv","lh_parstriangularis_gmv","lh_pericalcarine_gmv","lh_postcentral_gmv","lh_posteriorcingulate_gmv","lh_precentral_gmv","lh_precuneus_gmv","lh_rostralanteriorcingulate_gmv","lh_rostralmiddlefrontal_gmv","lh_superiorfrontal_gmv","lh_superiorparietal_gmv","lh_superiortemporal_gmv","lh_supramarginal_gmv","lh_frontalpole_gmv","lh_temporalpole_gmv","lh_transversetemporal_gmv","lh_insula_gmv"))%>%
  select(-one_of("rh_bankssts_gmv","rh_caudalanteriorcingulate_gmv","rh_caudalmiddlefrontal_gmv","rh_cuneus_gmv","rh_entorhinal_gmv","rh_fusiform_gmv","rh_inferiorparietal_gmv","rh_inferiortemporal_gmv","rh_isthmuscingulate_gmv","rh_lateraloccipital_gmv","rh_lateralorbitofrontal_gmv","rh_lingual_gmv","rh_medialorbitofrontal_gmv","rh_middletemporal_gmv","rh_parahippocampal_gmv","rh_paracentral_gmv","rh_parsopercularis_gmv","rh_parsorbitalis_gmv","rh_parstriangularis_gmv","rh_pericalcarine_gmv","rh_postcentral_gmv","rh_posteriorcingulate_gmv","rh_precentral_gmv","rh_precuneus_gmv","rh_rostralanteriorcingulate_gmv","rh_rostralmiddlefrontal_gmv","rh_superiorfrontal_gmv","rh_superiorparietal_gmv","rh_superiortemporal_gmv","rh_supramarginal_gmv","rh_frontalpole_gmv","rh_temporalpole_gmv","rh_transversetemporal_gmv","rh_insula_gmv"))

data_whitematter_lat <- data_whitematter %>%
  mutate(fornix_fa=(as.numeric(lh_fornix_fa)+as.numeric(rh_fornix_fa))/2,
         cingulatecingulum_fa=(as.numeric(lh_cingulatecingulum_fa)+as.numeric(rh_cingulatecingulum_fa))/2,
         parahippocampalcingulum_fa=(as.numeric(lh_parahippocampalcingulum_fa)+as.numeric(rh_parahippocampalcingulum_fa))/2,
         corticospinalpyramidal_fa=(as.numeric(lh_corticospinalpyramidal_fa)+as.numeric(rh_corticospinalpyramidal_fa))/2,
         anteriorthalamicradiations_fa=(as.numeric(lh_anteriorthalamicradiations_fa)+as.numeric(rh_anteriorthalamicradiations_fa))/2,
         uncinate_fa=(as.numeric(lh_uncinate_fa)+as.numeric(rh_uncinate_fa))/2,
         inferiorlongitudinalfasiculus_fa=(as.numeric(lh_inferiorlongitudinalfasiculus_fa)+as.numeric(rh_inferiorlongitudinalfasiculus_fa))/2,
         inferiorfrontooccipitalfasiculus_fa=(as.numeric(lh_inferiorfrontooccipitalfasiculus_fa)+as.numeric(rh_inferiorfrontooccipitalfasiculus_fa))/2,
         superiorlongitudinalfasiculus_fa=(as.numeric(lh_superiorlongitudinalfasiculus_fa)+as.numeric(rh_superiorlongitudinalfasiculus_fa))/2,
         temporalsuperiorlongitudinalfasiculus_fa=(as.numeric(lh_temporalsuperiorlongitudinalfasiculus_fa)+as.numeric(rh_temporalsuperiorlongitudinalfasiculus_fa))/2,
         parietalsuperiorlongitudinalfasiculus_fa=(as.numeric(lh_parietalsuperiorlongitudinalfasiculus_fa)+as.numeric(rh_parietalsuperiorlongitudinalfasiculus_fa))/2,
         superiorcorticostriate_fa=(as.numeric(lh_superiorcorticostriate_fa)+as.numeric(rh_superiorcorticostriate_fa))/2,
         superiorcorticostriatefrontalcortex_fa=(as.numeric(lh_superiorcorticostriatefrontalcortex_fa)+as.numeric(rh_superiorcorticostriatefrontalcortex_fa))/2,
         superiorcorticostriateparietalcortex_fa=(as.numeric(lh_superiorcorticostriateparietalcortex_fa)+as.numeric(rh_superiorcorticostriateparietalcortex_fa))/2,
         striatalinferiorfrontalcortex_fa=(as.numeric(lh_striatalinferiorfrontalcortex_fa)+as.numeric(rh_striatalinferiorfrontalcortex_fa))/2,
         inferiorfrontalsuperiorfrontalcortex_fa=(as.numeric(lh_inferiorfrontalsuperiorfrontalcortex_fa)+as.numeric(rh_inferiorfrontalsuperiorfrontalcortex_fa))/2,
         fornix_exfimbria_fa=(as.numeric(lh_fornix_exfimbria_fa)+as.numeric(rh_fornix_exfimbria_fa))/2) %>%
  select(-one_of("lh_fornix_fa","lh_cingulatecingulum_fa","lh_parahippocampalcingulum_fa","lh_corticospinalpyramidal_fa","lh_anteriorthalamicradiations_fa","lh_uncinate_fa","lh_inferiorlongitudinalfasiculus_fa","lh_inferiorfrontooccipitalfasiculus_fa","lh_superiorlongitudinalfasiculus_fa","lh_temporalsuperiorlongitudinalfasiculus_fa","lh_parietalsuperiorlongitudinalfasiculus_fa","lh_superiorcorticostriate_fa","lh_superiorcorticostriatefrontalcortex_fa","lh_superiorcorticostriateparietalcortex_fa","lh_striatalinferiorfrontalcortex_fa","lh_inferiorfrontalsuperiorfrontalcortex_fa","lh_fornix_exfimbria_fa"))%>%
  select(-one_of("rh_fornix_fa","rh_cingulatecingulum_fa","rh_parahippocampalcingulum_fa","rh_corticospinalpyramidal_fa","rh_anteriorthalamicradiations_fa","rh_uncinate_fa","rh_inferiorlongitudinalfasiculus_fa","rh_inferiorfrontooccipitalfasiculus_fa","rh_superiorlongitudinalfasiculus_fa","rh_temporalsuperiorlongitudinalfasiculus_fa","rh_parietalsuperiorlongitudinalfasiculus_fa","rh_superiorcorticostriate_fa","rh_superiorcorticostriatefrontalcortex_fa","rh_superiorcorticostriateparietalcortex_fa","rh_striatalinferiorfrontalcortex_fa","rh_inferiorfrontalsuperiorfrontalcortex_fa","rh_fornix_exfimbria_fa"))%>%
  mutate(fornix_md=(as.numeric(lh_fornix_md)+as.numeric(rh_fornix_md))/2,
         cingulatecingulum_md=(as.numeric(lh_cingulatecingulum_md)+as.numeric(rh_cingulatecingulum_md))/2,
         parahippocampalcingulum_md=(as.numeric(lh_parahippocampalcingulum_md)+as.numeric(rh_parahippocampalcingulum_md))/2,
         corticospinalpyramidal_md=(as.numeric(lh_corticospinalpyramidal_md)+as.numeric(rh_corticospinalpyramidal_md))/2,
         anteriorthalamicradiations_md=(as.numeric(lh_anteriorthalamicradiations_md)+as.numeric(rh_anteriorthalamicradiations_md))/2,
         uncinate_md=(as.numeric(lh_uncinate_md)+as.numeric(rh_uncinate_md))/2,
         inferiorlongitudinalfasiculus_md=(as.numeric(lh_inferiorlongitudinalfasiculus_md)+as.numeric(rh_inferiorlongitudinalfasiculus_md))/2,
         inferiorfrontooccipitalfasiculus_md=(as.numeric(lh_inferiorfrontooccipitalfasiculus_md)+as.numeric(rh_inferiorfrontooccipitalfasiculus_md))/2,
         superiorlongitudinalfasiculus_md=(as.numeric(lh_superiorlongitudinalfasiculus_md)+as.numeric(rh_superiorlongitudinalfasiculus_md))/2,
         temporalsuperiorlongitudinalfasiculus_md=(as.numeric(lh_temporalsuperiorlongitudinalfasiculus_md)+as.numeric(rh_temporalsuperiorlongitudinalfasiculus_md))/2,
         parietalsuperiorlongitudinalfasiculus_md=(as.numeric(lh_parietalsuperiorlongitudinalfasiculus_md)+as.numeric(rh_parietalsuperiorlongitudinalfasiculus_md))/2,
         superiorcorticostriate_md=(as.numeric(lh_superiorcorticostriate_md)+as.numeric(rh_superiorcorticostriate_md))/2,
         superiorcorticostriatefrontalcortex_md=(as.numeric(lh_superiorcorticostriatefrontalcortex_md)+as.numeric(rh_superiorcorticostriatefrontalcortex_md))/2,
         superiorcorticostriateparietalcortex_md=(as.numeric(lh_superiorcorticostriateparietalcortex_md)+as.numeric(rh_superiorcorticostriateparietalcortex_md))/2,
         striatalinferiorfrontalcortex_md=(as.numeric(lh_striatalinferiorfrontalcortex_md)+as.numeric(rh_striatalinferiorfrontalcortex_md))/2,
         inferiorfrontalsuperiorfrontalcortex_md=(as.numeric(lh_inferiorfrontalsuperiorfrontalcortex_md)+as.numeric(rh_inferiorfrontalsuperiorfrontalcortex_md))/2,
         fornix_exfimbria_md=(as.numeric(lh_fornix_exfimbria_md)+as.numeric(rh_fornix_exfimbria_md))/2) %>%
  select(-one_of("lh_fornix_md","lh_cingulatecingulum_md","lh_parahippocampalcingulum_md","lh_corticospinalpyramidal_md","lh_anteriorthalamicradiations_md","lh_uncinate_md","lh_inferiorlongitudinalfasiculus_md","lh_inferiorfrontooccipitalfasiculus_md","lh_superiorlongitudinalfasiculus_md","lh_temporalsuperiorlongitudinalfasiculus_md","lh_parietalsuperiorlongitudinalfasiculus_md","lh_superiorcorticostriate_md","lh_superiorcorticostriatefrontalcortex_md","lh_superiorcorticostriateparietalcortex_md","lh_striatalinferiorfrontalcortex_md","lh_inferiorfrontalsuperiorfrontalcortex_md","lh_fornix_exfimbria_md"))%>%
  select(-one_of("rh_fornix_md","rh_cingulatecingulum_md","rh_parahippocampalcingulum_md","rh_corticospinalpyramidal_md","rh_anteriorthalamicradiations_md","rh_uncinate_md","rh_inferiorlongitudinalfasiculus_md","rh_inferiorfrontooccipitalfasiculus_md","rh_superiorlongitudinalfasiculus_md","rh_temporalsuperiorlongitudinalfasiculus_md","rh_parietalsuperiorlongitudinalfasiculus_md","rh_superiorcorticostriate_md","rh_superiorcorticostriatefrontalcortex_md","rh_superiorcorticostriateparietalcortex_md","rh_striatalinferiorfrontalcortex_md","rh_inferiorfrontalsuperiorfrontalcortex_md","rh_fornix_exfimbria_md"))%>%
  mutate(fornix_wmv=(as.numeric(lh_fornix_wmv)+as.numeric(rh_fornix_wmv))/2,
         cingulatecingulum_wmv=(as.numeric(lh_cingulatecingulum_wmv)+as.numeric(rh_cingulatecingulum_wmv))/2,
         parahippocampalcingulum_wmv=(as.numeric(lh_parahippocampalcingulum_wmv)+as.numeric(rh_parahippocampalcingulum_wmv))/2,
         corticospinalpyramidal_wmv=(as.numeric(lh_corticospinalpyramidal_wmv)+as.numeric(rh_corticospinalpyramidal_wmv))/2,
         anteriorthalamicradiations_wmv=(as.numeric(lh_anteriorthalamicradiations_wmv)+as.numeric(rh_anteriorthalamicradiations_wmv))/2,
         uncinate_wmv=(as.numeric(lh_uncinate_wmv)+as.numeric(rh_uncinate_wmv))/2,
         inferiorlongitudinalfasiculus_wmv=(as.numeric(lh_inferiorlongitudinalfasiculus_wmv)+as.numeric(rh_inferiorlongitudinalfasiculus_wmv))/2,
         inferiorfrontooccipitalfasiculus_wmv=(as.numeric(lh_inferiorfrontooccipitalfasiculus_wmv)+as.numeric(rh_inferiorfrontooccipitalfasiculus_wmv))/2,
         superiorlongitudinalfasiculus_wmv=(as.numeric(lh_superiorlongitudinalfasiculus_wmv)+as.numeric(rh_superiorlongitudinalfasiculus_wmv))/2,
         temporalsuperiorlongitudinalfasiculus_wmv=(as.numeric(lh_temporalsuperiorlongitudinalfasiculus_wmv)+as.numeric(rh_temporalsuperiorlongitudinalfasiculus_wmv))/2,
         parietalsuperiorlongitudinalfasiculus_wmv=(as.numeric(lh_parietalsuperiorlongitudinalfasiculus_wmv)+as.numeric(rh_parietalsuperiorlongitudinalfasiculus_wmv))/2,
         superiorcorticostriate_wmv=(as.numeric(lh_superiorcorticostriate_wmv)+as.numeric(rh_superiorcorticostriate_wmv))/2,
         superiorcorticostriatefrontalcortex_wmv=(as.numeric(lh_superiorcorticostriatefrontalcortex_wmv)+as.numeric(rh_superiorcorticostriatefrontalcortex_wmv))/2,
         superiorcorticostriateparietalcortex_wmv=(as.numeric(lh_superiorcorticostriateparietalcortex_wmv)+as.numeric(rh_superiorcorticostriateparietalcortex_wmv))/2,
         striatalinferiorfrontalcortex_wmv=(as.numeric(lh_striatalinferiorfrontalcortex_wmv)+as.numeric(rh_striatalinferiorfrontalcortex_wmv))/2,
         inferiorfrontalsuperiorfrontalcortex_wmv=(as.numeric(lh_inferiorfrontalsuperiorfrontalcortex_wmv)+as.numeric(rh_inferiorfrontalsuperiorfrontalcortex_wmv))/2,
         fornix_exfimbria_wmv=(as.numeric(lh_fornix_exfimbria_wmv)+as.numeric(rh_fornix_exfimbria_wmv))/2) %>%
  select(-one_of("lh_fornix_wmv","lh_cingulatecingulum_wmv","lh_parahippocampalcingulum_wmv","lh_corticospinalpyramidal_wmv","lh_anteriorthalamicradiations_wmv","lh_uncinate_wmv","lh_inferiorlongitudinalfasiculus_wmv","lh_inferiorfrontooccipitalfasiculus_wmv","lh_superiorlongitudinalfasiculus_wmv","lh_temporalsuperiorlongitudinalfasiculus_wmv","lh_parietalsuperiorlongitudinalfasiculus_wmv","lh_superiorcorticostriate_wmv","lh_superiorcorticostriatefrontalcortex_wmv","lh_superiorcorticostriateparietalcortex_wmv","lh_striatalinferiorfrontalcortex_wmv","lh_inferiorfrontalsuperiorfrontalcortex_wmv","lh_fornix_exfimbria_wmv"))%>%
  select(-one_of("rh_fornix_wmv","rh_cingulatecingulum_wmv","rh_parahippocampalcingulum_wmv","rh_corticospinalpyramidal_wmv","rh_anteriorthalamicradiations_wmv","rh_uncinate_wmv","rh_inferiorlongitudinalfasiculus_wmv","rh_inferiorfrontooccipitalfasiculus_wmv","rh_superiorlongitudinalfasiculus_wmv","rh_temporalsuperiorlongitudinalfasiculus_wmv","rh_parietalsuperiorlongitudinalfasiculus_wmv","rh_superiorcorticostriate_wmv","rh_superiorcorticostriatefrontalcortex_wmv","rh_superiorcorticostriateparietalcortex_wmv","rh_striatalinferiorfrontalcortex_wmv","rh_inferiorfrontalsuperiorfrontalcortex_wmv","rh_fornix_exfimbria_wmv"))

-------------------------------------------------------------------------------------------------------------------------------
#Creation of the merge and processed datasets
-------------------------------------------------------------------------------------------------------------------------------

#Merge all the datasets together
data_brain<-merge(data_greymatter_lat,data_whitematter_lat,by=c("subjectkey","eventname","interview_age","sex"),all=T)
data_total<-merge(data_cognition_clean,data_brain,by=c("subjectkey","eventname","interview_age","sex"),all=T)

#Subsample 15% of the data to test the subsequent models
set.seed(2)
subsample<-sample(1:nrow(data_total),round(0.15*nrow(data_total)),replace=F)
data_total_subsample15<-data_total[c(subsample),]
data_total_subsample85<-data_total[-c(subsample),]

setwd("C:/Users/leamic/Documents/Project 1/Data/Processed data")
write.csv(data_total_subsample15,'data_total_subsample15.csv',row.names=F)
write.csv(data_total_subsample85,'data_total_subsample85.csv',row.names=F)

-------------------------------------------------------------------------------------------------------------------------------
#Descriptives
-------------------------------------------------------------------------------------------------------------------------------
  
#Visualize each variable
view(dfSummary(data_total[,c(1:9,10,112)]))  

#Correlations of the cognitive tasks (skip id, event, age and sex)
corrplot(cor(data_total[,c(5:9)],use="pairwise.complete.obs"))
cor(data_total[,c(5:9)],use="pairwise.complete.obs")
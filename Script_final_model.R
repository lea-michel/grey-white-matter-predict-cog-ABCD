#############################################################################################################################
#                                                                                                                           #
#                           Grey and white matter metrics demonstrate distinct and complementary                            #
#                                   prediction of differences in cognitive performance                                      #
#                                     in children: Findings from ABCD (N= 11 876)                                           #
#                                                                                                                           #
#                                            Michel et al, 2023,                                                            #
#                                                                                                                           #
#---------------------------------------------------------------------------------------------------------------------------#
#                                                                                                                           #
#                                      URL & DOI:                                                                           #
#                                 Supplementary materials: https://osf.io/eygwz                                             #
#                                                                                                                           #
#---------------------------------------------------------------------------------------------------------------------------#
#                                                                                                                           #
#                                         Script for computing the models                                                   #
#                                                                                                                           #
#############################################################################################################################



# -------------- Table of contents --------------------------------------------------------------------------------------------
# -------------- 1. Load libraries (line 56) ----------------------------------------------------------------------------------
# -------------- 2. Load the datasets (line 79) -------------------------------------------------------------------------------
# -------------- 3. Confirmatory Factor Analysis for the cognitive factor (line 87) -------------------------------------------

# -------------- 4. Sample 15%: Individual models grey matter (line 149) ------------------------------------------------------
# -------------- 5. Sample 15%: Models per metric - Cortical Thickness (line 1178) --------------------------------------------
# -------------- 6. Sample 15%: Models per metric - Surface Area (line 1292) --------------------------------------------------
# -------------- 7. Sample 15%: Models per metric - Grey Matter Volume (line 1403) --------------------------------------------
# -------------- 8. Sample 15%: Models including the three grey matter metrics (line 1517) ------------------------------------
# -------------- 9. Sample 15%: Individual models white matter (line 1718) ----------------------------------------------------
# -------------- 10. Sample 15%: Models per metric - Fractional Anisotropy (line 2320) ----------------------------------------
# -------------- 11. Sample 15%: Models per metric - Mean Diffusivity (line 2408) ---------------------------------------------
# -------------- 12. Sample 15%: Models per metric - White Matter Volume (line 2497) ------------------------------------------
# -------------- 13. Sample 15%: Models including the three white matter metrics (line 2587) ----------------------------------
# -------------- 14. Sample 15%: Models including the grey & white matter metrics (line 2732) ---------------------------------
# -------------- 15. Sample 15%: Compare the regularized models with grey and white matter (line 2814) ------------------------

# -------------- 16. Sample 85%: Individual models grey matter (line 2915) ----------------------------------------------------
# -------------- 17. Sample 85%: Models per metric - Cortical Thickness (line 4092) -------------------------------------------
# -------------- 18. Sample 85%: Models per metric - Surface Area (line 4203) -------------------------------------------------
# -------------- 19. Sample 85%: Models per metric - Grey Matter Volume (line 4317) -------------------------------------------
# -------------- 20. Sample 85%: Models including the three grey matter metrics (line 4431) -----------------------------------
# -------------- 21. Sample 85%: Individual models white matter (line 4688) ---------------------------------------------------
# -------------- 22. Sample 85%: Models per metric - Fractional Anisotropy (line 5394) ----------------------------------------
# -------------- 23. Sample 85%: Models per metric - Mean Diffusivity (line 5482) ---------------------------------------------
# -------------- 24. Sample 85%: Models per metric - White Matter Volume (line 5571) ------------------------------------------
# -------------- 25. Sample 85%: Models including the three white matter metrics (line 5661) ----------------------------------
# -------------- 26. Sample 85%: Models including the grey & white matter metrics (line 5810) ---------------------------------
# -------------- 27. Sample 85%: Compare the regularized models with grey and white matter (line 5956) ------------------------



#------------------------------------------------------------------------------------------------------------------------------
#Load libraries
#------------------------------------------------------------------------------------------------------------------------------
  
library(lavaan)
library(dplyr)
library(tidymodels)
library(psych)
library(tidyr)
library(plyr)
library(udunits2)
library(ggseg)
library(ggseg3d)
library(ggsegJHU)
library(glmnet)
library(semTools)
library(ggplot2)
library(lavaanPlot)
library(ggplot2)
library(qpcR)
library(viridis)
library(knitr)

#------------------------------------------------------------------------------------------------------------------------------
#Load the datasets
#------------------------------------------------------------------------------------------------------------------------------

setwd("C:/Users/leamic/Documents/Project 1/Data/Processed data")
data_total_subsample15<-read.csv(file="data_total_subsample15.csv", header=TRUE, sep=",", dec=".")
data_total_subsample85<-read.csv(file="data_total_subsample85.csv", header=TRUE, sep=",", dec=".")

#------------------------------------------------------------------------------------------------------------------------------
#Confirmatory Factor Analysis for the cognitive factor
#------------------------------------------------------------------------------------------------------------------------------
  
#Fit the cognitive factor for the 15% and 85% samples
Cognitive_model <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man'

fit_cfa_15 <- cfa(Cognitive_model, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_cfa_15, fit.measures=TRUE,rsquare=T,standardized=T)
fit_cfa_85 <- cfa(Cognitive_model, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_cfa_85, fit.measures=TRUE,rsquare=T,standardized=T)

#Visualize the confirmatory factor analysis path
labels <- list(Cognitive_factor = "Cognitive Factor", Picture_Vocabulary = "Picture Vocabulary", Flanker = "Flanker", Oral_Reading_Recognition = "Oral Reading Recognition", Rey_Auditory_Verbal_Learning = "Rey Auditory Verbal Learning", Little_Man = "Little Man")
path_85 <- lavaanPlot(model = fit_cfa_85, labels = labels, coefs = TRUE, stand = TRUE, sig = 0.05)  #standardized regression paths, showing only paths with p<= .05
save_png(path_85,"path_cfa_85.png") 

#Table to extract the parameters of the cfa
paramest_cfa_85 <- data.frame(lhs=rep(NA, times=22),
                           op=rep(NA, times=22),
                           rhs=rep(NA, times=22),
                           est=rep(NA, times=22),
                           se=rep(NA, times=22),
                           z=rep(NA, times=22),
                           pvalue=rep(NA, times=22),
                           ci.lower=rep(NA, times=22),
                           ci.upper=rep(NA, times=22),
                           std.lv=rep(NA, times=22),
                           std.all=rep(NA, times=22),
                           std.nox=rep(NA, times=22))

paramest_cfa_85[,] <-parameterEstimates(fit_cfa_85,standardized = TRUE, rsquare=TRUE)
write.csv(paramest_cfa_85, "paramest_CFA.csv", row.names=FALSE)

#Check out misfit
residuals(fit_cfa_85)
modificationindices(fit_cfa_85,sort = TRUE,maximum.number = 5)
hist(predict(fit_cfa_85),breaks = 15)

#Fit the cognitive factor with age
Cognitive_model_age <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                        Cognitive_factor ~ interview_age'

fit_cfa_85_age <- cfa(Cognitive_model_age, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_cfa_85_age, fit.measures=TRUE,rsquare=T,standardized=T)

#Check out misfit
residuals(fit_cfa_85_age)
modificationindices(fit_cfa_85_age, sort = TRUE,maximum.number = 5)
hist(predict(fit_cfa_85_age),breaks = 15)

#Add a column "Cognitive factor" for the regularization
data_total_subsample15 <- data_total_subsample15 %>%
  mutate(Cognitive_Factor=predict(fit_cfa_15))

#Exploratory Factor Analysis for the cognitive factor
fa.parallel(data_total_subsample15[,5:9], fm = 'minres', fa = 'fa')
factanal(~Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man, factors=1, data=data_total_subsample15, rotation="promax")
fa.parallel(data_total_subsample85[,5:9], fm = 'minres', fa = 'fa')
factanal(~Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man, factors=1, data=data_total_subsample85, rotation="promax")


#------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Individual models grey matter
#------------------------------------------------------------------------------------------------------------------------------

#Create a table to save the regression coefficients or the fit of the ROI model for each ROIs
data_subsample_GMlat_15_GGSEG <- data.frame(label=c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus","entorhinal","fusiform","inferiorparietal","inferiortemporal","isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual","medialorbitofrontal","middletemporal","parahippocampal","paracentral","parsopercularis","parsorbitalis","parstriangularis","pericalcarine","postcentral","posteriorcingulate","precentral","precuneus","rostralanteriorcingulate","rostralmiddlefrontal","superiorfrontal","superiorparietal","superiortemporal","supramarginal","frontalpole","temporalpole","transversetemporal","insula"),
                                              CT_fit=rep(NA, times=34),
                                              CT_p=rep(NA, times=34),
                                              CT_std=rep(NA, times=34),
                                              SA_fit=rep(NA, times=34),
                                              SA_p=rep(NA, times=34),
                                              SA_std=rep(NA, times=34),
                                              GMV_fit=rep(NA, times=34),
                                              GMV_p=rep(NA, times=34),
                                              GMV_std=rep(NA, times=34))

#Cortical Thickness (34 regions)
##bankssts
Model_CT_bankssts <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                      Cognitive_factor ~ bankssts_ct'
fit_sem_CT_bankssts <- sem(Model_CT_bankssts, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_bankssts, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[1,2]<-BIC(fit_sem_CT_bankssts)
data_subsample_GMlat_15_GGSEG[1,3]<-parameterEstimates(fit_sem_CT_bankssts)[6,7]
data_subsample_GMlat_15_GGSEG[1,4]<-lavInspect(fit_sem_CT_bankssts,what = "std.all")$beta[1,2]

##caudalanteriorcingulate
Model_CT_caudalanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ caudalanteriorcingulate_ct'
fit_sem_CT_caudalanteriorcingulate <- sem(Model_CT_caudalanteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_caudalanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[2,2]<-BIC(fit_sem_CT_caudalanteriorcingulate)
data_subsample_GMlat_15_GGSEG[2,3]<-parameterEstimates(fit_sem_CT_caudalanteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[2,4]<-lavInspect(fit_sem_CT_caudalanteriorcingulate,what = "std.all")$beta[1,2]

##caudalmiddlefrontal
Model_CT_caudalmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                 Cognitive_factor ~ caudalmiddlefrontal_ct'
fit_sem_CT_caudalmiddlefrontal <- sem(Model_CT_caudalmiddlefrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_caudalmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[3,2]<-BIC(fit_sem_CT_caudalmiddlefrontal)
data_subsample_GMlat_15_GGSEG[3,3]<-parameterEstimates(fit_sem_CT_caudalmiddlefrontal)[6,7]
data_subsample_GMlat_15_GGSEG[3,4]<-lavInspect(fit_sem_CT_caudalmiddlefrontal,what = "std.all")$beta[1,2]

##cuneus
Model_CT_cuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                    Cognitive_factor ~ cuneus_ct'
fit_sem_CT_cuneus <- sem(Model_CT_cuneus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_cuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[4,2]<-BIC(fit_sem_CT_cuneus)
data_subsample_GMlat_15_GGSEG[4,3]<-parameterEstimates(fit_sem_CT_cuneus)[6,7]
data_subsample_GMlat_15_GGSEG[4,4]<-lavInspect(fit_sem_CT_cuneus,what = "std.all")$beta[1,2]

##entorhinal
Model_CT_entorhinal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                        Cognitive_factor ~ entorhinal_ct'
fit_sem_CT_entorhinal <- sem(Model_CT_entorhinal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_entorhinal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[5,2]<-BIC(fit_sem_CT_entorhinal)
data_subsample_GMlat_15_GGSEG[5,3]<-parameterEstimates(fit_sem_CT_entorhinal)[6,7]
data_subsample_GMlat_15_GGSEG[5,4]<-lavInspect(fit_sem_CT_entorhinal,what = "std.all")$beta[1,2]

##fusiform
Model_CT_fusiform <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                      Cognitive_factor ~ fusiform_ct'
fit_sem_CT_fusiform <- sem(Model_CT_fusiform, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_fusiform, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[6,2]<-BIC(fit_sem_CT_fusiform)
data_subsample_GMlat_15_GGSEG[6,3]<-parameterEstimates(fit_sem_CT_fusiform)[6,7]
data_subsample_GMlat_15_GGSEG[6,4]<-lavInspect(fit_sem_CT_fusiform,what = "std.all")$beta[1,2]

##inferiorparietal
Model_CT_inferiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ inferiorparietal_ct'
fit_sem_CT_inferiorparietal <- sem(Model_CT_inferiorparietal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_inferiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[7,2]<-BIC(fit_sem_CT_inferiorparietal)
data_subsample_GMlat_15_GGSEG[7,3]<-parameterEstimates(fit_sem_CT_inferiorparietal)[6,7]
data_subsample_GMlat_15_GGSEG[7,4]<-lavInspect(fit_sem_CT_inferiorparietal,what = "std.all")$beta[1,2]

##inferiortemporal
Model_CT_inferiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ inferiortemporal_ct'
fit_sem_CT_inferiortemporal <- sem(Model_CT_inferiortemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_inferiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[8,2]<-BIC(fit_sem_CT_inferiortemporal)
data_subsample_GMlat_15_GGSEG[8,3]<-parameterEstimates(fit_sem_CT_inferiortemporal)[6,7]
data_subsample_GMlat_15_GGSEG[8,4]<-lavInspect(fit_sem_CT_inferiortemporal,what = "std.all")$beta[1,2]

##isthmuscingulate
Model_CT_isthmuscingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ isthmuscingulate_ct'
fit_sem_CT_isthmuscingulate <- sem(Model_CT_isthmuscingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_isthmuscingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[9,2]<-BIC(fit_sem_CT_isthmuscingulate)
data_subsample_GMlat_15_GGSEG[9,3]<-parameterEstimates(fit_sem_CT_isthmuscingulate)[6,7]
data_subsample_GMlat_15_GGSEG[9,4]<-lavInspect(fit_sem_CT_isthmuscingulate,what = "std.all")$beta[1,2]

##lateraloccipital
Model_CT_lateraloccipital <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ lateraloccipital_ct'
fit_sem_CT_lateraloccipital <- sem(Model_CT_lateraloccipital, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_lateraloccipital, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[10,2]<-BIC(fit_sem_CT_lateraloccipital)
data_subsample_GMlat_15_GGSEG[10,3]<-parameterEstimates(fit_sem_CT_lateraloccipital)[6,7]
data_subsample_GMlat_15_GGSEG[10,4]<-lavInspect(fit_sem_CT_lateraloccipital,what = "std.all")$beta[1,2]

##lateralorbitofrontal
Model_CT_lateralorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                  Cognitive_factor ~ lateralorbitofrontal_ct'
fit_sem_CT_lateralorbitofrontal <- sem(Model_CT_lateralorbitofrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_lateralorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[11,2]<-BIC(fit_sem_CT_lateralorbitofrontal)
data_subsample_GMlat_15_GGSEG[11,3]<-parameterEstimates(fit_sem_CT_lateralorbitofrontal)[6,7]
data_subsample_GMlat_15_GGSEG[11,4]<-lavInspect(fit_sem_CT_lateralorbitofrontal,what = "std.all")$beta[1,2]

##lingual
Model_CT_lingual <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                     Cognitive_factor ~ lingual_ct'
fit_sem_CT_lingual <- sem(Model_CT_lingual, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_lingual, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[12,2]<-BIC(fit_sem_CT_lingual)
data_subsample_GMlat_15_GGSEG[12,3]<-parameterEstimates(fit_sem_CT_lingual)[6,7]
data_subsample_GMlat_15_GGSEG[12,4]<-lavInspect(fit_sem_CT_lingual,what = "std.all")$beta[1,2]

##medialorbitofrontal
Model_CT_medialorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                 Cognitive_factor ~ medialorbitofrontal_ct'
fit_sem_CT_medialorbitofrontal <- sem(Model_CT_medialorbitofrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_medialorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[13,2]<-BIC(fit_sem_CT_medialorbitofrontal)
data_subsample_GMlat_15_GGSEG[13,3]<-parameterEstimates(fit_sem_CT_medialorbitofrontal)[6,7]
data_subsample_GMlat_15_GGSEG[13,4]<-lavInspect(fit_sem_CT_medialorbitofrontal,what = "std.all")$beta[1,2]

##middletemporal
Model_CT_middletemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ middletemporal_ct'
fit_sem_CT_middletemporal <- sem(Model_CT_middletemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_middletemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[14,2]<-BIC(fit_sem_CT_middletemporal)
data_subsample_GMlat_15_GGSEG[14,3]<-parameterEstimates(fit_sem_CT_middletemporal)[6,7]
data_subsample_GMlat_15_GGSEG[14,4]<-lavInspect(fit_sem_CT_middletemporal,what = "std.all")$beta[1,2]

##parahippocampal
Model_CT_parahippocampal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ parahippocampal_ct'
fit_sem_CT_parahippocampal <- sem(Model_CT_parahippocampal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_parahippocampal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[15,2]<-BIC(fit_sem_CT_parahippocampal)
data_subsample_GMlat_15_GGSEG[15,3]<-parameterEstimates(fit_sem_CT_parahippocampal)[6,7]
data_subsample_GMlat_15_GGSEG[15,4]<-lavInspect(fit_sem_CT_parahippocampal,what = "std.all")$beta[1,2]

##paracentral
Model_CT_paracentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ paracentral_ct'
fit_sem_CT_paracentral <- sem(Model_CT_paracentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_paracentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[16,2]<-BIC(fit_sem_CT_paracentral)
data_subsample_GMlat_15_GGSEG[16,3]<-parameterEstimates(fit_sem_CT_paracentral)[6,7]
data_subsample_GMlat_15_GGSEG[16,4]<-lavInspect(fit_sem_CT_paracentral,what = "std.all")$beta[1,2]

##parsopercularis
Model_CT_parsopercularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ parsopercularis_ct'
fit_sem_CT_parsopercularis <- sem(Model_CT_parsopercularis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_parsopercularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[17,2]<-BIC(fit_sem_CT_parsopercularis)
data_subsample_GMlat_15_GGSEG[17,3]<-parameterEstimates(fit_sem_CT_parsopercularis)[6,7]
data_subsample_GMlat_15_GGSEG[17,4]<-lavInspect(fit_sem_CT_parsopercularis,what = "std.all")$beta[1,2]

##parsorbitalis
Model_CT_parsorbitalis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ parsorbitalis_ct'
fit_sem_CT_parsorbitalis <- sem(Model_CT_parsorbitalis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_parsorbitalis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[18,2]<-BIC(fit_sem_CT_parsorbitalis)
data_subsample_GMlat_15_GGSEG[18,3]<-parameterEstimates(fit_sem_CT_parsorbitalis)[6,7]
data_subsample_GMlat_15_GGSEG[18,4]<-lavInspect(fit_sem_CT_parsorbitalis,what = "std.all")$beta[1,2]

##parstriangularis
Model_CT_parstriangularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ parstriangularis_ct'
fit_sem_CT_parstriangularis <- sem(Model_CT_parstriangularis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_parstriangularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[19,2]<-BIC(fit_sem_CT_parstriangularis)
data_subsample_GMlat_15_GGSEG[19,3]<-parameterEstimates(fit_sem_CT_parstriangularis)[6,7]
data_subsample_GMlat_15_GGSEG[19,4]<-lavInspect(fit_sem_CT_parstriangularis,what = "std.all")$beta[1,2]

##pericalcarine
Model_CT_pericalcarine <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ pericalcarine_ct'
fit_sem_CT_pericalcarine <- sem(Model_CT_pericalcarine, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_pericalcarine, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[20,2]<-BIC(fit_sem_CT_pericalcarine)
data_subsample_GMlat_15_GGSEG[20,3]<-parameterEstimates(fit_sem_CT_pericalcarine)[6,7]
data_subsample_GMlat_15_GGSEG[20,4]<-lavInspect(fit_sem_CT_pericalcarine,what = "std.all")$beta[1,2]

##postcentral
Model_CT_postcentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ postcentral_ct'
fit_sem_CT_postcentral <- sem(Model_CT_postcentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_postcentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[21,2]<-BIC(fit_sem_CT_postcentral)
data_subsample_GMlat_15_GGSEG[21,3]<-parameterEstimates(fit_sem_CT_postcentral)[6,7]
data_subsample_GMlat_15_GGSEG[21,4]<-lavInspect(fit_sem_CT_postcentral,what = "std.all")$beta[1,2]

##posteriorcingulate
Model_CT_posteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ posteriorcingulate_ct'
fit_sem_CT_posteriorcingulate <- sem(Model_CT_posteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_posteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[22,2]<-BIC(fit_sem_CT_posteriorcingulate)
data_subsample_GMlat_15_GGSEG[22,3]<-parameterEstimates(fit_sem_CT_posteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[22,4]<-lavInspect(fit_sem_CT_posteriorcingulate,what = "std.all")$beta[1,2]

##precentral
Model_CT_precentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                        Cognitive_factor ~ precentral_ct'
fit_sem_CT_precentral <- sem(Model_CT_precentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_precentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[23,2]<-BIC(fit_sem_CT_precentral)
data_subsample_GMlat_15_GGSEG[23,3]<-parameterEstimates(fit_sem_CT_precentral)[6,7]
data_subsample_GMlat_15_GGSEG[23,4]<-lavInspect(fit_sem_CT_precentral,what = "std.all")$beta[1,2]

##precuneus
Model_CT_precuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                       Cognitive_factor ~ precuneus_ct'
fit_sem_CT_precuneus <- sem(Model_CT_precuneus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_precuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[24,2]<-BIC(fit_sem_CT_precuneus)
data_subsample_GMlat_15_GGSEG[24,3]<-parameterEstimates(fit_sem_CT_precuneus)[6,7]
data_subsample_GMlat_15_GGSEG[24,4]<-lavInspect(fit_sem_CT_precuneus,what = "std.all")$beta[1,2]

##rostralanteriorcingulate
Model_CT_rostralanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                      Cognitive_factor ~ rostralanteriorcingulate_ct'
fit_sem_CT_rostralanteriorcingulate <- sem(Model_CT_rostralanteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_rostralanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[25,2]<-BIC(fit_sem_CT_rostralanteriorcingulate)
data_subsample_GMlat_15_GGSEG[25,3]<-parameterEstimates(fit_sem_CT_rostralanteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[25,4]<-lavInspect(fit_sem_CT_rostralanteriorcingulate,what = "std.all")$beta[1,2]

##rostralmiddlefrontal
Model_CT_rostralmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                  Cognitive_factor ~ rostralmiddlefrontal_ct'
fit_sem_CT_rostralmiddlefrontal <- sem(Model_CT_rostralmiddlefrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_rostralmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[26,2]<-BIC(fit_sem_CT_rostralmiddlefrontal)
data_subsample_GMlat_15_GGSEG[26,3]<-parameterEstimates(fit_sem_CT_rostralmiddlefrontal)[6,7]
data_subsample_GMlat_15_GGSEG[26,4]<-lavInspect(fit_sem_CT_rostralmiddlefrontal,what = "std.all")$beta[1,2]

##superiorfrontal
Model_CT_superiorfrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ superiorfrontal_ct'
fit_sem_CT_superiorfrontal <- sem(Model_CT_superiorfrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_superiorfrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[27,2]<-BIC(fit_sem_CT_superiorfrontal)
data_subsample_GMlat_15_GGSEG[27,3]<-parameterEstimates(fit_sem_CT_superiorfrontal)[6,7]
data_subsample_GMlat_15_GGSEG[27,4]<-lavInspect(fit_sem_CT_superiorfrontal,what = "std.all")$beta[1,2]

##superiorparietal
Model_CT_superiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ superiorparietal_ct'
fit_sem_CT_superiorparietal <- sem(Model_CT_superiorparietal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_superiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[28,2]<-BIC(fit_sem_CT_superiorparietal)
data_subsample_GMlat_15_GGSEG[28,3]<-parameterEstimates(fit_sem_CT_superiorparietal)[6,7]
data_subsample_GMlat_15_GGSEG[28,4]<-lavInspect(fit_sem_CT_superiorparietal,what = "std.all")$beta[1,2]

##superiortemporal
Model_CT_superiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ superiortemporal_ct'
fit_sem_CT_superiortemporal <- sem(Model_CT_superiortemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_superiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[29,2]<-BIC(fit_sem_CT_superiortemporal)
data_subsample_GMlat_15_GGSEG[29,3]<-parameterEstimates(fit_sem_CT_superiortemporal)[6,7]
data_subsample_GMlat_15_GGSEG[29,4]<-lavInspect(fit_sem_CT_superiortemporal,what = "std.all")$beta[1,2]

##supramarginal
Model_CT_supramarginal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ supramarginal_ct'
fit_sem_CT_supramarginal <- sem(Model_CT_supramarginal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_supramarginal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[30,2]<-BIC(fit_sem_CT_supramarginal)
data_subsample_GMlat_15_GGSEG[30,3]<-parameterEstimates(fit_sem_CT_supramarginal)[6,7]
data_subsample_GMlat_15_GGSEG[30,4]<-lavInspect(fit_sem_CT_supramarginal,what = "std.all")$beta[1,2]

##frontalpole
Model_CT_frontalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ frontalpole_ct'
fit_sem_CT_frontalpole <- sem(Model_CT_frontalpole, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_frontalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[31,2]<-BIC(fit_sem_CT_frontalpole)
data_subsample_GMlat_15_GGSEG[31,3]<-parameterEstimates(fit_sem_CT_frontalpole)[6,7]
data_subsample_GMlat_15_GGSEG[31,4]<-lavInspect(fit_sem_CT_frontalpole,what = "std.all")$beta[1,2]

##temporalpole
Model_CT_temporalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalpole_ct'
fit_sem_CT_temporalpole <- sem(Model_CT_temporalpole, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_temporalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[32,2]<-BIC(fit_sem_CT_temporalpole)
data_subsample_GMlat_15_GGSEG[32,3]<-parameterEstimates(fit_sem_CT_temporalpole)[6,7]
data_subsample_GMlat_15_GGSEG[32,4]<-lavInspect(fit_sem_CT_temporalpole,what = "std.all")$beta[1,2]

##transversetemporal
Model_CT_transversetemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ transversetemporal_ct'
fit_sem_CT_transversetemporal <- sem(Model_CT_transversetemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_transversetemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[33,2]<-BIC(fit_sem_CT_transversetemporal)
data_subsample_GMlat_15_GGSEG[33,3]<-parameterEstimates(fit_sem_CT_transversetemporal)[6,7]
data_subsample_GMlat_15_GGSEG[33,4]<-lavInspect(fit_sem_CT_transversetemporal,what = "std.all")$beta[1,2]

##insula
Model_CT_insula <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                    Cognitive_factor ~ insula_ct'
fit_sem_CT_insula <- sem(Model_CT_insula, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CT_insula, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[34,2]<-BIC(fit_sem_CT_insula)
data_subsample_GMlat_15_GGSEG[34,3]<-parameterEstimates(fit_sem_CT_insula)[6,7]
data_subsample_GMlat_15_GGSEG[34,4]<-lavInspect(fit_sem_CT_insula,what = "std.all")$beta[1,2]

#Surface Area (34 regions)
##bankssts
Model_SA_bankssts <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                      Cognitive_factor ~ bankssts_sa'
fit_sem_SA_bankssts <- sem(Model_SA_bankssts, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_bankssts, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[1,5]<-BIC(fit_sem_SA_bankssts)
data_subsample_GMlat_15_GGSEG[1,6]<-parameterEstimates(fit_sem_SA_bankssts)[6,7]
data_subsample_GMlat_15_GGSEG[1,7]<-lavInspect(fit_sem_SA_bankssts,what = "std.all")$beta[1,2]

##caudalanteriorcingulate
Model_SA_caudalanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ caudalanteriorcingulate_sa'
fit_sem_SA_caudalanteriorcingulate <- sem(Model_SA_caudalanteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_caudalanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[2,5]<-BIC(fit_sem_SA_caudalanteriorcingulate)
data_subsample_GMlat_15_GGSEG[2,6]<-parameterEstimates(fit_sem_SA_caudalanteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[2,7]<-lavInspect(fit_sem_SA_caudalanteriorcingulate,what = "std.all")$beta[1,2]

##caudalmiddlefrontal
Model_SA_caudalmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                 Cognitive_factor ~ caudalmiddlefrontal_sa'
fit_sem_SA_caudalmiddlefrontal <- sem(Model_SA_caudalmiddlefrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_caudalmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[3,5]<-BIC(fit_sem_SA_caudalmiddlefrontal)
data_subsample_GMlat_15_GGSEG[3,6]<-parameterEstimates(fit_sem_SA_caudalmiddlefrontal)[6,7]
data_subsample_GMlat_15_GGSEG[3,7]<-lavInspect(fit_sem_SA_caudalmiddlefrontal,what = "std.all")$beta[1,2]

##cuneus
Model_SA_cuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                    Cognitive_factor ~ cuneus_sa'
fit_sem_SA_cuneus <- sem(Model_SA_cuneus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_cuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[4,5]<-BIC(fit_sem_SA_cuneus)
data_subsample_GMlat_15_GGSEG[4,6]<-parameterEstimates(fit_sem_SA_cuneus)[6,7]
data_subsample_GMlat_15_GGSEG[4,7]<-lavInspect(fit_sem_SA_cuneus,what = "std.all")$beta[1,2]

##entorhinal
Model_SA_entorhinal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                        Cognitive_factor ~ entorhinal_sa'
fit_sem_SA_entorhinal <- sem(Model_SA_entorhinal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_entorhinal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[5,5]<-BIC(fit_sem_SA_entorhinal)
data_subsample_GMlat_15_GGSEG[5,6]<-parameterEstimates(fit_sem_SA_entorhinal)[6,7]
data_subsample_GMlat_15_GGSEG[5,7]<-lavInspect(fit_sem_SA_entorhinal,what = "std.all")$beta[1,2]

##fusiform
Model_SA_fusiform <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                      Cognitive_factor ~ fusiform_sa'
fit_sem_SA_fusiform <- sem(Model_SA_fusiform, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_fusiform, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[6,5]<-BIC(fit_sem_SA_fusiform)
data_subsample_GMlat_15_GGSEG[6,6]<-parameterEstimates(fit_sem_SA_fusiform)[6,7]
data_subsample_GMlat_15_GGSEG[6,7]<-lavInspect(fit_sem_SA_fusiform,what = "std.all")$beta[1,2]

##inferiorparietal
Model_SA_inferiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ inferiorparietal_sa'
fit_sem_SA_inferiorparietal <- sem(Model_SA_inferiorparietal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_inferiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[7,5]<-BIC(fit_sem_SA_inferiorparietal)
data_subsample_GMlat_15_GGSEG[7,6]<-parameterEstimates(fit_sem_SA_inferiorparietal)[6,7]
data_subsample_GMlat_15_GGSEG[7,7]<-lavInspect(fit_sem_SA_inferiorparietal,what = "std.all")$beta[1,2]

##inferiortemporal
Model_SA_inferiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ inferiortemporal_sa'
fit_sem_SA_inferiortemporal <- sem(Model_SA_inferiortemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_inferiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[8,5]<-BIC(fit_sem_SA_inferiortemporal)
data_subsample_GMlat_15_GGSEG[8,6]<-parameterEstimates(fit_sem_SA_inferiortemporal)[6,7]
data_subsample_GMlat_15_GGSEG[8,7]<-lavInspect(fit_sem_SA_inferiortemporal,what = "std.all")$beta[1,2]

##isthmuscingulate
Model_SA_isthmuscingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ isthmuscingulate_sa'
fit_sem_SA_isthmuscingulate <- sem(Model_SA_isthmuscingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_isthmuscingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[9,5]<-BIC(fit_sem_SA_isthmuscingulate)
data_subsample_GMlat_15_GGSEG[9,6]<-parameterEstimates(fit_sem_SA_isthmuscingulate)[6,7]
data_subsample_GMlat_15_GGSEG[9,7]<-lavInspect(fit_sem_SA_isthmuscingulate,what = "std.all")$beta[1,2]

##lateraloccipital
Model_SA_lateraloccipital <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ lateraloccipital_sa'
fit_sem_SA_lateraloccipital <- sem(Model_SA_lateraloccipital, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_lateraloccipital, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[10,5]<-BIC(fit_sem_SA_lateraloccipital)
data_subsample_GMlat_15_GGSEG[10,6]<-parameterEstimates(fit_sem_SA_lateraloccipital)[6,7]
data_subsample_GMlat_15_GGSEG[10,7]<-lavInspect(fit_sem_SA_lateraloccipital,what = "std.all")$beta[1,2]

##lateralorbitofrontal
Model_SA_lateralorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                  Cognitive_factor ~ lateralorbitofrontal_sa'
fit_sem_SA_lateralorbitofrontal <- sem(Model_SA_lateralorbitofrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_lateralorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[11,5]<-BIC(fit_sem_SA_lateralorbitofrontal)
data_subsample_GMlat_15_GGSEG[11,6]<-parameterEstimates(fit_sem_SA_lateralorbitofrontal)[6,7]
data_subsample_GMlat_15_GGSEG[11,7]<-lavInspect(fit_sem_SA_lateralorbitofrontal,what = "std.all")$beta[1,2]

##lingual
Model_SA_lingual <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                     Cognitive_factor ~ lingual_sa'
fit_sem_SA_lingual <- sem(Model_SA_lingual, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_lingual, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[12,5]<-BIC(fit_sem_SA_lingual)
data_subsample_GMlat_15_GGSEG[12,6]<-parameterEstimates(fit_sem_SA_lingual)[6,7]
data_subsample_GMlat_15_GGSEG[12,7]<-lavInspect(fit_sem_SA_lingual,what = "std.all")$beta[1,2]

##medialorbitofrontal
Model_SA_medialorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                 Cognitive_factor ~ medialorbitofrontal_sa'
fit_sem_SA_medialorbitofrontal <- sem(Model_SA_medialorbitofrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_medialorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[13,5]<-BIC(fit_sem_SA_medialorbitofrontal)
data_subsample_GMlat_15_GGSEG[13,6]<-parameterEstimates(fit_sem_SA_medialorbitofrontal)[6,7]
data_subsample_GMlat_15_GGSEG[13,7]<-lavInspect(fit_sem_SA_medialorbitofrontal,what = "std.all")$beta[1,2]

##middletemporal
Model_SA_middletemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ middletemporal_sa'
fit_sem_SA_middletemporal <- sem(Model_SA_middletemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_middletemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[14,5]<-BIC(fit_sem_SA_middletemporal)
data_subsample_GMlat_15_GGSEG[14,6]<-parameterEstimates(fit_sem_SA_middletemporal)[6,7]
data_subsample_GMlat_15_GGSEG[14,7]<-lavInspect(fit_sem_SA_middletemporal,what = "std.all")$beta[1,2]

##parahippocampal
Model_SA_parahippocampal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ parahippocampal_sa'
fit_sem_SA_parahippocampal <- sem(Model_SA_parahippocampal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_parahippocampal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[15,5]<-BIC(fit_sem_SA_parahippocampal)
data_subsample_GMlat_15_GGSEG[15,6]<-parameterEstimates(fit_sem_SA_parahippocampal)[6,7]
data_subsample_GMlat_15_GGSEG[15,7]<-lavInspect(fit_sem_SA_parahippocampal,what = "std.all")$beta[1,2]

##paracentral
Model_SA_paracentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ paracentral_sa'
fit_sem_SA_paracentral <- sem(Model_SA_paracentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_paracentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[16,5]<-BIC(fit_sem_SA_paracentral)
data_subsample_GMlat_15_GGSEG[16,6]<-parameterEstimates(fit_sem_SA_paracentral)[6,7]
data_subsample_GMlat_15_GGSEG[16,7]<-lavInspect(fit_sem_SA_paracentral,what = "std.all")$beta[1,2]

##parsopercularis
Model_SA_parsopercularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ parsopercularis_sa'
fit_sem_SA_parsopercularis <- sem(Model_SA_parsopercularis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_parsopercularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[17,5]<-BIC(fit_sem_SA_parsopercularis)
data_subsample_GMlat_15_GGSEG[17,6]<-parameterEstimates(fit_sem_SA_parsopercularis)[6,7]
data_subsample_GMlat_15_GGSEG[17,7]<-lavInspect(fit_sem_SA_parsopercularis,what = "std.all")$beta[1,2]

##parsorbitalis
Model_SA_parsorbitalis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ parsorbitalis_sa'
fit_sem_SA_parsorbitalis <- sem(Model_SA_parsorbitalis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_parsorbitalis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[18,5]<-BIC(fit_sem_SA_parsorbitalis)
data_subsample_GMlat_15_GGSEG[18,6]<-parameterEstimates(fit_sem_SA_parsorbitalis)[6,7]
data_subsample_GMlat_15_GGSEG[18,7]<-lavInspect(fit_sem_SA_parsorbitalis,what = "std.all")$beta[1,2]

##parstriangularis
Model_SA_parstriangularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ parstriangularis_sa'
fit_sem_SA_parstriangularis <- sem(Model_SA_parstriangularis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_parstriangularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[19,5]<-BIC(fit_sem_SA_parstriangularis)
data_subsample_GMlat_15_GGSEG[19,6]<-parameterEstimates(fit_sem_SA_parstriangularis)[6,7]
data_subsample_GMlat_15_GGSEG[19,7]<-lavInspect(fit_sem_SA_parstriangularis,what = "std.all")$beta[1,2]

##pericalcarine
Model_SA_pericalcarine <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ pericalcarine_sa'
fit_sem_SA_pericalcarine <- sem(Model_SA_pericalcarine, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_pericalcarine, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[20,5]<-BIC(fit_sem_SA_pericalcarine)
data_subsample_GMlat_15_GGSEG[20,6]<-parameterEstimates(fit_sem_SA_pericalcarine)[6,7]
data_subsample_GMlat_15_GGSEG[20,7]<-lavInspect(fit_sem_SA_pericalcarine,what = "std.all")$beta[1,2]

##postcentral
Model_SA_postcentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ postcentral_sa'
fit_sem_SA_postcentral <- sem(Model_SA_postcentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_postcentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[21,5]<-BIC(fit_sem_SA_postcentral)
data_subsample_GMlat_15_GGSEG[21,6]<-parameterEstimates(fit_sem_SA_postcentral)[6,7]
data_subsample_GMlat_15_GGSEG[21,7]<-lavInspect(fit_sem_SA_postcentral,what = "std.all")$beta[1,2]

##posteriorcingulate
Model_SA_posteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ posteriorcingulate_sa'
fit_sem_SA_posteriorcingulate <- sem(Model_SA_posteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_posteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[22,5]<-BIC(fit_sem_SA_posteriorcingulate)
data_subsample_GMlat_15_GGSEG[22,6]<-parameterEstimates(fit_sem_SA_posteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[22,7]<-lavInspect(fit_sem_SA_posteriorcingulate,what = "std.all")$beta[1,2]

##precentral
Model_SA_precentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                        Cognitive_factor ~ precentral_sa'
fit_sem_SA_precentral <- sem(Model_SA_precentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_precentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[23,5]<-BIC(fit_sem_SA_precentral)
data_subsample_GMlat_15_GGSEG[23,6]<-parameterEstimates(fit_sem_SA_precentral)[6,7]
data_subsample_GMlat_15_GGSEG[23,7]<-lavInspect(fit_sem_SA_precentral,what = "std.all")$beta[1,2]

##precuneus
Model_SA_precuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                       Cognitive_factor ~ precuneus_sa'
fit_sem_SA_precuneus <- sem(Model_SA_precuneus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_precuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[24,5]<-BIC(fit_sem_SA_precuneus)
data_subsample_GMlat_15_GGSEG[24,6]<-parameterEstimates(fit_sem_SA_precuneus)[6,7]
data_subsample_GMlat_15_GGSEG[24,7]<-lavInspect(fit_sem_SA_precuneus,what = "std.all")$beta[1,2]

##rostralanteriorcingulate
Model_SA_rostralanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                      Cognitive_factor ~ rostralanteriorcingulate_sa'
fit_sem_SA_rostralanteriorcingulate <- sem(Model_SA_rostralanteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_rostralanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[25,5]<-BIC(fit_sem_SA_rostralanteriorcingulate)
data_subsample_GMlat_15_GGSEG[25,6]<-parameterEstimates(fit_sem_SA_rostralanteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[25,7]<-lavInspect(fit_sem_SA_rostralanteriorcingulate,what = "std.all")$beta[1,2]

##rostralmiddlefrontal
Model_SA_rostralmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                  Cognitive_factor ~ rostralmiddlefrontal_sa'
fit_sem_SA_rostralmiddlefrontal <- sem(Model_SA_rostralmiddlefrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_rostralmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[26,5]<-BIC(fit_sem_SA_rostralmiddlefrontal)
data_subsample_GMlat_15_GGSEG[26,6]<-parameterEstimates(fit_sem_SA_rostralmiddlefrontal)[6,7]
data_subsample_GMlat_15_GGSEG[26,7]<-lavInspect(fit_sem_SA_rostralmiddlefrontal,what = "std.all")$beta[1,2]

##superiorfrontal
Model_SA_superiorfrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ superiorfrontal_sa'
fit_sem_SA_superiorfrontal <- sem(Model_SA_superiorfrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_superiorfrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[27,5]<-BIC(fit_sem_SA_superiorfrontal)
data_subsample_GMlat_15_GGSEG[27,6]<-parameterEstimates(fit_sem_SA_superiorfrontal)[6,7]
data_subsample_GMlat_15_GGSEG[27,7]<-lavInspect(fit_sem_SA_superiorfrontal,what = "std.all")$beta[1,2]

##superiorparietal
Model_SA_superiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ superiorparietal_sa'
fit_sem_SA_superiorparietal <- sem(Model_SA_superiorparietal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_superiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[28,5]<-BIC(fit_sem_SA_superiorparietal)
data_subsample_GMlat_15_GGSEG[28,6]<-parameterEstimates(fit_sem_SA_superiorparietal)[6,7]
data_subsample_GMlat_15_GGSEG[28,7]<-lavInspect(fit_sem_SA_superiorparietal,what = "std.all")$beta[1,2]

##superiortemporal
Model_SA_superiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ superiortemporal_sa'
fit_sem_SA_superiortemporal <- sem(Model_SA_superiortemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_superiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[29,5]<-BIC(fit_sem_SA_superiortemporal)
data_subsample_GMlat_15_GGSEG[29,6]<-parameterEstimates(fit_sem_SA_superiortemporal)[6,7]
data_subsample_GMlat_15_GGSEG[29,7]<-lavInspect(fit_sem_SA_superiortemporal,what = "std.all")$beta[1,2]

##supramarginal
Model_SA_supramarginal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ supramarginal_sa'
fit_sem_SA_supramarginal <- sem(Model_SA_supramarginal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_supramarginal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[30,5]<-BIC(fit_sem_SA_supramarginal)
data_subsample_GMlat_15_GGSEG[30,6]<-parameterEstimates(fit_sem_SA_supramarginal)[6,7]
data_subsample_GMlat_15_GGSEG[30,7]<-lavInspect(fit_sem_SA_supramarginal,what = "std.all")$beta[1,2]

##frontalpole
Model_SA_frontalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ frontalpole_sa'
fit_sem_SA_frontalpole <- sem(Model_SA_frontalpole, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_frontalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[31,5]<-BIC(fit_sem_SA_frontalpole)
data_subsample_GMlat_15_GGSEG[31,6]<-parameterEstimates(fit_sem_SA_frontalpole)[6,7]
data_subsample_GMlat_15_GGSEG[31,7]<-lavInspect(fit_sem_SA_frontalpole,what = "std.all")$beta[1,2]

##temporalpole
Model_SA_temporalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalpole_sa'
fit_sem_SA_temporalpole <- sem(Model_SA_temporalpole, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_temporalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[32,5]<-BIC(fit_sem_SA_temporalpole)
data_subsample_GMlat_15_GGSEG[32,6]<-parameterEstimates(fit_sem_SA_temporalpole)[6,7]
data_subsample_GMlat_15_GGSEG[32,7]<-lavInspect(fit_sem_SA_temporalpole,what = "std.all")$beta[1,2]

##transversetemporal
Model_SA_transversetemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ transversetemporal_sa'
fit_sem_SA_transversetemporal <- sem(Model_SA_transversetemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_transversetemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[33,5]<-BIC(fit_sem_SA_transversetemporal)
data_subsample_GMlat_15_GGSEG[33,6]<-parameterEstimates(fit_sem_SA_transversetemporal)[6,7]
data_subsample_GMlat_15_GGSEG[33,7]<-lavInspect(fit_sem_SA_transversetemporal,what = "std.all")$beta[1,2]

##insula
Model_SA_insula <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                    Cognitive_factor ~ insula_sa'
fit_sem_SA_insula <- sem(Model_SA_insula, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SA_insula, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[34,5]<-BIC(fit_sem_SA_insula)
data_subsample_GMlat_15_GGSEG[34,6]<-parameterEstimates(fit_sem_SA_insula)[6,7]
data_subsample_GMlat_15_GGSEG[34,7]<-lavInspect(fit_sem_SA_insula,what = "std.all")$beta[1,2]

#Grey Matter Volume (34 regions)
##bankssts
Model_GMV_bankssts <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                       Cognitive_factor ~ bankssts_gmv'
fit_sem_GMV_bankssts <- sem(Model_GMV_bankssts, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_bankssts, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[1,8]<-BIC(fit_sem_GMV_bankssts)
data_subsample_GMlat_15_GGSEG[1,9]<-parameterEstimates(fit_sem_GMV_bankssts)[6,7]
data_subsample_GMlat_15_GGSEG[1,10]<-lavInspect(fit_sem_GMV_bankssts,what = "std.all")$beta[1,2]

##caudalanteriorcingulate
Model_GMV_caudalanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                      Cognitive_factor ~ caudalanteriorcingulate_gmv'
fit_sem_GMV_caudalanteriorcingulate <- sem(Model_GMV_caudalanteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_caudalanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[2,8]<-BIC(fit_sem_GMV_caudalanteriorcingulate)
data_subsample_GMlat_15_GGSEG[2,9]<-parameterEstimates(fit_sem_GMV_caudalanteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[2,10]<-lavInspect(fit_sem_GMV_caudalanteriorcingulate,what = "std.all")$beta[1,2]

##caudalmiddlefrontal
Model_GMV_caudalmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                  Cognitive_factor ~ caudalmiddlefrontal_gmv'
fit_sem_GMV_caudalmiddlefrontal <- sem(Model_GMV_caudalmiddlefrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_caudalmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[3,8]<-BIC(fit_sem_GMV_caudalmiddlefrontal)
data_subsample_GMlat_15_GGSEG[3,9]<-parameterEstimates(fit_sem_GMV_caudalmiddlefrontal)[6,7]
data_subsample_GMlat_15_GGSEG[3,10]<-lavInspect(fit_sem_GMV_caudalmiddlefrontal,what = "std.all")$beta[1,2]

##cuneus
Model_GMV_cuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                     Cognitive_factor ~ cuneus_gmv'
fit_sem_GMV_cuneus <- sem(Model_GMV_cuneus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_cuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[4,8]<-BIC(fit_sem_GMV_cuneus)
data_subsample_GMlat_15_GGSEG[4,9]<-parameterEstimates(fit_sem_GMV_cuneus)[6,7]
data_subsample_GMlat_15_GGSEG[4,10]<-lavInspect(fit_sem_GMV_cuneus,what = "std.all")$beta[1,2]

##entorhinal
Model_GMV_entorhinal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ entorhinal_gmv'
fit_sem_GMV_entorhinal <- sem(Model_GMV_entorhinal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_entorhinal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[5,8]<-BIC(fit_sem_GMV_entorhinal)
data_subsample_GMlat_15_GGSEG[5,9]<-parameterEstimates(fit_sem_GMV_entorhinal)[6,7]
data_subsample_GMlat_15_GGSEG[5,10]<-lavInspect(fit_sem_GMV_entorhinal,what = "std.all")$beta[1,2]

##fusiform
Model_GMV_fusiform <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                       Cognitive_factor ~ fusiform_gmv'
fit_sem_GMV_fusiform <- sem(Model_GMV_fusiform, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_fusiform, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[6,8]<-BIC(fit_sem_GMV_fusiform)
data_subsample_GMlat_15_GGSEG[6,9]<-parameterEstimates(fit_sem_GMV_fusiform)[6,7]
data_subsample_GMlat_15_GGSEG[6,10]<-lavInspect(fit_sem_GMV_fusiform,what = "std.all")$beta[1,2]

##inferiorparietal
Model_GMV_inferiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                               Cognitive_factor ~ inferiorparietal_gmv'
fit_sem_GMV_inferiorparietal <- sem(Model_GMV_inferiorparietal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_inferiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[7,8]<-BIC(fit_sem_GMV_inferiorparietal)
data_subsample_GMlat_15_GGSEG[7,9]<-parameterEstimates(fit_sem_GMV_inferiorparietal)[6,7]
data_subsample_GMlat_15_GGSEG[7,10]<-lavInspect(fit_sem_GMV_inferiorparietal,what = "std.all")$beta[1,2]

##inferiortemporal
Model_GMV_inferiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                               Cognitive_factor ~ inferiortemporal_gmv'
fit_sem_GMV_inferiortemporal <- sem(Model_GMV_inferiortemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_inferiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[8,8]<-BIC(fit_sem_GMV_inferiortemporal)
data_subsample_GMlat_15_GGSEG[8,9]<-parameterEstimates(fit_sem_GMV_inferiortemporal)[6,7]
data_subsample_GMlat_15_GGSEG[8,10]<-lavInspect(fit_sem_GMV_inferiortemporal,what = "std.all")$beta[1,2]

##isthmuscingulate
Model_GMV_isthmuscingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                               Cognitive_factor ~ isthmuscingulate_gmv'
fit_sem_GMV_isthmuscingulate <- sem(Model_GMV_isthmuscingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_isthmuscingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[9,8]<-BIC(fit_sem_GMV_isthmuscingulate)
data_subsample_GMlat_15_GGSEG[9,9]<-parameterEstimates(fit_sem_GMV_isthmuscingulate)[6,7]
data_subsample_GMlat_15_GGSEG[9,10]<-lavInspect(fit_sem_GMV_isthmuscingulate,what = "std.all")$beta[1,2]

##lateraloccipital
Model_GMV_lateraloccipital <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                               Cognitive_factor ~ lateraloccipital_gmv'
fit_sem_GMV_lateraloccipital <- sem(Model_GMV_lateraloccipital, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_lateraloccipital, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[10,8]<-BIC(fit_sem_GMV_lateraloccipital)
data_subsample_GMlat_15_GGSEG[10,9]<-parameterEstimates(fit_sem_GMV_lateraloccipital)[6,7]
data_subsample_GMlat_15_GGSEG[10,10]<-lavInspect(fit_sem_GMV_lateraloccipital,what = "std.all")$beta[1,2]

##lateralorbitofrontal
Model_GMV_lateralorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ lateralorbitofrontal_gmv'
fit_sem_GMV_lateralorbitofrontal <- sem(Model_GMV_lateralorbitofrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_lateralorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[11,8]<-BIC(fit_sem_GMV_lateralorbitofrontal)
data_subsample_GMlat_15_GGSEG[11,9]<-parameterEstimates(fit_sem_GMV_lateralorbitofrontal)[6,7]
data_subsample_GMlat_15_GGSEG[11,10]<-lavInspect(fit_sem_GMV_lateralorbitofrontal,what = "std.all")$beta[1,2]

##lingual
Model_GMV_lingual <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lingual_gmv'
fit_sem_GMV_lingual <- sem(Model_GMV_lingual, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_lingual, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[12,8]<-BIC(fit_sem_GMV_lingual)
data_subsample_GMlat_15_GGSEG[12,9]<-parameterEstimates(fit_sem_GMV_lingual)[6,7]
data_subsample_GMlat_15_GGSEG[12,10]<-lavInspect(fit_sem_GMV_lingual,what = "std.all")$beta[1,2]

##medialorbitofrontal
Model_GMV_medialorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ medialorbitofrontal_gmv'
fit_sem_GMV_medialorbitofrontal <- sem(Model_GMV_medialorbitofrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_medialorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[13,8]<-BIC(fit_sem_GMV_medialorbitofrontal)
data_subsample_GMlat_15_GGSEG[13,9]<-parameterEstimates(fit_sem_GMV_medialorbitofrontal)[6,7]
data_subsample_GMlat_15_GGSEG[13,10]<-lavInspect(fit_sem_GMV_medialorbitofrontal,what = "std.all")$beta[1,2]

##middletemporal
Model_GMV_middletemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ middletemporal_gmv'
fit_sem_GMV_middletemporal <- sem(Model_GMV_middletemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_middletemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[14,8]<-BIC(fit_sem_GMV_middletemporal)
data_subsample_GMlat_15_GGSEG[14,9]<-parameterEstimates(fit_sem_GMV_middletemporal)[6,7]
data_subsample_GMlat_15_GGSEG[14,10]<-lavInspect(fit_sem_GMV_middletemporal,what = "std.all")$beta[1,2]

##parahippocampal
Model_GMV_parahippocampal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parahippocampal_gmv'
fit_sem_GMV_parahippocampal <- sem(Model_GMV_parahippocampal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_parahippocampal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[15,8]<-BIC(fit_sem_GMV_parahippocampal)
data_subsample_GMlat_15_GGSEG[15,9]<-parameterEstimates(fit_sem_GMV_parahippocampal)[6,7]
data_subsample_GMlat_15_GGSEG[15,10]<-lavInspect(fit_sem_GMV_parahippocampal,what = "std.all")$beta[1,2]

##paracentral
Model_GMV_paracentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ paracentral_gmv'
fit_sem_GMV_paracentral <- sem(Model_GMV_paracentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_paracentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[16,8]<-BIC(fit_sem_GMV_paracentral)
data_subsample_GMlat_15_GGSEG[16,9]<-parameterEstimates(fit_sem_GMV_paracentral)[6,7]
data_subsample_GMlat_15_GGSEG[16,10]<-lavInspect(fit_sem_GMV_paracentral,what = "std.all")$beta[1,2]

##parsopercularis
Model_GMV_parsopercularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parsopercularis_gmv'
fit_sem_GMV_parsopercularis <- sem(Model_GMV_parsopercularis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_parsopercularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[17,8]<-BIC(fit_sem_GMV_parsopercularis)
data_subsample_GMlat_15_GGSEG[17,9]<-parameterEstimates(fit_sem_GMV_parsopercularis)[6,7]
data_subsample_GMlat_15_GGSEG[17,10]<-lavInspect(fit_sem_GMV_parsopercularis,what = "std.all")$beta[1,2]

##parsorbitalis
Model_GMV_parsorbitalis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parsorbitalis_gmv'
fit_sem_GMV_parsorbitalis <- sem(Model_GMV_parsorbitalis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_parsorbitalis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[18,8]<-BIC(fit_sem_GMV_parsorbitalis)
data_subsample_GMlat_15_GGSEG[18,9]<-parameterEstimates(fit_sem_GMV_parsorbitalis)[6,7]
data_subsample_GMlat_15_GGSEG[18,10]<-lavInspect(fit_sem_GMV_parsorbitalis,what = "std.all")$beta[1,2]

##parstriangularis
Model_GMV_parstriangularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parstriangularis_gmv'
fit_sem_GMV_parstriangularis <- sem(Model_GMV_parstriangularis, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_parstriangularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[19,8]<-BIC(fit_sem_GMV_parstriangularis)
data_subsample_GMlat_15_GGSEG[19,9]<-parameterEstimates(fit_sem_GMV_parstriangularis)[6,7]
data_subsample_GMlat_15_GGSEG[19,10]<-lavInspect(fit_sem_GMV_parstriangularis,what = "std.all")$beta[1,2]

##pericalcarine
Model_GMV_pericalcarine <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ pericalcarine_gmv'
fit_sem_GMV_pericalcarine <- sem(Model_GMV_pericalcarine, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_pericalcarine, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[20,8]<-BIC(fit_sem_GMV_pericalcarine)
data_subsample_GMlat_15_GGSEG[20,9]<-parameterEstimates(fit_sem_GMV_pericalcarine)[6,7]
data_subsample_GMlat_15_GGSEG[20,10]<-lavInspect(fit_sem_GMV_pericalcarine,what = "std.all")$beta[1,2]

##postcentral
Model_GMV_postcentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ postcentral_gmv'
fit_sem_GMV_postcentral <- sem(Model_GMV_postcentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_postcentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[21,8]<-BIC(fit_sem_GMV_postcentral)
data_subsample_GMlat_15_GGSEG[21,9]<-parameterEstimates(fit_sem_GMV_postcentral)[6,7]
data_subsample_GMlat_15_GGSEG[21,10]<-lavInspect(fit_sem_GMV_postcentral,what = "std.all")$beta[1,2]

##posteriorcingulate
Model_GMV_posteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ posteriorcingulate_gmv'
fit_sem_GMV_posteriorcingulate <- sem(Model_GMV_posteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_posteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[22,8]<-BIC(fit_sem_GMV_posteriorcingulate)
data_subsample_GMlat_15_GGSEG[22,9]<-parameterEstimates(fit_sem_GMV_posteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[22,10]<-lavInspect(fit_sem_GMV_posteriorcingulate,what = "std.all")$beta[1,2]

##precentral
Model_GMV_precentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ precentral_gmv'
fit_sem_GMV_precentral <- sem(Model_GMV_precentral, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_precentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[23,8]<-BIC(fit_sem_GMV_precentral)
data_subsample_GMlat_15_GGSEG[23,9]<-parameterEstimates(fit_sem_GMV_precentral)[6,7]
data_subsample_GMlat_15_GGSEG[23,10]<-lavInspect(fit_sem_GMV_precentral,what = "std.all")$beta[1,2]

##precuneus
Model_GMV_precuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ precuneus_gmv'
fit_sem_GMV_precuneus <- sem(Model_GMV_precuneus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_precuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[24,8]<-BIC(fit_sem_GMV_precuneus)
data_subsample_GMlat_15_GGSEG[24,9]<-parameterEstimates(fit_sem_GMV_precuneus)[6,7]
data_subsample_GMlat_15_GGSEG[24,10]<-lavInspect(fit_sem_GMV_precuneus,what = "std.all")$beta[1,2]

##rostralanteriorcingulate
Model_GMV_rostralanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ rostralanteriorcingulate_gmv'
fit_sem_GMV_rostralanteriorcingulate <- sem(Model_GMV_rostralanteriorcingulate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_rostralanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[25,8]<-BIC(fit_sem_GMV_rostralanteriorcingulate)
data_subsample_GMlat_15_GGSEG[25,9]<-parameterEstimates(fit_sem_GMV_rostralanteriorcingulate)[6,7]
data_subsample_GMlat_15_GGSEG[25,10]<-lavInspect(fit_sem_GMV_rostralanteriorcingulate,what = "std.all")$beta[1,2]

##rostralmiddlefrontal
Model_GMV_rostralmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ rostralmiddlefrontal_gmv'
fit_sem_GMV_rostralmiddlefrontal <- sem(Model_GMV_rostralmiddlefrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_rostralmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[26,8]<-BIC(fit_sem_GMV_rostralmiddlefrontal)
data_subsample_GMlat_15_GGSEG[26,9]<-parameterEstimates(fit_sem_GMV_rostralmiddlefrontal)[6,7]
data_subsample_GMlat_15_GGSEG[26,10]<-lavInspect(fit_sem_GMV_rostralmiddlefrontal,what = "std.all")$beta[1,2]

##superiorfrontal
Model_GMV_superiorfrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorfrontal_gmv'
fit_sem_GMV_superiorfrontal <- sem(Model_GMV_superiorfrontal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_superiorfrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[27,8]<-BIC(fit_sem_GMV_superiorfrontal)
data_subsample_GMlat_15_GGSEG[27,9]<-parameterEstimates(fit_sem_GMV_superiorfrontal)[6,7]
data_subsample_GMlat_15_GGSEG[27,10]<-lavInspect(fit_sem_GMV_superiorfrontal,what = "std.all")$beta[1,2]

##superiorparietal
Model_GMV_superiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorparietal_gmv'
fit_sem_GMV_superiorparietal <- sem(Model_GMV_superiorparietal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_superiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[28,8]<-BIC(fit_sem_GMV_superiorparietal)
data_subsample_GMlat_15_GGSEG[28,9]<-parameterEstimates(fit_sem_GMV_superiorparietal)[6,7]
data_subsample_GMlat_15_GGSEG[28,10]<-lavInspect(fit_sem_GMV_superiorparietal,what = "std.all")$beta[1,2]

##superiortemporal
Model_GMV_superiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiortemporal_gmv'
fit_sem_GMV_superiortemporal <- sem(Model_GMV_superiortemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_superiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[29,8]<-BIC(fit_sem_GMV_superiortemporal)
data_subsample_GMlat_15_GGSEG[29,9]<-parameterEstimates(fit_sem_GMV_superiortemporal)[6,7]
data_subsample_GMlat_15_GGSEG[29,10]<-lavInspect(fit_sem_GMV_superiortemporal,what = "std.all")$beta[1,2]

##supramarginal
Model_GMV_supramarginal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ supramarginal_gmv'
fit_sem_GMV_supramarginal <- sem(Model_GMV_supramarginal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_supramarginal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[30,8]<-BIC(fit_sem_GMV_supramarginal)
data_subsample_GMlat_15_GGSEG[30,9]<-parameterEstimates(fit_sem_GMV_supramarginal)[6,7]
data_subsample_GMlat_15_GGSEG[30,10]<-lavInspect(fit_sem_GMV_supramarginal,what = "std.all")$beta[1,2]

##frontalpole
Model_GMV_frontalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ frontalpole_gmv'
fit_sem_GMV_frontalpole <- sem(Model_GMV_frontalpole, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_frontalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[31,8]<-BIC(fit_sem_GMV_frontalpole)
data_subsample_GMlat_15_GGSEG[31,9]<-parameterEstimates(fit_sem_GMV_frontalpole)[6,7]
data_subsample_GMlat_15_GGSEG[31,10]<-lavInspect(fit_sem_GMV_frontalpole,what = "std.all")$beta[1,2]

##temporalpole
Model_GMV_temporalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalpole_gmv'
fit_sem_GMV_temporalpole <- sem(Model_GMV_temporalpole, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_temporalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[32,8]<-BIC(fit_sem_GMV_temporalpole)
data_subsample_GMlat_15_GGSEG[32,9]<-parameterEstimates(fit_sem_GMV_temporalpole)[6,7]
data_subsample_GMlat_15_GGSEG[32,10]<-lavInspect(fit_sem_GMV_temporalpole,what = "std.all")$beta[1,2]

##transversetemporal
Model_GMV_transversetemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ transversetemporal_gmv'
fit_sem_GMV_transversetemporal <- sem(Model_GMV_transversetemporal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_transversetemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[33,8]<-BIC(fit_sem_GMV_transversetemporal)
data_subsample_GMlat_15_GGSEG[33,9]<-parameterEstimates(fit_sem_GMV_transversetemporal)[6,7]
data_subsample_GMlat_15_GGSEG[33,10]<-lavInspect(fit_sem_GMV_transversetemporal,what = "std.all")$beta[1,2]

##insula
Model_GMV_insula <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ insula_gmv'
fit_sem_GMV_insula <- sem(Model_GMV_insula, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_insula, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_15_GGSEG[34,8]<-BIC(fit_sem_GMV_insula)
data_subsample_GMlat_15_GGSEG[34,9]<-parameterEstimates(fit_sem_GMV_insula)[6,7]
data_subsample_GMlat_15_GGSEG[34,10]<-lavInspect(fit_sem_GMV_insula,what = "std.all")$beta[1,2]

#Data Visualization for the individual models in the three grey matter metrics
data_subsample_GMlat_15_rl_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                 CT_fit=rep(NA, times=68),
                                                 CT_p=rep(NA, times=68),
                                                 CT_std=rep(NA, times=68),
                                                 SA_fit=rep(NA, times=68),
                                                 SA_p=rep(NA, times=68),
                                                 SA_std=rep(NA, times=68),
                                                 GMV_fit=rep(NA, times=68),
                                                 GMV_p=rep(NA, times=68),
                                                 GMV_std=rep(NA, times=68))

for (j in 2:10) {  #The plots are made with GGSEG which only work with right and left regions so we need to duplicate the value in the left and right regions
  a=1
  for (i in seq(from = 1, to = 68, by = 2)) {
    data_subsample_GMlat_15_rl_GGSEG[i,j]<-data_subsample_GMlat_15_GGSEG[a,j]
    data_subsample_GMlat_15_rl_GGSEG[i+1,j]<-data_subsample_GMlat_15_GGSEG[a,j]
    a=a+1
  }
}

data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_fit)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="Fit")
data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_fit)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="Fit")
data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_fit)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="Fit")
data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

#Select all the std for the grey matter metrics
data_subsample_GMlat_15_GGSEG_CT <- data_subsample_GMlat_15_GGSEG %>%
  dplyr::select(-c("CT_fit","CT_p","SA_fit","SA_p","SA_std","GMV_fit","GMV_p","GMV_std")) %>%
  rename(c("CT_std" = "Std")) %>%
  add_column(Metric = "CT")

data_subsample_GMlat_15_GGSEG_SA <- data_subsample_GMlat_15_GGSEG %>%
  dplyr::select(-c("CT_fit","CT_p","CT_std","SA_fit","SA_p","GMV_fit","GMV_p","GMV_std")) %>%
  rename(c("SA_std" = "Std")) %>%
  add_column(Metric = "SA")

data_subsample_GMlat_15_GGSEG_GMV <- data_subsample_GMlat_15_GGSEG %>%
  dplyr::select(-c("CT_fit","CT_p","CT_std","SA_fit","SA_p","SA_std","GMV_fit","GMV_p")) %>%
  rename(c("GMV_std" = "Std")) %>%
  add_column(Metric = "GMV")

#Compute the maximum and the minimum
min(data_subsample_GMlat_15_GGSEG_CT[,2])
max(data_subsample_GMlat_15_GGSEG_CT[,2])
min(data_subsample_GMlat_15_GGSEG_SA[,2])
max(data_subsample_GMlat_15_GGSEG_SA[,2])
min(data_subsample_GMlat_15_GGSEG_GMV[,2])
max(data_subsample_GMlat_15_GGSEG_GMV[,2])

data_subsample_STD_GMlat_15 <- rbind(data_subsample_GMlat_15_GGSEG_CT,data_subsample_GMlat_15_GGSEG_SA,by=c("label"))
data_subsample_STD_GMlat_15 <- rbind(data_subsample_STD_GMlat_15,data_subsample_GMlat_15_GGSEG_GMV,by=c("label"))
data_subsample_STD_GMlat_15 <- data_subsample_STD_GMlat_15[-c(69,104),] 

data_subsample_STD_GMlat_15 <- data_subsample_STD_GMlat_15 %>% group_by(Metric) %>%  mutate(mean = mean(as.numeric(Std)))

ggplot(data_subsample_STD_GMlat_15,aes(as.numeric(Std),fill=Metric)) + 
  geom_histogram(binwidth=0.005)+
  geom_density(adjust=2,alpha=0.2) +
  facet_grid(Metric~.) +
  geom_vline(aes(xintercept = mean, group = Metric), linetype="dashed", size=1,alpha=0.4) +
  theme_classic(base_size = 30)+
  xlab('Standardized estimated model parameters') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_fill_viridis(discrete=TRUE, option="rocket",begin=0.3,end=0.9,name="Metrics",labels=c("Cortical Thickness","Surface Area","Grey Matter Volume"))
ggsave("plot_GMlat_15_std.png")

#------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Model per metric - cortical thickness
#------------------------------------------------------------------------------------------------------------------------------
  
#Lasso regularization of the cortical thickness measures from all 34 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_ct <- data_total_subsample15[,c(10:43,172)]
df_ct <- makeX(df_ct, na.impute = TRUE)
x_ct <- df_ct[,1:34]
y_ct <- df_ct[,35]
cv_lasso_CTlat <- cv.glmnet(x_ct,y_ct, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_CTlat)
new_coef_ct <- coef(cv_lasso_CTlat, s = "lambda.min")
new_coef_ct <- as.data.frame(summary(new_coef_ct))
new_coef_ct <- new_coef_ct[,c(1,3)]
names(new_coef_ct) <- c("index","coef")

ROIs_coef_ct <- data.frame(label=c("bankssts_ct","caudalanteriorcingulate_ct","caudalmiddlefrontal_ct","cuneus_ct","entorhinal_ct","fusiform_ct","inferiorparietal_ct","inferiortemporal_ct","isthmuscingulate_ct","lateraloccipital_ct","lateralorbitofrontal_ct","lingual_ct","medialorbitofrontal_ct","middletemporal_ct","parahippocampal_ct","paracentral_ct","parsopercularis_ct","parsorbitalis_ct","parstriangularis_ct","pericalcarine_ct","postcentral_ct","posteriorcingulate_ct","precentral_ct","precuneus_ct","rostralanteriorcingulate_ct","rostralmiddlefrontal_ct","superiorfrontal_ct","superiorparietal_ct","superiortemporal_ct","supramarginal_ct","frontalpole_ct","temporalpole_ct","transversetemporal_ct","insula_ct"),
                           coef=rep("NA",times=34))

ROIs_coef_ct[new_coef_ct[1,1],2] <- new_coef_ct[1,2]

for (i in 1:nrow(new_coef_ct)) {
  ROIs_coef_ct[new_coef_ct[i,1],2] <- new_coef_ct[i,2]
}
ROI_ct<-c()
for (i in 1:34) {
  if (ROIs_coef_ct[i,2] != "NA") ROI_ct[i]<- ROIs_coef_ct[i,1]
}
ROI_ct<-na.omit(ROI_ct)

ROI_ct_plus<-c(ROI_ct[1])    #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_ct)) {
  ROI_ct_plus<-paste(ROI_ct_plus,"+",ROI_ct[i])
}
ROI_ct_plus_2<-ROI_ct_plus
cat(ROI_ct_plus_2)

#Regularized cortical thickness model
Model_CTlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + fusiform_ct + isthmuscingulate_ct + lateraloccipital_ct + lateralorbitofrontal_ct + medialorbitofrontal_ct + parahippocampal_ct + paracentral_ct + parsopercularis_ct + parsorbitalis_ct + pericalcarine_ct + postcentral_ct + posteriorcingulate_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiorparietal_ct + superiortemporal_ct + transversetemporal_ct + insula_ct'

fit_sem_CTlat_15_reg_free <- sem(Model_CTlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CTlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized cortical thickness model
Model_CTlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*fusiform_ct + a*isthmuscingulate_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*medialorbitofrontal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsopercularis_ct + a*parsorbitalis_ct + a*pericalcarine_ct + a*postcentral_ct + a*posteriorcingulate_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*superiortemporal_ct + a*transversetemporal_ct + a*insula_ct'

fit_sem_CTlat_15_reg_constrained <- sem(Model_CTlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CTlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
anova(fit_sem_CTlat_15_reg_free,fit_sem_CTlat_15_reg_constrained)

#Data Visualization for the regularized cortical thickness model
BIC(fit_sem_CTlat_15_reg_free)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
data_subsample_CTlat_15_rl_reg_GGSEG <- data.frame(label=c("lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_fusiform","rh_fusiform","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                   CT_p=rep(NA, times=44),
                                                   CT_std=rep(NA, times=44))
a=1
for (i in seq(from = 1, to = 44, by = 2)) {
  data_subsample_CTlat_15_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_CTlat_15_reg_free)[a+5,7]
  data_subsample_CTlat_15_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_CTlat_15_reg_free)[a+5,7]
  data_subsample_CTlat_15_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_CTlat_15_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_CTlat_15_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_CTlat_15_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_CTlat_15_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_CTlat_15_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
  

#Model for cortical thickness with all the regions of interest (34 ROIs)
Model_CTlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ bankssts_ct + caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + cuneus_ct + entorhinal_ct + fusiform_ct + inferiorparietal_ct + inferiortemporal_ct + isthmuscingulate_ct + lateraloccipital_ct + lateralorbitofrontal_ct + lingual_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsopercularis_ct + parsorbitalis_ct + parstriangularis_ct + pericalcarine_ct + postcentral_ct + posteriorcingulate_ct + precentral_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiorparietal_ct + superiortemporal_ct + supramarginal_ct + frontalpole_ct + temporalpole_ct + transversetemporal_ct + insula_ct'

fit_sem_CTlat_15_all_free <- sem(Model_CTlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CTlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the cortical thickness model with all the regions
Model_CTlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ a*bankssts_ct + a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*cuneus_ct + a*entorhinal_ct + a*fusiform_ct + a*inferiorparietal_ct + a*inferiortemporal_ct + a*isthmuscingulate_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*lingual_ct + a*medialorbitofrontal_ct + a*middletemporal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsopercularis_ct + a*parsorbitalis_ct + a*parstriangularis_ct + a*pericalcarine_ct + a*postcentral_ct + a*posteriorcingulate_ct + a*precentral_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*superiortemporal_ct + a*supramarginal_ct + a*frontalpole_ct + a*temporalpole_ct + a*transversetemporal_ct + a*insula_ct'

fit_sem_CTlat_15_all_constrained <- sem(Model_CTlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_CTlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
anova(fit_sem_CTlat_15_all_free,fit_sem_CTlat_15_all_constrained)

#Data Visualization for the cortical thickness model with all the regions
BIC(fit_sem_CTlat_15_all_free)
data_subsample_CTlat_15_rl_all_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                   CT_p=rep(NA, times=68),
                                                   CT_std=rep(NA, times=68))
a=1
for (i in seq(from = 1, to = 68, by = 2)) {
  data_subsample_CTlat_15_rl_all_GGSEG[i,2]<-parameterEstimates(fit_sem_CTlat_15_all_free)[a+5,7]
  data_subsample_CTlat_15_rl_all_GGSEG[i+1,2]<-parameterEstimates(fit_sem_CTlat_15_all_free)[a+5,7]
  data_subsample_CTlat_15_rl_all_GGSEG[i,3]<-lavInspect(fit_sem_CTlat_15_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_CTlat_15_rl_all_GGSEG[i+1,3]<-lavInspect(fit_sem_CTlat_15_all_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_CTlat_15_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_CTlat_15_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Model per metric - Surface Area
#------------------------------------------------------------------------------------------------------------------------------
  
#Lasso regularization of the surface area measures from all 34 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_sa <- data_total_subsample15[,c(44:77,172)]
df_sa <- makeX(df_sa, na.impute = TRUE)
x_sa <- df_sa[,1:34]
y_sa <- df_sa[,35]
cv_lasso_SAlat <- cv.glmnet(x_sa,y_sa, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_SAlat)
new_coef_sa <- coef(cv_lasso_SAlat, s = "lambda.min")
new_coef_sa <- as.data.frame(summary(new_coef_sa))
new_coef_sa <- new_coef_sa[,c(1,3)]
names(new_coef_sa) <- c("index","coef")

ROIs_coef_sa <- data.frame(label=c("bankssts_sa","caudalanteriorcingulate_sa","caudalmiddlefrontal_sa","cuneus_sa","entorhinal_sa","fusiform_sa","inferiorparietal_sa","inferiortemporal_sa","isthmuscingulate_sa","lateraloccipital_sa","lateralorbitofrontal_sa","lingual_sa","medialorbitofrontal_sa","middletemporal_sa","parahippocampal_sa","paracentral_sa","parsopercularis_sa","parsorbitalis_sa","parstriangularis_sa","pericalcarine_sa","postcentral_sa","posteriorcingulate_sa","precentral_sa","precuneus_sa","rostralanteriorcingulate_sa","rostralmiddlefrontal_sa","superiorfrontal_sa","superiorparietal_sa","superiortemporal_sa","supramarginal_sa","frontalpole_sa","temporalpole_sa","transversetemporal_sa","insula_sa"),
                           coef=rep("NA",times=34))

ROIs_coef_sa[new_coef_sa[1,1],2] <- new_coef_sa[1,2]

for (i in 1:nrow(new_coef_sa)) {
  ROIs_coef_sa[new_coef_sa[i,1],2] <- new_coef_sa[i,2]
}
ROI_sa<-c()
for (i in 1:34) {
  if (ROIs_coef_sa[i,2] != "NA") ROI_sa[i]<- ROIs_coef_sa[i,1]
}
ROI_sa<-na.omit(ROI_sa)

ROI_sa_plus<-c(ROI_sa[1])       #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_sa)) {
  ROI_sa_plus<-paste(ROI_sa_plus,"+",ROI_sa[i])
}
cat(ROI_sa_plus)
  
#Regularized surface area model
Model_SAlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ caudalanteriorcingulate_sa + cuneus_sa + entorhinal_sa + fusiform_sa + inferiorparietal_sa + isthmuscingulate_sa + lateraloccipital_sa + lingual_sa + parahippocampal_sa + paracentral_sa + parsopercularis_sa + parsorbitalis_sa + pericalcarine_sa + postcentral_sa + precentral_sa + precuneus_sa + rostralmiddlefrontal_sa + transversetemporal_sa'

fit_sem_SAlat_15_reg_free <- sem(Model_SAlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SAlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized surface area model
Model_SAlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*caudalanteriorcingulate_sa + a*cuneus_sa + a*entorhinal_sa + a*fusiform_sa + a*inferiorparietal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*lingual_sa + a*parahippocampal_sa + a*paracentral_sa + a*parsopercularis_sa + a*parsorbitalis_sa + a*pericalcarine_sa + a*postcentral_sa + a*precentral_sa + a*precuneus_sa + a*rostralmiddlefrontal_sa + a*transversetemporal_sa'

fit_sem_SAlat_15_reg_constrained <- sem(Model_SAlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SAlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
anova(fit_sem_SAlat_15_reg_free,fit_sem_SAlat_15_reg_constrained)

#Data Visualization for the regularized surface area model
BIC(fit_sem_SAlat_15_reg_free) 
data_subsample_SAlat_15_rl_reg_GGSEG <- data.frame(label=c("lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lingual","rh_lingual","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_transversetemporal","rh_transversetemporal"),
                                                   SA_p=rep(NA, times=36),
                                                   SA_std=rep(NA, times=36))
a=1
for (i in seq(from = 1, to = 36, by = 2)) {
  data_subsample_SAlat_15_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_SAlat_15_reg_free)[a+5,7]
  data_subsample_SAlat_15_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_SAlat_15_reg_free)[a+5,7]
  data_subsample_SAlat_15_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_SAlat_15_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_SAlat_15_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_SAlat_15_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_SAlat_15_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_SAlat_15_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#Model for surface area with all the regions of interest (34 ROIs)
Model_SAlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ bankssts_sa + caudalanteriorcingulate_sa + caudalmiddlefrontal_sa + cuneus_sa + entorhinal_sa + fusiform_sa + inferiorparietal_sa + inferiortemporal_sa + isthmuscingulate_sa + lateraloccipital_sa + lateralorbitofrontal_sa + lingual_sa + medialorbitofrontal_sa + middletemporal_sa + parahippocampal_sa + paracentral_sa + parsopercularis_sa + parsorbitalis_sa + parstriangularis_sa + pericalcarine_sa + postcentral_sa + posteriorcingulate_sa + precentral_sa + precuneus_sa + rostralanteriorcingulate_sa + rostralmiddlefrontal_sa + superiorfrontal_sa + superiorparietal_sa + superiortemporal_sa + supramarginal_sa + frontalpole_sa + temporalpole_sa + transversetemporal_sa + insula_sa'

fit_sem_SAlat_15_all_free <- sem(Model_SAlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SAlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the surface area model with all the regions
Model_SAlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*bankssts_sa + a*caudalanteriorcingulate_sa + a*caudalmiddlefrontal_sa + a*cuneus_sa + a*entorhinal_sa + a*fusiform_sa + a*inferiorparietal_sa + a*inferiortemporal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*lateralorbitofrontal_sa + a*lingual_sa + a*medialorbitofrontal_sa + a*middletemporal_sa + a*parahippocampal_sa + a*paracentral_sa + a*parsopercularis_sa + a*parsorbitalis_sa + a*parstriangularis_sa + a*pericalcarine_sa + a*postcentral_sa + a*posteriorcingulate_sa + a*precentral_sa + a*precuneus_sa + a*rostralanteriorcingulate_sa + a*rostralmiddlefrontal_sa + a*superiorfrontal_sa + a*superiorparietal_sa + a*superiortemporal_sa + a*supramarginal_sa + a*frontalpole_sa + a*temporalpole_sa + a*transversetemporal_sa + a*insula_sa'

fit_sem_SAlat_15_all_constrained <- sem(Model_SAlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_SAlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
anova(fit_sem_SAlat_15_all_free,fit_sem_SAlat_15_all_constrained)

#Data Visualization for the surface area model with all the regions
BIC(fit_sem_SAlat_15_all_free)
data_subsample_SAlat_15_rl_all_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                   SA_p=rep(NA, times=68),
                                                   SA_std=rep(NA, times=68))
a=1
for (i in seq(from = 1, to = 68, by = 2)) {
  data_subsample_SAlat_15_rl_all_GGSEG[i,2]<-parameterEstimates(fit_sem_SAlat_15_all_free)[a+5,7]
  data_subsample_SAlat_15_rl_all_GGSEG[i+1,2]<-parameterEstimates(fit_sem_SAlat_15_all_free)[a+5,7]
  data_subsample_SAlat_15_rl_all_GGSEG[i,3]<-lavInspect(fit_sem_SAlat_15_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_SAlat_15_rl_all_GGSEG[i+1,3]<-lavInspect(fit_sem_SAlat_15_all_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_SAlat_15_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_SAlat_15_rl_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

#------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Model per metric - Grey Matter Volume
#------------------------------------------------------------------------------------------------------------------------------
  
#Lasso regularization of the grey matter volume measures from all 34 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_gmv <- data_total_subsample15[,c(78:111,172)]
df_gmv <- makeX(df_gmv, na.impute = TRUE)
x_gmv <- df_gmv[,1:34]
y_gmv <- df_gmv[,35]

cv_lasso_GMVlat <- cv.glmnet(x_gmv,y_gmv, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_GMVlat)
new_coef_gmv <- coef(cv_lasso_GMVlat, s = "lambda.min")
new_coef_gmv <- as.data.frame(summary(new_coef_gmv))
new_coef_gmv <- new_coef_gmv[,c(1,3)]
names(new_coef_gmv) <- c("index","coef")

ROIs_coef_gmv <- data.frame(label=c("bankssts_gmv","caudalanteriorcingulate_gmv","caudalmiddlefrontal_gmv","cuneus_gmv","entorhinal_gmv","fusiform_gmv","inferiorparietal_gmv","inferiortemporal_gmv","isthmuscingulate_gmv","lateraloccipital_gmv","lateralorbitofrontal_gmv","lingual_gmv","medialorbitofrontal_gmv","middletemporal_gmv","parahippocampal_gmv","paracentral_gmv","parsopercularis_gmv","parsorbitalis_gmv","parstriangularis_gmv","pericalcarine_gmv","postcentral_gmv","posteriorcingulate_gmv","precentral_gmv","precuneus_gmv","rostralanteriorcingulate_gmv","rostralmiddlefrontal_gmv","superiorfrontal_gmv","superiorparietal_gmv","superiortemporal_gmv","supramarginal_gmv","frontalpole_gmv","temporalpole_gmv","transversetemporal_gmv","insula_gmv"),
                            coef=rep("NA",times=34))

ROIs_coef_gmv[new_coef_gmv[1,1],2] <- new_coef_gmv[1,2]

for (i in 1:nrow(new_coef_gmv)) {
  ROIs_coef_gmv[new_coef_gmv[i,1],2] <- new_coef_gmv[i,2]
}
ROI_gmv<-c()
for (i in 1:34) {
  if (ROIs_coef_gmv[i,2] != "NA") ROI_gmv[i]<- ROIs_coef_gmv[i,1]
}
ROI_gmv<-na.omit(ROI_gmv)

ROI_gmv_plus<-c(ROI_gmv[1])   #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_gmv)) {
  ROI_gmv_plus<-paste(ROI_gmv_plus,"+",ROI_gmv[i])
}
cat(ROI_gmv_plus)
  
#Regularized grey matter volume model 
Model_GMVlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ entorhinal_gmv + fusiform_gmv + inferiorparietal_gmv + lateraloccipital_gmv + parahippocampal_gmv + paracentral_gmv + pericalcarine_gmv + precentral_gmv + precuneus_gmv + supramarginal_gmv + insula_gmv'

fit_sem_GMVlat_15_reg_free <- sem(Model_GMVlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMVlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized grey matter volume model
Model_GMVlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ a*entorhinal_gmv + a*fusiform_gmv + a*inferiorparietal_gmv + a*lateraloccipital_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*pericalcarine_gmv + a*precentral_gmv + a*precuneus_gmv + a*supramarginal_gmv + a*insula_gmv'

fit_sem_GMVlat_15_reg_constrained <- sem(Model_GMVlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMVlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_GMVlat_15_reg_free,fit_sem_GMVlat_15_reg_constrained)

#Data Visualization for the regularized grey matter volume model
BIC(fit_sem_GMVlat_15_reg_free)                      
data_subsample_GMVlat_15_rl_reg_GGSEG <- data.frame(label=c("lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_lateraloccipital","rh_lateraloccipital","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_pericalcarine","rh_pericalcarine","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_supramarginal","rh_supramarginal","lh_insula","rh_insula"),
                                                    GMV_p=rep(NA, times=22),
                                                    GMV_std=rep(NA, times=22))
a=1
for (i in seq(from = 1, to = 22, by = 2)) {
  data_subsample_GMVlat_15_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_GMVlat_15_reg_free)[a+5,7]
  data_subsample_GMVlat_15_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMVlat_15_reg_free)[a+5,7]
  data_subsample_GMVlat_15_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_GMVlat_15_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMVlat_15_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_GMVlat_15_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_GMVlat_15_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMVlat_15_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

  
#Model for grey matter volume with all the regions of interest (34 ROIs)
Model_GMVlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ bankssts_gmv + caudalanteriorcingulate_gmv + caudalmiddlefrontal_gmv + cuneus_gmv + entorhinal_gmv + fusiform_gmv + inferiorparietal_gmv + inferiortemporal_gmv + isthmuscingulate_gmv + lateraloccipital_gmv + lateralorbitofrontal_gmv + lingual_gmv + medialorbitofrontal_gmv + middletemporal_gmv + parahippocampal_gmv + paracentral_gmv + parsopercularis_gmv + parsorbitalis_gmv + parstriangularis_gmv + pericalcarine_gmv + postcentral_gmv + posteriorcingulate_gmv + precentral_gmv + precuneus_gmv + rostralanteriorcingulate_gmv + rostralmiddlefrontal_gmv + superiorfrontal_gmv + superiorparietal_gmv + superiortemporal_gmv + supramarginal_gmv + frontalpole_gmv + temporalpole_gmv + transversetemporal_gmv + insula_gmv'

fit_sem_GMVlat_15_all_free <- sem(Model_GMVlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMVlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the grey matter volume model with all the regions
Model_GMVlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ a*bankssts_gmv + a*caudalanteriorcingulate_gmv + a*caudalmiddlefrontal_gmv + a*cuneus_gmv + a*entorhinal_gmv + a*fusiform_gmv + a*inferiorparietal_gmv + a*inferiortemporal_gmv + a*isthmuscingulate_gmv + a*lateraloccipital_gmv + a*lateralorbitofrontal_gmv + a*lingual_gmv + a*medialorbitofrontal_gmv + a*middletemporal_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*parsopercularis_gmv + a*parsorbitalis_gmv + a*parstriangularis_gmv + a*pericalcarine_gmv + a*postcentral_gmv + a*posteriorcingulate_gmv + a*precentral_gmv + a*precuneus_gmv + a*rostralanteriorcingulate_gmv + a*rostralmiddlefrontal_gmv + a*superiorfrontal_gmv + a*superiorparietal_gmv + a*superiortemporal_gmv + a*supramarginal_gmv + a*frontalpole_gmv + a*temporalpole_gmv + a*transversetemporal_gmv + a*insula_gmv'

fit_sem_GMVlat_15_all_constrained <- sem(Model_GMVlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMVlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
anova(fit_sem_GMVlat_15_all_free,fit_sem_GMVlat_15_all_constrained)

#Data Visualization for the grey matter volume model with all the regions
BIC(fit_sem_GMVlat_15_all_free)
data_subsample_GMVlat_15_rl_all_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                    GMV_p=rep(NA, times=68),
                                                    GMV_std=rep(NA, times=68))
a=1
for (i in seq(from = 1, to = 68, by = 2)) {
  data_subsample_GMVlat_15_rl_all_GGSEG[i,2]<-parameterEstimates(fit_sem_GMVlat_15_all_free)[a+5,7]
  data_subsample_GMVlat_15_rl_all_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMVlat_15_all_free)[a+5,7]
  data_subsample_GMVlat_15_rl_all_GGSEG[i,3]<-lavInspect(fit_sem_GMVlat_15_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMVlat_15_rl_all_GGSEG[i+1,3]<-lavInspect(fit_sem_GMVlat_15_all_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_GMVlat_15_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMVlat_15_rl_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Models including the three grey matter metrics
#------------------------------------------------------------------------------------------------------------------------------
  
#Lasso regularization of the grey matter measures from all 102 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_greym <- data_total_subsample15[,c(10:111,171)]
df_greym <- makeX(df_greym, na.impute = TRUE)
x_greym <- df_greym[,1:102]
y_greym <- df_greym[,103]

cv_lasso_greymlat <- cv.glmnet(x_greym,y_greym, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_greymlat)
new_coef_greym <- coef(cv_lasso_greymlat, s = "lambda.min")
new_coef_greym <- as.data.frame(summary(new_coef_greym))
new_coef_greym <- new_coef_greym[,c(1,3)]
names(new_coef_greym) <- c("index","coef")

ROIs_coef_greym <- data.frame(label=c("bankssts_ct","caudalanteriorcingulate_ct","caudalmiddlefrontal_ct","cuneus_ct","entorhinal_ct","fusiform_ct","inferiorparietal_ct","inferiortemporal_ct","isthmuscingulate_ct","lateraloccipital_ct","lateralorbitofrontal_ct","lingual_ct","medialorbitofrontal_ct","middletemporal_ct","parahippocampal_ct","paracentral_ct","parsopercularis_ct","parsorbitalis_ct","parstriangularis_ct","pericalcarine_ct","postcentral_ct","posteriorcingulate_ct","precentral_ct","precuneus_ct","rostralanteriorcingulate_ct","rostralmiddlefrontal_ct","superiorfrontal_ct","superiorparietal_ct","superiortemporal_ct","supramarginal_ct","frontalpole_ct","temporalpole_ct","transversetemporal_ct","insula_ct","bankssts_sa","caudalanteriorcingulate_sa","caudalmiddlefrontal_sa","cuneus_sa","entorhinal_sa","fusiform_sa","inferiorparietal_sa","inferiortemporal_sa","isthmuscingulate_sa","lateraloccipital_sa","lateralorbitofrontal_sa","lingual_sa","medialorbitofrontal_sa","middletemporal_sa","parahippocampal_sa","paracentral_sa","parsopercularis_sa","parsorbitalis_sa","parstriangularis_sa","pericalcarine_sa","postcentral_sa","posteriorcingulate_sa","precentral_sa","precuneus_sa","rostralanteriorcingulate_sa","rostralmiddlefrontal_sa","superiorfrontal_sa","superiorparietal_sa","superiortemporal_sa","supramarginal_sa","frontalpole_sa","temporalpole_sa","transversetemporal_sa","insula_sa","bankssts_gmv","caudalanteriorcingulate_gmv","caudalmiddlefrontal_gmv","cuneus_gmv","entorhinal_gmv","fusiform_gmv","inferiorparietal_gmv","inferiortemporal_gmv","isthmuscingulate_gmv","lateraloccipital_gmv","lateralorbitofrontal_gmv","lingual_gmv","medialorbitofrontal_gmv","middletemporal_gmv","parahippocampal_gmv","paracentral_gmv","parsopercularis_gmv","parsorbitalis_gmv","parstriangularis_gmv","pericalcarine_gmv","postcentral_gmv","posteriorcingulate_gmv","precentral_gmv","precuneus_gmv","rostralanteriorcingulate_gmv","rostralmiddlefrontal_gmv","superiorfrontal_gmv","superiorparietal_gmv","superiortemporal_gmv","supramarginal_gmv","frontalpole_gmv","temporalpole_gmv","transversetemporal_gmv","insula_gmv"),
                              coef=rep("NA",times=102))

ROIs_coef_greym[new_coef_greym[1,1],2] <- new_coef_greym[1,2]

for (i in 1:nrow(new_coef_greym)) {
  ROIs_coef_greym[new_coef_greym[i,1],2] <- new_coef_greym[i,2]
}
ROI_greym<-c()
for (i in 1:102) {
  if (ROIs_coef_greym[i,2] != "NA") ROI_greym[i]<- ROIs_coef_greym[i,1]
}
ROI_greym<-na.omit(ROI_greym)

ROI_greym_plus<-c(ROI_greym[1])    #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_greym)) {
  ROI_greym_plus<-paste(ROI_greym_plus,"+",ROI_greym[i])
}
cat(ROI_greym_plus)

#Regularized grey matter metrics model 
Model_GMlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + entorhinal_ct + fusiform_ct + inferiorparietal_ct + lateraloccipital_ct + lateralorbitofrontal_ct + lingual_ct + middletemporal_ct + parahippocampal_ct + parsopercularis_ct + pericalcarine_ct + superiorfrontal_ct + superiorparietal_ct + frontalpole_ct + temporalpole_ct + insula_ct + caudalmiddlefrontal_sa + entorhinal_sa + inferiorparietal_sa + isthmuscingulate_sa + lateraloccipital_sa + medialorbitofrontal_sa + middletemporal_sa + paracentral_sa + postcentral_sa + precentral_sa + precuneus_sa + superiorfrontal_sa + superiorparietal_sa + superiortemporal_sa + transversetemporal_sa + insula_sa + bankssts_gmv + cuneus_gmv + parahippocampal_gmv + paracentral_gmv + parstriangularis_gmv + precuneus_gmv + rostralmiddlefrontal_gmv + superiorfrontal_gmv + temporalpole_gmv'

fit_sem_GMlat_15_reg_free <- sem(Model_GMlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized grey matter metrics model
Model_GMlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*entorhinal_ct + a*fusiform_ct + a*inferiorparietal_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*lingual_ct + a*middletemporal_ct + a*parahippocampal_ct + a*parsopercularis_ct + a*pericalcarine_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*frontalpole_ct + a*temporalpole_ct + a*insula_ct + a*caudalmiddlefrontal_sa + a*entorhinal_sa + a*inferiorparietal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*medialorbitofrontal_sa + a*middletemporal_sa + a*paracentral_sa + a*postcentral_sa + a*precentral_sa + a*precuneus_sa + a*superiorfrontal_sa + a*superiorparietal_sa + a*superiortemporal_sa + a*transversetemporal_sa + a*insula_sa + a*bankssts_gmv + a*cuneus_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*parstriangularis_gmv + a*precuneus_gmv + a*rostralmiddlefrontal_gmv + a*superiorfrontal_gmv + a*temporalpole_gmv'

fit_sem_GMlat_15_reg_constrained <- sem(Model_GMlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
anova(fit_sem_GMlat_15_reg_free,fit_sem_GMlat_15_reg_constrained)

#Data Visualization for the regularized grey matter metrics model
BIC(fit_sem_GMlat_15_reg_free)           
data_subsample_GMlat_15_regCT_GGSEG <- data.frame(label=c("lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_parsopercularis","rh_parsopercularis","lh_pericalcarine","rh_pericalcarine","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_insula","rh_insula"),
                                                  CT_p=rep(NA, times=34),
                                                  CT_std=rep(NA, times=34))
a=1
for (i in seq(from = 1, to = 34, by = 2)) {
  data_subsample_GMlat_15_regCT_GGSEG[i,2]<-parameterEstimates(fit_sem_GMlat_15_reg_free)[a+5,7]
  data_subsample_GMlat_15_regCT_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMlat_15_reg_free)[a+5,7]
  data_subsample_GMlat_15_regCT_GGSEG[i,3]<-lavInspect(fit_sem_GMlat_15_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMlat_15_regCT_GGSEG[i+1,3]<-lavInspect(fit_sem_GMlat_15_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_GMlat_15_regSA_GGSEG <- data.frame(label=c("lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_entorhinal","rh_entorhinal","lh_inferiorparietal","rh_inferiorparietal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_paracentral","rh_paracentral","lh_postcentral","rh_postcentral","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                  SA_p=rep(NA, times=32),
                                                  SA_std=rep(NA, times=32))
a=1
for (i in seq(from = 1, to = 32, by = 2)) {
  data_subsample_GMlat_15_regSA_GGSEG[i,2]<-parameterEstimates(fit_sem_GMlat_15_reg_free)[a+19,7]
  data_subsample_GMlat_15_regSA_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMlat_15_reg_free)[a+19,7]
  data_subsample_GMlat_15_regSA_GGSEG[i,3]<-lavInspect(fit_sem_GMlat_15_reg_free,what = "std.all")$beta[1,a+15]
  data_subsample_GMlat_15_regSA_GGSEG[i+1,3]<-lavInspect(fit_sem_GMlat_15_reg_free,what = "std.all")$beta[1,a+15]
  a=a+1
}

                                                                                                                                                                                                                                                                                                                                                
data_subsample_GMlat_15_regGMV_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_cuneus","rh_cuneus","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parstriangularis","rh_parstriangularis","lh_precuneus","rh_precuneus","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_temporalpole","rh_temporalpole"),
                                                   GMV_p=rep(NA, times=18),
                                                   GMV_std=rep(NA, times=18))
a=1
for (i in seq(from = 1, to = 18, by = 2)) {
  data_subsample_GMlat_15_regGMV_GGSEG[i,2]<-parameterEstimates(fit_sem_GMlat_15_reg_free)[a+30,7]
  data_subsample_GMlat_15_regGMV_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMlat_15_reg_free)[a+30,7]
  data_subsample_GMlat_15_regGMV_GGSEG[i,3]<-lavInspect(fit_sem_GMlat_15_reg_free,what = "std.all")$beta[1,a+26]
  data_subsample_GMlat_15_regGMV_GGSEG[i+1,3]<-lavInspect(fit_sem_GMlat_15_reg_free,what = "std.all")$beta[1,a+26]
  a=a+1
}

data_subsample_GMlat_15_regCT_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_regCT_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
data_subsample_GMlat_15_regSA_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_regSA_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
data_subsample_GMlat_15_regGMV_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_regGMV_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

#Compare the information provide by the regularized regions in the three grey matter metrics 
data_subsample_GMlat_15_reg_GGSEG<-merge(data_subsample_GMlat_15_regCT_GGSEG[,c(1,3)],data_subsample_GMlat_15_regSA_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_GMlat_15_reg_GGSEG<-merge(data_subsample_GMlat_15_reg_GGSEG,data_subsample_GMlat_15_regGMV_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_GMlat_15_reg_GGSEG, "STD_GMlat_15_reg.csv", row.names=FALSE)

data_subsample_GMlat_15_reg_GGSEG_long <- data_subsample_GMlat_15_reg_GGSEG %>%
  rename(c("CT_std" = "CT")) %>%
  rename(c("SA_std" = "SA")) %>%
  rename(c("GMV_std" = "GMV")) %>%
  gather(Metric, STD, CT, SA, GMV) 
data_subsample_GMlat_15_reg_GGSEG_long$Metric<-factor(data_subsample_GMlat_15_reg_GGSEG_long$Metric,levels=c("CT","SA","GMV"))
metric_names <- c(CT = "Cortical Thickness",SA = "Surface Area",GMV = "Grey Matter Volume")

data_subsample_GMlat_15_reg_GGSEG_long %>%
  group_by(Metric) %>%
  ggseg(hemisphere="left",colour="black",size=.6,
        mapping=aes(fill=STD), show.legend = TRUE) +
  theme_classic(base_size = 20)+
  theme(axis.title.x = element_text(size = 15))+
  scale_fill_gradient(low="goldenrod",high="firebrick")+
  facet_wrap(~Metric, nrow=3,labeller = as_labeller(metric_names))+
  labs(fill="Standardized\nestimation\nparameters")
ggsave("plot_GMlat_15_reg.png")
  
#Model for grey matter metrics with all the regions of interest (102 ROIs)
Model_GMlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ bankssts_ct + caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + cuneus_ct + entorhinal_ct + fusiform_ct + inferiorparietal_ct + inferiortemporal_ct + isthmuscingulate_ct + lateraloccipital_ct + lateralorbitofrontal_ct + lingual_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsopercularis_ct + parsorbitalis_ct + parstriangularis_ct + pericalcarine_ct + postcentral_ct + posteriorcingulate_ct + precentral_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiorparietal_ct + superiortemporal_ct + supramarginal_ct + frontalpole_ct + temporalpole_ct + transversetemporal_ct + insula_ct + bankssts_sa + caudalanteriorcingulate_sa + caudalmiddlefrontal_sa + cuneus_sa + entorhinal_sa + fusiform_sa + inferiorparietal_sa + inferiortemporal_sa + isthmuscingulate_sa + lateraloccipital_sa + lateralorbitofrontal_sa + lingual_sa + medialorbitofrontal_sa + middletemporal_sa + parahippocampal_sa + paracentral_sa + parsopercularis_sa + parsorbitalis_sa + parstriangularis_sa + pericalcarine_sa + postcentral_sa + posteriorcingulate_sa + precentral_sa + precuneus_sa + rostralanteriorcingulate_sa + rostralmiddlefrontal_sa + superiorfrontal_sa + superiorparietal_sa + superiortemporal_sa + supramarginal_sa + frontalpole_sa + temporalpole_sa + transversetemporal_sa + insula_sa + bankssts_gmv + caudalanteriorcingulate_gmv + caudalmiddlefrontal_gmv + cuneus_gmv + entorhinal_gmv + fusiform_gmv + inferiorparietal_gmv + inferiortemporal_gmv + isthmuscingulate_gmv + lateraloccipital_gmv + lateralorbitofrontal_gmv + lingual_gmv + medialorbitofrontal_gmv + middletemporal_gmv + parahippocampal_gmv + paracentral_gmv + parsopercularis_gmv + parsorbitalis_gmv + parstriangularis_gmv + pericalcarine_gmv + postcentral_gmv + posteriorcingulate_gmv + precentral_gmv + precuneus_gmv + rostralanteriorcingulate_gmv + rostralmiddlefrontal_gmv + superiorfrontal_gmv + superiorparietal_gmv + superiortemporal_gmv + supramarginal_gmv + frontalpole_gmv + temporalpole_gmv + transversetemporal_gmv + insula_gmv'

fit_sem_GMlat_15_all_free <- sem(Model_GMlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the grey matter metrics model with all the regions
Model_GMlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*bankssts_ct + a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*cuneus_ct + a*entorhinal_ct + a*fusiform_ct + a*inferiorparietal_ct + a*inferiortemporal_ct + a*isthmuscingulate_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*lingual_ct + a*medialorbitofrontal_ct + a*middletemporal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsopercularis_ct + a*parsorbitalis_ct + a*parstriangularis_ct + a*pericalcarine_ct + a*postcentral_ct + a*posteriorcingulate_ct + a*precentral_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*superiortemporal_ct + a*supramarginal_ct + a*frontalpole_ct + a*temporalpole_ct + a*transversetemporal_ct + a*insula_ct + a*bankssts_sa + a*caudalanteriorcingulate_sa + a*caudalmiddlefrontal_sa + a*cuneus_sa + a*entorhinal_sa + a*fusiform_sa + a*inferiorparietal_sa + a*inferiortemporal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*lateralorbitofrontal_sa + a*lingual_sa + a*medialorbitofrontal_sa + a*middletemporal_sa + a*parahippocampal_sa + a*paracentral_sa + a*parsopercularis_sa + a*parsorbitalis_sa + a*parstriangularis_sa + a*pericalcarine_sa + a*postcentral_sa + a*posteriorcingulate_sa + a*precentral_sa + a*precuneus_sa + a*rostralanteriorcingulate_sa + a*rostralmiddlefrontal_sa + a*superiorfrontal_sa + a*superiorparietal_sa + a*superiortemporal_sa + a*supramarginal_sa + a*frontalpole_sa + a*temporalpole_sa + a*transversetemporal_sa + a*insula_sa + a*bankssts_gmv + a*caudalanteriorcingulate_gmv + a*caudalmiddlefrontal_gmv + a*cuneus_gmv + a*entorhinal_gmv + a*fusiform_gmv + a*inferiorparietal_gmv + a*inferiortemporal_gmv + a*isthmuscingulate_gmv + a*lateraloccipital_gmv + a*lateralorbitofrontal_gmv + a*lingual_gmv + a*medialorbitofrontal_gmv + a*middletemporal_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*parsopercularis_gmv + a*parsorbitalis_gmv + a*parstriangularis_gmv + a*pericalcarine_gmv + a*postcentral_gmv + a*posteriorcingulate_gmv + a*precentral_gmv + a*precuneus_gmv + a*rostralanteriorcingulate_gmv + a*rostralmiddlefrontal_gmv + a*superiorfrontal_gmv + a*superiorparietal_gmv + a*superiortemporal_gmv + a*supramarginal_gmv + a*frontalpole_gmv + a*temporalpole_gmv + a*transversetemporal_gmv + a*insula_gmv'

fit_sem_GMlat_15_all_constrained <- sem(Model_GMlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_GMlat_15_all_free,fit_sem_GMlat_15_all_constrained)

#Data Visualization for the grey matter volume model with all the regions
BIC(fit_sem_GMlat_15_all_free)
data_subsample_GMlat_15_rl_all_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                   CT_p=rep(NA, times=68),
                                                   CT_std=rep(NA, times=68),
                                                   SA_p=rep(NA, times=68),
                                                   SA_std=rep(NA, times=68),
                                                   GMV_p=rep(NA, times=68),
                                                   GMV_std=rep(NA, times=68))

a=1
for (i in seq(from = 1, to = 68, by = 2)) {
  data_subsample_GMlat_15_rl_all_GGSEG[i,2]<-parameterEstimates(fit_sem_GMlat_15_all_free)[a+5,7]
  data_subsample_GMlat_15_rl_all_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMlat_15_all_free)[a+5,7]
  data_subsample_GMlat_15_rl_all_GGSEG[i,4]<-parameterEstimates(fit_sem_GMlat_15_all_free)[a+39,7]
  data_subsample_GMlat_15_rl_all_GGSEG[i+1,4]<-parameterEstimates(fit_sem_GMlat_15_all_free)[a+39,7]
  data_subsample_GMlat_15_rl_all_GGSEG[i,6]<-parameterEstimates(fit_sem_GMlat_15_all_free)[a+73,7]
  data_subsample_GMlat_15_rl_all_GGSEG[i+1,6]<-parameterEstimates(fit_sem_GMlat_15_all_free)[a+73,7]
  data_subsample_GMlat_15_rl_all_GGSEG[i,3]<-lavInspect(fit_sem_GMlat_15_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMlat_15_rl_all_GGSEG[i+1,3]<-lavInspect(fit_sem_GMlat_15_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMlat_15_rl_all_GGSEG[i,5]<-lavInspect(fit_sem_GMlat_15_all_free,what = "std.all")$beta[1,a+35]
  data_subsample_GMlat_15_rl_all_GGSEG[i+1,5]<-lavInspect(fit_sem_GMlat_15_all_free,what = "std.all")$beta[1,a+35]
  data_subsample_GMlat_15_rl_all_GGSEG[i,7]<-lavInspect(fit_sem_GMlat_15_all_free,what = "std.all")$beta[1,a+69]
  data_subsample_GMlat_15_rl_all_GGSEG[i+1,7]<-lavInspect(fit_sem_GMlat_15_all_free,what = "std.all")$beta[1,a+69]
  a=a+1
}

data_subsample_GMlat_15_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_rl_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
data_subsample_GMlat_15_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_rl_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
data_subsample_GMlat_15_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_15_rl_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")




#------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Individual models white matter
#------------------------------------------------------------------------------------------------------------------------------
  
#Create a table to save the regression coefficients or the fit of the ROI model for each ROIs
data_subsample_WMlat_15_GGSEG <- data.frame(label=c("fornix","cingulatecingulum","parahippocampalcingulum","corticospinalpyramidal","anteriorthalamicradiations","uncinate","inferiorlongitudinalfasiculus","inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","corpuscallosum","superiorlongitudinalfasiculus","temporalsuperiorlongitudinalfasiculus","parietalsuperiorlongitudinalfasiculus","superiorcorticostriate","superiorcorticostriatefrontalcortex","superiorcorticostriateparietalcortex","striatalinferiorfrontalcortex","inferiorfrontalsuperiorfrontalcortex","fornix_exfimbria"),
                                              FA_fit=rep(NA, times=20),
                                              FA_p=rep(NA, times=20),
                                              FA_std=rep(NA, times=20),
                                              MD_fit=rep(NA, times=20),
                                              MD_p=rep(NA, times=20),
                                              MD_std=rep(NA, times=20),
                                              WMV_fit=rep(NA, times=20),
                                              WMV_p=rep(NA, times=20),
                                              WMV_std=rep(NA, times=20))

#Fractional Anisotropy (20 regions)
##fornix
Model_FA_fornix <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                    Cognitive_factor ~ fornix_fa'
fit_sem_FA_fornix <- sem(Model_FA_fornix, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_fornix, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[1,2]<-BIC(fit_sem_FA_fornix)
data_subsample_WMlat_15_GGSEG[1,3]<-parameterEstimates(fit_sem_FA_fornix)[6,7]
data_subsample_WMlat_15_GGSEG[1,4]<-lavInspect(fit_sem_FA_fornix,what = "std.all")$beta[1,2]

##cingulatecingulum
Model_FA_cingulatecingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                               Cognitive_factor ~ cingulatecingulum_fa'
fit_sem_FA_cingulatecingulum <- sem(Model_FA_cingulatecingulum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_cingulatecingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[2,2]<-BIC(fit_sem_FA_cingulatecingulum)
data_subsample_WMlat_15_GGSEG[2,3]<-parameterEstimates(fit_sem_FA_cingulatecingulum)[6,7]
data_subsample_WMlat_15_GGSEG[2,4]<-lavInspect(fit_sem_FA_cingulatecingulum,what = "std.all")$beta[1,2]

##parahippocampalcingulum
Model_FA_parahippocampalcingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ parahippocampalcingulum_fa'
fit_sem_FA_parahippocampalcingulum <- sem(Model_FA_parahippocampalcingulum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_parahippocampalcingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[3,2]<-BIC(fit_sem_FA_parahippocampalcingulum)
data_subsample_WMlat_15_GGSEG[3,3]<-parameterEstimates(fit_sem_FA_parahippocampalcingulum)[6,7]
data_subsample_WMlat_15_GGSEG[3,4]<-lavInspect(fit_sem_FA_parahippocampalcingulum,what = "std.all")$beta[1,2]

##corticospinalpyramidal
Model_FA_corticospinalpyramidal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ corticospinalpyramidal_fa'
fit_sem_FA_corticospinalpyramidal <- sem(Model_FA_corticospinalpyramidal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_corticospinalpyramidal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[4,2]<-BIC(fit_sem_FA_corticospinalpyramidal)
data_subsample_WMlat_15_GGSEG[4,3]<-parameterEstimates(fit_sem_FA_corticospinalpyramidal)[6,7]
data_subsample_WMlat_15_GGSEG[4,4]<-lavInspect(fit_sem_FA_corticospinalpyramidal,what = "std.all")$beta[1,2]

##anteriorthalamicradiations
Model_FA_anteriorthalamicradiations <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                        Cognitive_factor ~ anteriorthalamicradiations_fa'
fit_sem_FA_anteriorthalamicradiations <- sem(Model_FA_anteriorthalamicradiations, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_anteriorthalamicradiations, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[5,2]<-BIC(fit_sem_FA_anteriorthalamicradiations)
data_subsample_WMlat_15_GGSEG[5,3]<-parameterEstimates(fit_sem_FA_anteriorthalamicradiations)[6,7]
data_subsample_WMlat_15_GGSEG[5,4]<-lavInspect(fit_sem_FA_anteriorthalamicradiations,what = "std.all")$beta[1,2]

##uncinate
Model_FA_uncinate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                      Cognitive_factor ~ uncinate_fa'
fit_sem_FA_uncinate <- sem(Model_FA_uncinate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_uncinate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[6,2]<-BIC(fit_sem_FA_uncinate)
data_subsample_WMlat_15_GGSEG[6,3]<-parameterEstimates(fit_sem_FA_uncinate)[6,7]
data_subsample_WMlat_15_GGSEG[6,4]<-lavInspect(fit_sem_FA_uncinate,what = "std.all")$beta[1,2]

##inferiorlongitudinalfasiculus
Model_FA_inferiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                           Cognitive_factor ~ inferiorlongitudinalfasiculus_fa'
fit_sem_FA_inferiorlongitudinalfasiculus <- sem(Model_FA_inferiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_inferiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[7,2]<-BIC(fit_sem_FA_inferiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[7,3]<-parameterEstimates(fit_sem_FA_inferiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[7,4]<-lavInspect(fit_sem_FA_inferiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##inferiorfrontooccipitalfasiculus
Model_FA_inferiorfrontooccipitalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                              Cognitive_factor ~ inferiorfrontooccipitalfasiculus_fa'
fit_sem_FA_inferiorfrontooccipitalfasiculus <- sem(Model_FA_inferiorfrontooccipitalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_inferiorfrontooccipitalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[8,2]<-BIC(fit_sem_FA_inferiorfrontooccipitalfasiculus)
data_subsample_WMlat_15_GGSEG[8,3]<-parameterEstimates(fit_sem_FA_inferiorfrontooccipitalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[8,4]<-lavInspect(fit_sem_FA_inferiorfrontooccipitalfasiculus,what = "std.all")$beta[1,2]

##forcepsmajor
Model_FA_forcepsmajor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ forcepsmajor_fa'
fit_sem_FA_forcepsmajor <- sem(Model_FA_forcepsmajor, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_forcepsmajor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[9,2]<-BIC(fit_sem_FA_forcepsmajor)
data_subsample_WMlat_15_GGSEG[9,3]<-parameterEstimates(fit_sem_FA_forcepsmajor)[6,7]
data_subsample_WMlat_15_GGSEG[9,4]<-lavInspect(fit_sem_FA_forcepsmajor,what = "std.all")$beta[1,2]

##forcepsminor
Model_FA_forcepsminor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ forcepsminor_fa'
fit_sem_FA_forcepsminor <- sem(Model_FA_forcepsminor, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_forcepsminor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[10,2]<-BIC(fit_sem_FA_forcepsminor)
data_subsample_WMlat_15_GGSEG[10,3]<-parameterEstimates(fit_sem_FA_forcepsminor)[6,7]
data_subsample_WMlat_15_GGSEG[10,4]<-lavInspect(fit_sem_FA_forcepsminor,what = "std.all")$beta[1,2]

##corpuscallosum
Model_FA_corpuscallosum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ corpuscallosum_fa'
fit_sem_FA_corpuscallosum <- sem(Model_FA_corpuscallosum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_corpuscallosum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[11,2]<-BIC(fit_sem_FA_corpuscallosum)
data_subsample_WMlat_15_GGSEG[11,3]<-parameterEstimates(fit_sem_FA_corpuscallosum)[6,7]
data_subsample_WMlat_15_GGSEG[11,4]<-lavInspect(fit_sem_FA_corpuscallosum,what = "std.all")$beta[1,2]

##superiorlongitudinalfasiculus
Model_FA_superiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                           Cognitive_factor ~ superiorlongitudinalfasiculus_fa'
fit_sem_FA_superiorlongitudinalfasiculus <- sem(Model_FA_superiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_superiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[12,2]<-BIC(fit_sem_FA_superiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[12,3]<-parameterEstimates(fit_sem_FA_superiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[12,4]<-lavInspect(fit_sem_FA_superiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##temporalsuperiorlongitudinalfasiculus
Model_FA_temporalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                   Cognitive_factor ~ temporalsuperiorlongitudinalfasiculus_fa'
fit_sem_FA_temporalsuperiorlongitudinalfasiculus <- sem(Model_FA_temporalsuperiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_temporalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[13,2]<-BIC(fit_sem_FA_temporalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[13,3]<-parameterEstimates(fit_sem_FA_temporalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[13,4]<-lavInspect(fit_sem_FA_temporalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##parietalsuperiorlongitudinalfasiculus
Model_FA_parietalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                   Cognitive_factor ~ parietalsuperiorlongitudinalfasiculus_fa'
fit_sem_FA_parietalsuperiorlongitudinalfasiculus <- sem(Model_FA_parietalsuperiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_parietalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[14,2]<-BIC(fit_sem_FA_parietalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[14,3]<-parameterEstimates(fit_sem_FA_parietalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[14,4]<-lavInspect(fit_sem_FA_parietalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##superiorcorticostriate
Model_FA_superiorcorticostriate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ superiorcorticostriate_fa'
fit_sem_FA_superiorcorticostriate <- sem(Model_FA_superiorcorticostriate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_superiorcorticostriate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[15,2]<-BIC(fit_sem_FA_superiorcorticostriate)
data_subsample_WMlat_15_GGSEG[15,3]<-parameterEstimates(fit_sem_FA_superiorcorticostriate)[6,7]
data_subsample_WMlat_15_GGSEG[15,4]<-lavInspect(fit_sem_FA_superiorcorticostriate,what = "std.all")$beta[1,2]

##superiorcorticostriatefrontalcortex
Model_FA_superiorcorticostriatefrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                 Cognitive_factor ~ superiorcorticostriatefrontalcortex_fa'
fit_sem_FA_superiorcorticostriatefrontalcortex <- sem(Model_FA_superiorcorticostriatefrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_superiorcorticostriatefrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[16,2]<-BIC(fit_sem_FA_superiorcorticostriatefrontalcortex)
data_subsample_WMlat_15_GGSEG[16,3]<-parameterEstimates(fit_sem_FA_superiorcorticostriatefrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[16,4]<-lavInspect(fit_sem_FA_superiorcorticostriatefrontalcortex,what = "std.all")$beta[1,2]

##superiorcorticostriateparietalcortex
Model_FA_superiorcorticostriateparietalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                  Cognitive_factor ~ superiorcorticostriateparietalcortex_fa'
fit_sem_FA_superiorcorticostriateparietalcortex <- sem(Model_FA_superiorcorticostriateparietalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_superiorcorticostriateparietalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[17,2]<-BIC(fit_sem_FA_superiorcorticostriateparietalcortex)
data_subsample_WMlat_15_GGSEG[17,3]<-parameterEstimates(fit_sem_FA_superiorcorticostriateparietalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[17,4]<-lavInspect(fit_sem_FA_superiorcorticostriateparietalcortex,what = "std.all")$beta[1,2]

##striatalinferiorfrontalcortex
Model_FA_striatalinferiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                           Cognitive_factor ~ striatalinferiorfrontalcortex_fa'
fit_sem_FA_striatalinferiorfrontalcortex <- sem(Model_FA_striatalinferiorfrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_striatalinferiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[18,2]<-BIC(fit_sem_FA_striatalinferiorfrontalcortex)
data_subsample_WMlat_15_GGSEG[18,3]<-parameterEstimates(fit_sem_FA_striatalinferiorfrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[18,4]<-lavInspect(fit_sem_FA_striatalinferiorfrontalcortex,what = "std.all")$beta[1,2]

##inferiorfrontalsuperiorfrontalcortex
Model_FA_inferiorfrontalsuperiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                  Cognitive_factor ~ inferiorfrontalsuperiorfrontalcortex_fa'
fit_sem_FA_inferiorfrontalsuperiorfrontalcortex <- sem(Model_FA_inferiorfrontalsuperiorfrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[19,2]<-BIC(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex)
data_subsample_WMlat_15_GGSEG[19,3]<-parameterEstimates(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[19,4]<-lavInspect(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex,what = "std.all")$beta[1,2]

##fornix_exfimbria
Model_FA_fornix_exfimbria <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ fornix_exfimbria_fa'
fit_sem_FA_fornix_exfimbria <- sem(Model_FA_fornix_exfimbria, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FA_fornix_exfimbria, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[20,2]<-BIC(fit_sem_FA_fornix_exfimbria)
data_subsample_WMlat_15_GGSEG[20,3]<-parameterEstimates(fit_sem_FA_fornix_exfimbria)[6,7]
data_subsample_WMlat_15_GGSEG[20,4]<-lavInspect(fit_sem_FA_fornix_exfimbria,what = "std.all")$beta[1,2]

#Mean Diffusivity (20 regions)
##fornix
Model_MD_fornix <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                    Cognitive_factor ~ fornix_md'
fit_sem_MD_fornix <- sem(Model_MD_fornix, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_fornix, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[1,5]<-BIC(fit_sem_MD_fornix)
data_subsample_WMlat_15_GGSEG[1,6]<-parameterEstimates(fit_sem_MD_fornix)[6,7]
data_subsample_WMlat_15_GGSEG[1,7]<-lavInspect(fit_sem_MD_fornix,what = "std.all")$beta[1,2]

##cingulatecingulum
Model_MD_cingulatecingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                               Cognitive_factor ~ cingulatecingulum_md'
fit_sem_MD_cingulatecingulum <- sem(Model_MD_cingulatecingulum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_cingulatecingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[2,5]<-BIC(fit_sem_MD_cingulatecingulum)
data_subsample_WMlat_15_GGSEG[2,6]<-parameterEstimates(fit_sem_MD_cingulatecingulum)[6,7]
data_subsample_WMlat_15_GGSEG[2,7]<-lavInspect(fit_sem_MD_cingulatecingulum,what = "std.all")$beta[1,2]

##parahippocampalcingulum
Model_MD_parahippocampalcingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ parahippocampalcingulum_md'
fit_sem_MD_parahippocampalcingulum <- sem(Model_MD_parahippocampalcingulum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_parahippocampalcingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[3,5]<-BIC(fit_sem_MD_parahippocampalcingulum)
data_subsample_WMlat_15_GGSEG[3,6]<-parameterEstimates(fit_sem_MD_parahippocampalcingulum)[6,7]
data_subsample_WMlat_15_GGSEG[3,7]<-lavInspect(fit_sem_MD_parahippocampalcingulum,what = "std.all")$beta[1,2]

##corticospinalpyramidal
Model_MD_corticospinalpyramidal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ corticospinalpyramidal_md'
fit_sem_MD_corticospinalpyramidal <- sem(Model_MD_corticospinalpyramidal, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_corticospinalpyramidal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[4,5]<-BIC(fit_sem_MD_corticospinalpyramidal)
data_subsample_WMlat_15_GGSEG[4,6]<-parameterEstimates(fit_sem_MD_corticospinalpyramidal)[6,7]
data_subsample_WMlat_15_GGSEG[4,7]<-lavInspect(fit_sem_MD_corticospinalpyramidal,what = "std.all")$beta[1,2]

##anteriorthalamicradiations
Model_MD_anteriorthalamicradiations <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                        Cognitive_factor ~ anteriorthalamicradiations_md'
fit_sem_MD_anteriorthalamicradiations <- sem(Model_MD_anteriorthalamicradiations, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_anteriorthalamicradiations, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[5,5]<-BIC(fit_sem_MD_anteriorthalamicradiations)
data_subsample_WMlat_15_GGSEG[5,6]<-parameterEstimates(fit_sem_MD_anteriorthalamicradiations)[6,7]
data_subsample_WMlat_15_GGSEG[5,7]<-lavInspect(fit_sem_MD_anteriorthalamicradiations,what = "std.all")$beta[1,2]

##uncinate
Model_MD_uncinate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                      Cognitive_factor ~ uncinate_md'
fit_sem_MD_uncinate <- sem(Model_MD_uncinate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_uncinate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[6,5]<-BIC(fit_sem_MD_uncinate)
data_subsample_WMlat_15_GGSEG[6,6]<-parameterEstimates(fit_sem_MD_uncinate)[6,7]
data_subsample_WMlat_15_GGSEG[6,7]<-lavInspect(fit_sem_MD_uncinate,what = "std.all")$beta[1,2]

##inferiorlongitudinalfasiculus
Model_MD_inferiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                           Cognitive_factor ~ inferiorlongitudinalfasiculus_md'
fit_sem_MD_inferiorlongitudinalfasiculus <- sem(Model_MD_inferiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_inferiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[7,5]<-BIC(fit_sem_MD_inferiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[7,6]<-parameterEstimates(fit_sem_MD_inferiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[7,7]<-lavInspect(fit_sem_MD_inferiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##inferiorfrontooccipitalfasiculus
Model_MD_inferiorfrontooccipitalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                              Cognitive_factor ~ inferiorfrontooccipitalfasiculus_md'
fit_sem_MD_inferiorfrontooccipitalfasiculus <- sem(Model_MD_inferiorfrontooccipitalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_inferiorfrontooccipitalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[8,5]<-BIC(fit_sem_MD_inferiorfrontooccipitalfasiculus)
data_subsample_WMlat_15_GGSEG[8,6]<-parameterEstimates(fit_sem_MD_inferiorfrontooccipitalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[8,7]<-lavInspect(fit_sem_MD_inferiorfrontooccipitalfasiculus,what = "std.all")$beta[1,2]

##forcepsmajor
Model_MD_forcepsmajor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ forcepsmajor_md'
fit_sem_MD_forcepsmajor <- sem(Model_MD_forcepsmajor, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_forcepsmajor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[9,5]<-BIC(fit_sem_MD_forcepsmajor)
data_subsample_WMlat_15_GGSEG[9,6]<-parameterEstimates(fit_sem_MD_forcepsmajor)[6,7]
data_subsample_WMlat_15_GGSEG[9,7]<-lavInspect(fit_sem_MD_forcepsmajor,what = "std.all")$beta[1,2]

##forcepsminor
Model_MD_forcepsminor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                           Cognitive_factor ~ forcepsminor_md'
fit_sem_MD_forcepsminor <- sem(Model_MD_forcepsminor, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_forcepsminor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[10,5]<-BIC(fit_sem_MD_forcepsminor)
data_subsample_WMlat_15_GGSEG[10,6]<-parameterEstimates(fit_sem_MD_forcepsminor)[6,7]
data_subsample_WMlat_15_GGSEG[10,7]<-lavInspect(fit_sem_MD_forcepsminor,what = "std.all")$beta[1,2]

##corpuscallosum
Model_MD_corpuscallosum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ corpuscallosum_md'
fit_sem_MD_corpuscallosum <- sem(Model_MD_corpuscallosum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_corpuscallosum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[11,5]<-BIC(fit_sem_MD_corpuscallosum)
data_subsample_WMlat_15_GGSEG[11,6]<-parameterEstimates(fit_sem_MD_corpuscallosum)[6,7]
data_subsample_WMlat_15_GGSEG[11,7]<-lavInspect(fit_sem_MD_corpuscallosum,what = "std.all")$beta[1,2]

##superiorlongitudinalfasiculus
Model_MD_superiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                           Cognitive_factor ~ superiorlongitudinalfasiculus_md'
fit_sem_MD_superiorlongitudinalfasiculus <- sem(Model_MD_superiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_superiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[12,5]<-BIC(fit_sem_MD_superiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[12,6]<-parameterEstimates(fit_sem_MD_superiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[12,7]<-lavInspect(fit_sem_MD_superiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##temporalsuperiorlongitudinalfasiculus
Model_MD_temporalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                   Cognitive_factor ~ temporalsuperiorlongitudinalfasiculus_md'
fit_sem_MD_temporalsuperiorlongitudinalfasiculus <- sem(Model_MD_temporalsuperiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_temporalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[13,5]<-BIC(fit_sem_MD_temporalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[13,6]<-parameterEstimates(fit_sem_MD_temporalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[13,7]<-lavInspect(fit_sem_MD_temporalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##parietalsuperiorlongitudinalfasiculus
Model_MD_parietalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                   Cognitive_factor ~ parietalsuperiorlongitudinalfasiculus_md'
fit_sem_MD_parietalsuperiorlongitudinalfasiculus <- sem(Model_MD_parietalsuperiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_parietalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[14,5]<-BIC(fit_sem_MD_parietalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[14,6]<-parameterEstimates(fit_sem_MD_parietalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[14,7]<-lavInspect(fit_sem_MD_parietalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##superiorcorticostriate
Model_MD_superiorcorticostriate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ superiorcorticostriate_md'
fit_sem_MD_superiorcorticostriate <- sem(Model_MD_superiorcorticostriate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_superiorcorticostriate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[15,5]<-BIC(fit_sem_MD_superiorcorticostriate)
data_subsample_WMlat_15_GGSEG[15,6]<-parameterEstimates(fit_sem_MD_superiorcorticostriate)[6,7]
data_subsample_WMlat_15_GGSEG[15,7]<-lavInspect(fit_sem_MD_superiorcorticostriate,what = "std.all")$beta[1,2]

##superiorcorticostriatefrontalcortex
Model_MD_superiorcorticostriatefrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                 Cognitive_factor ~ superiorcorticostriatefrontalcortex_md'
fit_sem_MD_superiorcorticostriatefrontalcortex <- sem(Model_MD_superiorcorticostriatefrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_superiorcorticostriatefrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[16,5]<-BIC(fit_sem_MD_superiorcorticostriatefrontalcortex)
data_subsample_WMlat_15_GGSEG[16,6]<-parameterEstimates(fit_sem_MD_superiorcorticostriatefrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[16,7]<-lavInspect(fit_sem_MD_superiorcorticostriatefrontalcortex,what = "std.all")$beta[1,2]

##superiorcorticostriateparietalcortex
Model_MD_superiorcorticostriateparietalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                  Cognitive_factor ~ superiorcorticostriateparietalcortex_md'
fit_sem_MD_superiorcorticostriateparietalcortex <- sem(Model_MD_superiorcorticostriateparietalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_superiorcorticostriateparietalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[17,5]<-BIC(fit_sem_MD_superiorcorticostriateparietalcortex)
data_subsample_WMlat_15_GGSEG[17,6]<-parameterEstimates(fit_sem_MD_superiorcorticostriateparietalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[17,7]<-lavInspect(fit_sem_MD_superiorcorticostriateparietalcortex,what = "std.all")$beta[1,2]

##striatalinferiorfrontalcortex
Model_MD_striatalinferiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                           Cognitive_factor ~ striatalinferiorfrontalcortex_md'
fit_sem_MD_striatalinferiorfrontalcortex <- sem(Model_MD_striatalinferiorfrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_striatalinferiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[18,5]<-BIC(fit_sem_MD_striatalinferiorfrontalcortex)
data_subsample_WMlat_15_GGSEG[18,6]<-parameterEstimates(fit_sem_MD_striatalinferiorfrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[18,7]<-lavInspect(fit_sem_MD_striatalinferiorfrontalcortex,what = "std.all")$beta[1,2]

##inferiorfrontalsuperiorfrontalcortex
Model_MD_inferiorfrontalsuperiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                  Cognitive_factor ~ inferiorfrontalsuperiorfrontalcortex_md'
fit_sem_MD_inferiorfrontalsuperiorfrontalcortex <- sem(Model_MD_inferiorfrontalsuperiorfrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[19,5]<-BIC(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex)
data_subsample_WMlat_15_GGSEG[19,6]<-parameterEstimates(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[19,7]<-lavInspect(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex,what = "std.all")$beta[1,2]

##fornix_exfimbria
Model_MD_fornix_exfimbria <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ fornix_exfimbria_md'
fit_sem_MD_fornix_exfimbria <- sem(Model_MD_fornix_exfimbria, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MD_fornix_exfimbria, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[20,5]<-BIC(fit_sem_MD_fornix_exfimbria)
data_subsample_WMlat_15_GGSEG[20,6]<-parameterEstimates(fit_sem_MD_fornix_exfimbria)[6,7]
data_subsample_WMlat_15_GGSEG[20,7]<-lavInspect(fit_sem_MD_fornix_exfimbria,what = "std.all")$beta[1,2]

#White Matter Volume (20 regions)
##fornix
Model_WMV_fornix <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                     Cognitive_factor ~ fornix_wmv'
fit_sem_WMV_fornix <- sem(Model_WMV_fornix, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_fornix, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[1,8]<-BIC(fit_sem_WMV_fornix)
data_subsample_WMlat_15_GGSEG[1,9]<-parameterEstimates(fit_sem_WMV_fornix)[6,7]
data_subsample_WMlat_15_GGSEG[1,10]<-lavInspect(fit_sem_WMV_fornix,what = "std.all")$beta[1,2]

##cingulatecingulum
Model_WMV_cingulatecingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ cingulatecingulum_wmv'
fit_sem_WMV_cingulatecingulum <- sem(Model_WMV_cingulatecingulum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_cingulatecingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[2,8]<-BIC(fit_sem_WMV_cingulatecingulum)
data_subsample_WMlat_15_GGSEG[2,9]<-parameterEstimates(fit_sem_WMV_cingulatecingulum)[6,7]
data_subsample_WMlat_15_GGSEG[2,10]<-lavInspect(fit_sem_WMV_cingulatecingulum,what = "std.all")$beta[1,2]

##parahippocampalcingulum
Model_WMV_parahippocampalcingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                      Cognitive_factor ~ parahippocampalcingulum_wmv'
fit_sem_WMV_parahippocampalcingulum <- sem(Model_WMV_parahippocampalcingulum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_parahippocampalcingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[3,8]<-BIC(fit_sem_WMV_parahippocampalcingulum)
data_subsample_WMlat_15_GGSEG[3,9]<-parameterEstimates(fit_sem_WMV_parahippocampalcingulum)[6,7]
data_subsample_WMlat_15_GGSEG[3,10]<-lavInspect(fit_sem_WMV_parahippocampalcingulum,what = "std.all")$beta[1,2]

##corticospinalpyramidal
Model_WMV_corticospinalpyramidal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ corticospinalpyramidal_wmv'
  fit_sem_WMV_corticospinalpyramidal <- sem(Model_WMV_corticospinalpyramidal, data=data_total_subsample15,estimator="mlr",missing="fiml")
  summary(fit_sem_WMV_corticospinalpyramidal, fit.measures=TRUE,rsquare=T,standardized=T)
  data_subsample_WMlat_15_GGSEG[4,8]<-BIC(fit_sem_WMV_corticospinalpyramidal)
  data_subsample_WMlat_15_GGSEG[4,9]<-parameterEstimates(fit_sem_WMV_corticospinalpyramidal)[6,7]
data_subsample_WMlat_15_GGSEG[4,10]<-lavInspect(fit_sem_WMV_corticospinalpyramidal,what = "std.all")$beta[1,2]

##anteriorthalamicradiations
Model_WMV_anteriorthalamicradiations <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                         Cognitive_factor ~ anteriorthalamicradiations_wmv'
fit_sem_WMV_anteriorthalamicradiations <- sem(Model_WMV_anteriorthalamicradiations, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_anteriorthalamicradiations, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[5,8]<-BIC(fit_sem_WMV_anteriorthalamicradiations)
data_subsample_WMlat_15_GGSEG[5,9]<-parameterEstimates(fit_sem_WMV_anteriorthalamicradiations)[6,7]
data_subsample_WMlat_15_GGSEG[5,10]<-lavInspect(fit_sem_WMV_anteriorthalamicradiations,what = "std.all")$beta[1,2]

##uncinate
Model_WMV_uncinate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                       Cognitive_factor ~ uncinate_wmv'
fit_sem_WMV_uncinate <- sem(Model_WMV_uncinate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_uncinate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[6,8]<-BIC(fit_sem_WMV_uncinate)
data_subsample_WMlat_15_GGSEG[6,9]<-parameterEstimates(fit_sem_WMV_uncinate)[6,7]
data_subsample_WMlat_15_GGSEG[6,10]<-lavInspect(fit_sem_WMV_uncinate,what = "std.all")$beta[1,2]

##inferiorlongitudinalfasiculus
Model_WMV_inferiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                            Cognitive_factor ~ inferiorlongitudinalfasiculus_wmv'
fit_sem_WMV_inferiorlongitudinalfasiculus <- sem(Model_WMV_inferiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_inferiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[7,8]<-BIC(fit_sem_WMV_inferiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[7,9]<-parameterEstimates(fit_sem_WMV_inferiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[7,10]<-lavInspect(fit_sem_WMV_inferiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##inferiorfrontooccipitalfasiculus
Model_WMV_inferiorfrontooccipitalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                               Cognitive_factor ~ inferiorfrontooccipitalfasiculus_wmv'
fit_sem_WMV_inferiorfrontooccipitalfasiculus <- sem(Model_WMV_inferiorfrontooccipitalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_inferiorfrontooccipitalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[8,8]<-BIC(fit_sem_WMV_inferiorfrontooccipitalfasiculus)
data_subsample_WMlat_15_GGSEG[8,9]<-parameterEstimates(fit_sem_WMV_inferiorfrontooccipitalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[8,10]<-lavInspect(fit_sem_WMV_inferiorfrontooccipitalfasiculus,what = "std.all")$beta[1,2]

##forcepsmajor
Model_WMV_forcepsmajor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ forcepsmajor_wmv'
fit_sem_WMV_forcepsmajor <- sem(Model_WMV_forcepsmajor, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_forcepsmajor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[9,8]<-BIC(fit_sem_WMV_forcepsmajor)
data_subsample_WMlat_15_GGSEG[9,9]<-parameterEstimates(fit_sem_WMV_forcepsmajor)[6,7]
data_subsample_WMlat_15_GGSEG[9,10]<-lavInspect(fit_sem_WMV_forcepsmajor,what = "std.all")$beta[1,2]

##forcepsminor
Model_WMV_forcepsminor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ forcepsminor_wmv'
fit_sem_WMV_forcepsminor <- sem(Model_WMV_forcepsminor, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_forcepsminor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[10,8]<-BIC(fit_sem_WMV_forcepsminor)
data_subsample_WMlat_15_GGSEG[10,9]<-parameterEstimates(fit_sem_WMV_forcepsminor)[6,7]
data_subsample_WMlat_15_GGSEG[10,10]<-lavInspect(fit_sem_WMV_forcepsminor,what = "std.all")$beta[1,2]

##corpuscallosum
Model_WMV_corpuscallosum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ corpuscallosum_wmv'
fit_sem_WMV_corpuscallosum <- sem(Model_WMV_corpuscallosum, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_corpuscallosum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[11,8]<-BIC(fit_sem_WMV_corpuscallosum)
data_subsample_WMlat_15_GGSEG[11,9]<-parameterEstimates(fit_sem_WMV_corpuscallosum)[6,7]
data_subsample_WMlat_15_GGSEG[11,10]<-lavInspect(fit_sem_WMV_corpuscallosum,what = "std.all")$beta[1,2]

##superiorlongitudinalfasiculus
Model_WMV_superiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                            Cognitive_factor ~ superiorlongitudinalfasiculus_wmv'
fit_sem_WMV_superiorlongitudinalfasiculus <- sem(Model_WMV_superiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_superiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[12,8]<-BIC(fit_sem_WMV_superiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[12,9]<-parameterEstimates(fit_sem_WMV_superiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[12,10]<-lavInspect(fit_sem_WMV_superiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##temporalsuperiorlongitudinalfasiculus
Model_WMV_temporalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                    Cognitive_factor ~ temporalsuperiorlongitudinalfasiculus_wmv'
fit_sem_WMV_temporalsuperiorlongitudinalfasiculus <- sem(Model_WMV_temporalsuperiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[13,8]<-BIC(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[13,9]<-parameterEstimates(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[13,10]<-lavInspect(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##parietalsuperiorlongitudinalfasiculus
Model_WMV_parietalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                    Cognitive_factor ~ parietalsuperiorlongitudinalfasiculus_wmv'
fit_sem_WMV_parietalsuperiorlongitudinalfasiculus <- sem(Model_WMV_parietalsuperiorlongitudinalfasiculus, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[14,8]<-BIC(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_15_GGSEG[14,9]<-parameterEstimates(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_15_GGSEG[14,10]<-lavInspect(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]

##superiorcorticostriate
Model_WMV_superiorcorticostriate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ superiorcorticostriate_wmv'
fit_sem_WMV_superiorcorticostriate <- sem(Model_WMV_superiorcorticostriate, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_superiorcorticostriate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[15,8]<-BIC(fit_sem_WMV_superiorcorticostriate)
data_subsample_WMlat_15_GGSEG[15,9]<-parameterEstimates(fit_sem_WMV_superiorcorticostriate)[6,7]
data_subsample_WMlat_15_GGSEG[15,10]<-lavInspect(fit_sem_WMV_superiorcorticostriate,what = "std.all")$beta[1,2]

##superiorcorticostriatefrontalcortex
Model_WMV_superiorcorticostriatefrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                  Cognitive_factor ~ superiorcorticostriatefrontalcortex_wmv'
fit_sem_WMV_superiorcorticostriatefrontalcortex <- sem(Model_WMV_superiorcorticostriatefrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_superiorcorticostriatefrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[16,8]<-BIC(fit_sem_WMV_superiorcorticostriatefrontalcortex)
data_subsample_WMlat_15_GGSEG[16,9]<-parameterEstimates(fit_sem_WMV_superiorcorticostriatefrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[16,10]<-lavInspect(fit_sem_WMV_superiorcorticostriatefrontalcortex,what = "std.all")$beta[1,2]

##superiorcorticostriateparietalcortex
Model_WMV_superiorcorticostriateparietalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                   Cognitive_factor ~ superiorcorticostriateparietalcortex_wmv'
fit_sem_WMV_superiorcorticostriateparietalcortex <- sem(Model_WMV_superiorcorticostriateparietalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_superiorcorticostriateparietalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[17,8]<-BIC(fit_sem_WMV_superiorcorticostriateparietalcortex)
data_subsample_WMlat_15_GGSEG[17,9]<-parameterEstimates(fit_sem_WMV_superiorcorticostriateparietalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[17,10]<-lavInspect(fit_sem_WMV_superiorcorticostriateparietalcortex,what = "std.all")$beta[1,2]

##striatalinferiorfrontalcortex
Model_WMV_striatalinferiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                            Cognitive_factor ~ striatalinferiorfrontalcortex_wmv'
fit_sem_WMV_striatalinferiorfrontalcortex <- sem(Model_WMV_striatalinferiorfrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_striatalinferiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[18,8]<-BIC(fit_sem_WMV_striatalinferiorfrontalcortex)
data_subsample_WMlat_15_GGSEG[18,9]<-parameterEstimates(fit_sem_WMV_striatalinferiorfrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[18,10]<-lavInspect(fit_sem_WMV_striatalinferiorfrontalcortex,what = "std.all")$beta[1,2]

##inferiorfrontalsuperiorfrontalcortex
Model_WMV_inferiorfrontalsuperiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                                   Cognitive_factor ~ inferiorfrontalsuperiorfrontalcortex_wmv'
fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex <- sem(Model_WMV_inferiorfrontalsuperiorfrontalcortex, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[19,8]<-BIC(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex)
data_subsample_WMlat_15_GGSEG[19,9]<-parameterEstimates(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex)[6,7]
data_subsample_WMlat_15_GGSEG[19,10]<-lavInspect(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex,what = "std.all")$beta[1,2]

##fornix_exfimbria
Model_WMV_fornix_exfimbria <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                               Cognitive_factor ~ fornix_exfimbria_wmv'
fit_sem_WMV_fornix_exfimbria <- sem(Model_WMV_fornix_exfimbria, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_fornix_exfimbria, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_15_GGSEG[20,8]<-BIC(fit_sem_WMV_fornix_exfimbria)
data_subsample_WMlat_15_GGSEG[20,9]<-parameterEstimates(fit_sem_WMV_fornix_exfimbria)[6,7]
data_subsample_WMlat_15_GGSEG[20,10]<-lavInspect(fit_sem_WMV_fornix_exfimbria,what = "std.all")$beta[1,2]


#Data Visualization for the individual models in the three white matter metrics
#Select all the std for the white matter metrics
data_subsample_WMlat_15_GGSEG_FA <- data_subsample_WMlat_15_GGSEG %>%
  dplyr::select(-c("FA_fit","FA_p","MD_fit","MD_p","MD_std","WMV_fit","WMV_p","WMV_std")) %>%
  rename(c("FA_std" = "Std")) %>%
  add_column(Metric = "FA")

data_subsample_WMlat_15_GGSEG_MD <- data_subsample_WMlat_15_GGSEG %>%
  dplyr::select(-c("FA_fit","FA_p","FA_std","MD_fit","MD_p","WMV_fit","WMV_p","WMV_std")) %>%
  rename(c("MD_std" = "Std")) %>%
  add_column(Metric = "MD")

data_subsample_WMlat_15_GGSEG_WMV <- data_subsample_WMlat_15_GGSEG %>%
  dplyr::select(-c("FA_fit","FA_p","FA_std","MD_fit","MD_p","MD_std","WMV_fit","WMV_p")) %>%
  rename(c("WMV_std" = "Std")) %>%
  add_column(Metric = "WMV")

#Compute the maximum and the minimum
min(data_subsample_WMlat_15_GGSEG_FA[,2])
max(data_subsample_WMlat_15_GGSEG_FA[,2])
min(data_subsample_WMlat_15_GGSEG_MD[,2])
max(data_subsample_WMlat_15_GGSEG_MD[,2])
min(data_subsample_WMlat_15_GGSEG_WMV[,2])
max(data_subsample_WMlat_15_GGSEG_WMV[,2])

data_subsample_STD_WMlat_15 <- rbind(data_subsample_WMlat_15_GGSEG_FA,data_subsample_WMlat_15_GGSEG_MD,by=c("label"))
data_subsample_STD_WMlat_15 <- rbind(data_subsample_STD_WMlat_15,data_subsample_WMlat_15_GGSEG_WMV,by=c("label"))
data_subsample_STD_WMlat_15<- data_subsample_STD_WMlat_15[-c(41,62),] 

data_subsample_STD_WMlat_15 <- data_subsample_STD_WMlat_15 %>% group_by(Metric) %>%  mutate(mean = mean(as.numeric(Std)))

ggplot(data_subsample_STD_WMlat_15,aes(as.numeric(Std),fill=Metric)) + 
  geom_histogram(binwidth=0.005)+
  geom_density(adjust=2,alpha=0.2) +
  facet_grid(Metric~.) +
  geom_vline(aes(xintercept = mean, group = Metric), linetype="dashed", size=1,alpha=0.4) +
  theme_classic(base_size = 30)+
  xlab('Standardized estimated model parameters') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_fill_viridis(discrete=TRUE, option="mako",begin=0.3,end=0.9,name="Metrics",labels=c("Fractional Anisotropy","Mean Diffusivity","White Matter Volume"))
ggsave("plot_WMlat_15_std.png")

#------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Model per metric - fractional anisotropy
#------------------------------------------------------------------------------------------------------------------------------
  
#Lasso regularization of the fractional anisotropy measure from all 20 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_fa <- data_total_subsample15[,c(112:114,121:137,172)]
df_fa <- makeX(df_fa, na.impute = TRUE)
x_fa <- df_fa[,1:20]
y_fa <- df_fa[,21]
cv_lasso_FAlat <- cv.glmnet(x_fa,y_fa, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_FAlat)
new_coef_fa <- coef(cv_lasso_FAlat, s = "lambda.min")
new_coef_fa <- as.data.frame(summary(new_coef_fa))
new_coef_fa <- new_coef_fa[,c(1,3)]
names(new_coef_fa) <- c("index","coef")

ROIs_coef_fa <- data.frame(label=c("forcepsmajor_fa","forcepsminor_fa","corpuscallosum_fa","fornix_fa","cingulatecingulum_fa","parahippocampalcingulum_fa","corticospinalpyramidal_fa","anteriorthalamicradiations_fa","uncinate_fa","inferiorlongitudinalfasiculus_fa","inferiorfrontooccipitalfasiculus_fa","superiorlongitudinalfasiculus_fa","temporalsuperiorlongitudinalfasiculus_fa","parietalsuperiorlongitudinalfasiculus_fa","superiorcorticostriate_fa","superiorcorticostriatefrontalcortex_fa","superiorcorticostriateparietalcortex_fa","striatalinferiorfrontalcortex_fa","inferiorfrontalsuperiorfrontalcortex_fa","fornix_exfimbria_fa"),
                           coef=rep("NA",times=20))

ROIs_coef_fa[new_coef_fa[1,1],2] <- new_coef_fa[1,2]

for (i in 1:nrow(new_coef_fa)) {
  ROIs_coef_fa[new_coef_fa[i,1],2] <- new_coef_fa[i,2]
}
ROI_fa<-c()
for (i in 1:20) {
  if (ROIs_coef_fa[i,2] != "NA") ROI_fa[i]<- ROIs_coef_fa[i,1]
}
ROI_fa<-na.omit(ROI_fa)
ROI_fa_plus<-c(ROI_fa[1])      #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_fa)) {
  ROI_fa_plus<-paste(ROI_fa_plus,"+",ROI_fa[i])
}
cat(ROI_fa_plus)
  
#Regularized fractional anisotropy model
Model_FAlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ forcepsminor_fa + corpuscallosum_fa + cingulatecingulum_fa + corticospinalpyramidal_fa + anteriorthalamicradiations_fa + inferiorfrontooccipitalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa'

fit_sem_FAlat_15_reg_free <- sem(Model_FAlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FAlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized fractional anisotropy model
Model_FAlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*forcepsminor_fa + a*corpuscallosum_fa + a*cingulatecingulum_fa + a*corticospinalpyramidal_fa + a*anteriorthalamicradiations_fa + a*inferiorfrontooccipitalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa'

fit_sem_FAlat_15_reg_constrained <- sem(Model_FAlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FAlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_FAlat_15_reg_free,fit_sem_FAlat_15_reg_constrained)

#Table for the regularized fractional anisotropy model
BIC(fit_sem_FAlat_15_reg_free)
data_subsample_FAlat_15_rl_reg_GGSEG <- data.frame(label=c("forcepsminor","corpuscallosum","lh_cingulatecingulum","rh_cingulatecingulum","lh_corticospinalpyramidal","rh_corticospinalpyramidal","lh_anteriorthalamicradiations","rh_anteriorthalamicradiations","lh_inferiorfrontooccipitalfasiculus","rh_inferiorfrontooccipitalfasiculus","lh_temporalsuperiorlongitudinalfasiculus","rh_temporalsuperiorlongitudinalfasiculus"),
                                                   FA_p=rep(NA, times=12),
                                                   FA_std=rep(NA, times=12))
data_subsample_FAlat_15_rl_reg_GGSEG[1,2]<-parameterEstimates(fit_sem_FAlat_15_reg_free)[6,7]
data_subsample_FAlat_15_rl_reg_GGSEG[1,3]<-lavInspect(fit_sem_FAlat_15_reg_free,what = "std.all")$beta[1,2]
data_subsample_FAlat_15_rl_reg_GGSEG[2,2]<-parameterEstimates(fit_sem_FAlat_15_reg_free)[7,7]
data_subsample_FAlat_15_rl_reg_GGSEG[2,3]<-lavInspect(fit_sem_FAlat_15_reg_free,what = "std.all")$beta[1,3]

a=3
for (i in seq(from = 3, to = 12, by = 2)) {
  data_subsample_FAlat_15_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_FAlat_15_reg_free)[a+5,7]
  data_subsample_FAlat_15_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_FAlat_15_reg_free)[a+5,7]
  data_subsample_FAlat_15_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_FAlat_15_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_FAlat_15_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_FAlat_15_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

#Model for fractional anisotropy with all the regions of interest (20 ROIs)
Model_FAlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ fornix_fa + cingulatecingulum_fa + parahippocampalcingulum_fa + corticospinalpyramidal_fa + anteriorthalamicradiations_fa + uncinate_fa + inferiorlongitudinalfasiculus_fa + inferiorfrontooccipitalfasiculus_fa + forcepsmajor_fa + forcepsminor_fa + corpuscallosum_fa + superiorlongitudinalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + parietalsuperiorlongitudinalfasiculus_fa + superiorcorticostriate_fa + superiorcorticostriatefrontalcortex_fa + superiorcorticostriateparietalcortex_fa + striatalinferiorfrontalcortex_fa + inferiorfrontalsuperiorfrontalcortex_fa + fornix_exfimbria_fa'

fit_sem_FAlat_15_all_free <- sem(Model_FAlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FAlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the fractional anisotropy model with all the regions
Model_FAlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*fornix_fa + a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*corticospinalpyramidal_fa + a*anteriorthalamicradiations_fa + a*uncinate_fa + a*inferiorlongitudinalfasiculus_fa + a*inferiorfrontooccipitalfasiculus_fa + a*forcepsmajor_fa + a*forcepsminor_fa + a*corpuscallosum_fa + a*superiorlongitudinalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*parietalsuperiorlongitudinalfasiculus_fa + a*superiorcorticostriate_fa + a*superiorcorticostriatefrontalcortex_fa + a*superiorcorticostriateparietalcortex_fa + a*striatalinferiorfrontalcortex_fa + a*inferiorfrontalsuperiorfrontalcortex_fa + a*fornix_exfimbria_fa'

fit_sem_FAlat_15_all_constrained <- sem(Model_FAlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_FAlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_FAlat_15_all_free,fit_sem_FAlat_15_all_constrained)

#-------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Model per metric - mean diffusivity
#-------------------------------------------------------------------------------------------------------------------------------
  
#Lasso regularization of the mean diffusivity measure from all 20 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_md <- data_total_subsample15[,c(115:117,138:154,172)]
df_md <- makeX(df_md, na.impute = TRUE)
x_md <- df_md[,1:20]
y_md <- df_md[,21]
cv_lasso_MDlat <- cv.glmnet(x_md,y_md, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_MDlat)
new_coef_md <- coef(cv_lasso_MDlat, s = "lambda.min")
new_coef_md <- as.data.frame(summary(new_coef_md))
new_coef_md <- new_coef_md[,c(1,3)]
names(new_coef_md) <- c("index","coef")

ROIs_coef_md <- data.frame(label=c("forcepsmajor_md","forcepsminor_md","corpuscallosum_md","fornix_md","cingulatecingulum_md","parahippocampalcingulum_md","corticospinalpyramidal_md","anteriorthalamicradiations_md","uncinate_md","inferiorlongitudinalfasiculus_md","inferiorfrontooccipitalfasiculus_md","superiorlongitudinalfasiculus_md","temporalsuperiorlongitudinalfasiculus_md","parietalsuperiorlongitudinalfasiculus_md","superiorcorticostriate_md","superiorcorticostriatefrontalcortex_md","superiorcorticostriateparietalcortex_md","striatalinferiorfrontalcortex_md","inferiorfrontalsuperiorfrontalcortex_md","fornix_exfimbria_md"),
                           coef=rep("NA",times=20))

ROIs_coef_md[new_coef_md[1,1],2] <- new_coef_md[1,2]

for (i in 1:nrow(new_coef_md)) {
  ROIs_coef_md[new_coef_md[i,1],2] <- new_coef_md[i,2]
}

ROI_md<-c()
for (i in 1:20) {
  if (ROIs_coef_md[i,2] != "NA") ROI_md[i]<- ROIs_coef_md[i,1]
}
ROI_md<-na.omit(ROI_md)
ROI_md_plus<-c(ROI_md[1])      #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_md)) {
  ROI_md_plus<-paste(ROI_md_plus,"+",ROI_md[i])
}
cat(ROI_md_plus)  
  
#Regularized mean diffusivity model
Model_MDlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ forcepsminor_md + corpuscallosum_md + parahippocampalcingulum_md + anteriorthalamicradiations_md + superiorcorticostriate_md'

fit_sem_MDlat_15_reg_free <- sem(Model_MDlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MDlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized mean diffusivity model
Model_MDlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*forcepsminor_md + a*corpuscallosum_md + a*parahippocampalcingulum_md + a*anteriorthalamicradiations_md + a*superiorcorticostriate_md'

fit_sem_MDlat_15_reg_constrained <- sem(Model_MDlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MDlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
anova(fit_sem_MDlat_15_reg_free,fit_sem_MDlat_15_reg_constrained)

#Table for the regularized mean diffusivity model
BIC(fit_sem_MDlat_15_reg_free)
data_subsample_MDlat_15_rl_reg_GGSEG <- data.frame(label=c("forcepsminor","corpuscallosum","lh_parahippocampalcingulum","rh_parahippocampalcingulum","lh_anteriorthalamicradiations","rh_anteriorthalamicradiations","lh_superiorcorticostriate","rh_superiorcorticostriate"),
                                                   MD_p=rep(NA, times=8),
                                                   MD_std=rep(NA, times=8))
data_subsample_MDlat_15_rl_reg_GGSEG[1,2]<-parameterEstimates(fit_sem_MDlat_15_reg_free)[6,7]
data_subsample_MDlat_15_rl_reg_GGSEG[1,3]<-lavInspect(fit_sem_MDlat_15_reg_free,what = "std.all")$beta[1,2]
data_subsample_MDlat_15_rl_reg_GGSEG[2,2]<-parameterEstimates(fit_sem_MDlat_15_reg_free)[7,7]
data_subsample_MDlat_15_rl_reg_GGSEG[2,3]<-lavInspect(fit_sem_MDlat_15_reg_free,what = "std.all")$beta[1,3]

a=3
for (i in seq(from = 3, to = 8, by = 2)) {
  data_subsample_MDlat_15_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_MDlat_15_reg_free)[a+5,7]
  data_subsample_MDlat_15_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_MDlat_15_reg_free)[a+5,7]
  data_subsample_MDlat_15_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_MDlat_15_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_MDlat_15_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_MDlat_15_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

#Model for mean diffusivity with all the regions of interest (20 ROIs)
Model_MDlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ fornix_md + cingulatecingulum_md + parahippocampalcingulum_md + corticospinalpyramidal_md + anteriorthalamicradiations_md + uncinate_md + inferiorlongitudinalfasiculus_md + inferiorfrontooccipitalfasiculus_md + forcepsmajor_md + forcepsminor_md + corpuscallosum_md + superiorlongitudinalfasiculus_md + temporalsuperiorlongitudinalfasiculus_md + parietalsuperiorlongitudinalfasiculus_md + superiorcorticostriate_md + superiorcorticostriatefrontalcortex_md + superiorcorticostriateparietalcortex_md + striatalinferiorfrontalcortex_md + inferiorfrontalsuperiorfrontalcortex_md + fornix_exfimbria_md'

fit_sem_MDlat_15_all_free <- sem(Model_MDlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MDlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the mean diffusivity model with all the regions
Model_MDlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*fornix_md + a*cingulatecingulum_md + a*parahippocampalcingulum_md + a*corticospinalpyramidal_md + a*anteriorthalamicradiations_md + a*uncinate_md + a*inferiorlongitudinalfasiculus_md + a*inferiorfrontooccipitalfasiculus_md + a*forcepsmajor_md + a*forcepsminor_md + a*corpuscallosum_md + a*superiorlongitudinalfasiculus_md + a*temporalsuperiorlongitudinalfasiculus_md + a*parietalsuperiorlongitudinalfasiculus_md + a*superiorcorticostriate_md + a*superiorcorticostriatefrontalcortex_md + a*superiorcorticostriateparietalcortex_md + a*striatalinferiorfrontalcortex_md + a*inferiorfrontalsuperiorfrontalcortex_md + a*fornix_exfimbria_md'

fit_sem_MDlat_15_all_constrained <- sem(Model_MDlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_MDlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_MDlat_15_all_free,fit_sem_MDlat_15_all_constrained)


#-------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Model per metric - white matter volume
#-------------------------------------------------------------------------------------------------------------------------------

#Lasso regularization of the white matter volume measure from all 20 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_wmv <- data_total_subsample15[,c(118:120,155:172)]
df_wmv <- makeX(df_wmv, na.impute = TRUE)
x_wmv <- df_wmv[,1:20]
y_wmv <- df_wmv[,21]
cv_lasso_WMVlat <- cv.glmnet(x_wmv,y_wmv, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_WMVlat)
new_coef_wmv <- coef(cv_lasso_WMVlat, s = "lambda.min")
new_coef_wmv <- as.data.frame(summary(new_coef_wmv))
new_coef_wmv <- new_coef_wmv[,c(1,3)]
names(new_coef_wmv) <- c("index","coef")

ROIs_coef_wmv <- data.frame(label=c("forcepsmajor_wmv","forcepsminor_wmv","corpuscallosum_wmv","fornix_wmv","cingulatecingulum_wmv","parahippocampalcingulum_wmv","corticospinalpyramidal_wmv","anteriorthalamicradiations_wmv","uncinate_wmv","inferiorlongitudinalfasiculus_wmv","inferiorfrontooccipitalfasiculus_wmv","superiorlongitudinalfasiculus_wmv","temporalsuperiorlongitudinalfasiculus_wmv","parietalsuperiorlongitudinalfasiculus_wmv","superiorcorticostriate_wmv","superiorcorticostriatefrontalcortex_wmv","superiorcorticostriateparietalcortex_wmv","striatalinferiorfrontalcortex_wmv","inferiorfrontalsuperiorfrontalcortex_wmv","fornix_exfimbria_wmv"),
                            coef=rep("NA",times=20))

ROIs_coef_wmv[new_coef_wmv[1,1],2] <- new_coef_wmv[1,2]

for (i in 1:nrow(new_coef_wmv)) {
  ROIs_coef_wmv[new_coef_wmv[i,1],2] <- new_coef_wmv[i,2]
}
ROI_wmv<-c()
for (i in 1:20) {
  if (ROIs_coef_wmv[i,2] != "NA") ROI_wmv[i]<- ROIs_coef_wmv[i,1]
}
ROI_wmv<-na.omit(ROI_wmv)
ROI_wmv_plus<-c(ROI_wmv[1])      #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_wmv)) {
  ROI_wmv_plus<-paste(ROI_wmv_plus,"+",ROI_wmv[i])
}
cat(ROI_wmv_plus)  

#Regularized white matter volume model
Model_WMVlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ forcepsminor_wmv + corpuscallosum_wmv + cingulatecingulum_wmv + parahippocampalcingulum_wmv + corticospinalpyramidal_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + superiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + superiorcorticostriateparietalcortex_wmv + striatalinferiorfrontalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv + fornix_exfimbria_wmv'

fit_sem_WMVlat_15_reg_free <- sem(Model_WMVlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMVlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)
 
#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized white matter volume model
Model_WMVlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ a*forcepsminor_wmv + a*corpuscallosum_wmv + a*cingulatecingulum_wmv + a*parahippocampalcingulum_wmv + a*corticospinalpyramidal_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*superiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*superiorcorticostriateparietalcortex_wmv + a*striatalinferiorfrontalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv + a*fornix_exfimbria_wmv'

fit_sem_WMVlat_15_reg_constrained <- sem(Model_WMVlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMVlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_WMVlat_15_reg_free,fit_sem_WMVlat_15_reg_constrained) 


#Table for the regularized white matter volume model
BIC(fit_sem_WMVlat_15_reg_free)
data_subsample_WMVlat_15_rl_reg_GGSEG <- data.frame(label=c("forcepsminor","corpuscallosum","lh_cingulatecingulum","rh_cingulatecingulum","lh_parahippocampalcingulum","rh_parahippocampalcingulum","lh_corticospinalpyramidal","rh_corticospinalpyramidal","lh_anteriorthalamicradiations","rh_anteriorthalamicradiations","lh_uncinate","rh_uncinate","lh_inferiorlongitudinalfasiculus","rh_inferiorlongitudinalfasiculus","lh_inferiorfrontooccipitalfasiculus","rh_inferiorfrontooccipitalfasiculus","lh_superiorlongitudinalfasiculus","rh_superiorlongitudinalfasiculus","lh_parietalsuperiorlongitudinalfasiculus","rh_parietalsuperiorlongitudinalfasiculus","lh_superiorcorticostriateparietalcortex","rh_superiorcorticostriateparietalcortex","lh_striatalinferiorfrontalcortex","rh_striatalinferiorfrontalcortex","lh_inferiorfrontalsuperiorfrontalcortex","rh_inferiorfrontalsuperiorfrontalcortex","lh_fornix_exfimbria","rh_fornix_exfimbria"),
                                                   MD_p=rep(NA, times=28),
                                                   MD_std=rep(NA, times=28))
data_subsample_WMVlat_15_rl_reg_GGSEG[1,2]<-parameterEstimates(fit_sem_WMVlat_15_reg_free)[6,7]
data_subsample_WMVlat_15_rl_reg_GGSEG[1,3]<-lavInspect(fit_sem_WMVlat_15_reg_free,what = "std.all")$beta[1,2]
data_subsample_WMVlat_15_rl_reg_GGSEG[2,2]<-parameterEstimates(fit_sem_WMVlat_15_reg_free)[7,7]
data_subsample_WMVlat_15_rl_reg_GGSEG[2,3]<-lavInspect(fit_sem_WMVlat_15_reg_free,what = "std.all")$beta[1,3]

a=3
for (i in seq(from = 3, to = 28, by = 2)) {
  data_subsample_WMVlat_15_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_WMVlat_15_reg_free)[a+5,7]
  data_subsample_WMVlat_15_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_WMVlat_15_reg_free)[a+5,7]
  data_subsample_WMVlat_15_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_WMVlat_15_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_WMVlat_15_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_WMVlat_15_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}


#Model for white matter volume with all the regions of interest (20 ROIs)
Model_WMVlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ fornix_wmv + cingulatecingulum_wmv + parahippocampalcingulum_wmv + corticospinalpyramidal_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + forcepsmajor_wmv + forcepsminor_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + temporalsuperiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + superiorcorticostriate_wmv + superiorcorticostriatefrontalcortex_wmv + superiorcorticostriateparietalcortex_wmv + striatalinferiorfrontalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv + fornix_exfimbria_wmv'

fit_sem_WMVlat_15_all_free <- sem(Model_WMVlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMVlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the white matter volume model with all the regions
Model_WMVlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ a*fornix_wmv + a*cingulatecingulum_wmv + a*parahippocampalcingulum_wmv + a*corticospinalpyramidal_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*forcepsmajor_wmv + a*forcepsminor_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*temporalsuperiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*superiorcorticostriate_wmv + a*superiorcorticostriatefrontalcortex_wmv + a*superiorcorticostriateparietalcortex_wmv + a*striatalinferiorfrontalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv + a*fornix_exfimbria_wmv'

fit_sem_WMVlat_15_all_constrained <- sem(Model_WMVlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMVlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_WMVlat_15_all_free,fit_sem_WMVlat_15_all_constrained)

#-------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Models including the three white matter metrics
#-------------------------------------------------------------------------------------------------------------------------------

#Lasso regularization of the white matter measures from all 60 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_whitem <- data_total_subsample15[,c(112:172)]
df_whitem <- makeX(df_whitem, na.impute = TRUE)
x_whitem <- df_whitem[,1:60]
y_whitem <- df_whitem[,61]
cv_lasso_whitemlat <- cv.glmnet(x_whitem,y_whitem, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_whitemlat)
new_coef_whitem <- coef(cv_lasso_whitemlat, s = "lambda.min")
new_coef_whitem <- as.data.frame(summary(new_coef_whitem))
new_coef_whitem <- new_coef_whitem[,c(1,3)]
names(new_coef_whitem) <- c("index","coef")

ROIs_coef_whitem <- data.frame(label=c("fornix_fa","cingulatecingulum_fa","parahippocampalcingulum_fa","corticospinalpyramidal_fa","anteriorthalamicradiations_fa","uncinate_fa","inferiorlongitudinalfasiculus_fa","inferiorfrontooccipitalfasiculus_fa","forcepsmajor_fa","forcepsminor_fa","corpuscallosum_fa","superiorlongitudinalfasiculus_fa","temporalsuperiorlongitudinalfasiculus_fa","parietalsuperiorlongitudinalfasiculus_fa","superiorcorticostriate_fa","superiorcorticostriatefrontalcortex_fa","superiorcorticostriateparietalcortex_fa","striatalinferiorfrontalcortex_fa","inferiorfrontalsuperiorfrontalcortex_fa","fornix_exfimbria_fa","fornix_md","cingulatecingulum_md","parahippocampalcingulum_md","corticospinalpyramidal_md","anteriorthalamicradiations_md","uncinate_md","inferiorlongitudinalfasiculus_md","inferiorfrontooccipitalfasiculus_md","forcepsmajor_md","forcepsminor_md","corpuscallosum_md","superiorlongitudinalfasiculus_md","temporalsuperiorlongitudinalfasiculus_md","parietalsuperiorlongitudinalfasiculus_md","superiorcorticostriate_md","superiorcorticostriatefrontalcortex_md","superiorcorticostriateparietalcortex_md","striatalinferiorfrontalcortex_md","inferiorfrontalsuperiorfrontalcortex_md","fornix_exfimbria_md","fornix_wmv","cingulatecingulum_wmv","parahippocampalcingulum_wmv","corticospinalpyramidal_wmv","anteriorthalamicradiations_wmv","uncinate_wmv","inferiorlongitudinalfasiculus_wmv","inferiorfrontooccipitalfasiculus_wmv","forcepsmajor_wmv","forcepsminor_wmv","corpuscallosum_wmv","superiorlongitudinalfasiculus_wmv","temporalsuperiorlongitudinalfasiculus_wmv","parietalsuperiorlongitudinalfasiculus_wmv","superiorcorticostriate_wmv","superiorcorticostriatefrontalcortex_wmv","superiorcorticostriateparietalcortex_wmv","striatalinferiorfrontalcortex_wmv","inferiorfrontalsuperiorfrontalcortex_wmv","fornix_exfimbria_wmv"),
                               coef=rep("NA",times=60))

ROIs_coef_whitem[new_coef_whitem[1,1],2] <- new_coef_whitem[1,2]

for (i in 1:nrow(new_coef_whitem)) {
  ROIs_coef_whitem[new_coef_whitem[i,1],2] <- new_coef_whitem[i,2]
}
ROI_whitem<-c()
for (i in 1:60) {
  if (ROIs_coef_whitem[i,2] != "NA") ROI_whitem[i]<- ROIs_coef_whitem[i,1]
}
ROI_whitem<-na.omit(ROI_whitem)
ROI_whitem_plus<-c(ROI_whitem[1])     #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_whitem)) {
  ROI_whitem_plus<-paste(ROI_whitem_plus,"+",ROI_whitem[i])
}
cat(ROI_whitem_plus)

#Regularized white matter metrics model 
Model_WMlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ cingulatecingulum_fa + parahippocampalcingulum_fa + uncinate_fa + inferiorfrontooccipitalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + inferiorfrontalsuperiorfrontalcortex_fa + fornix_exfimbria_fa + uncinate_md + inferiorfrontooccipitalfasiculus_md + forcepsmajor_md + forcepsminor_md + striatalinferiorfrontalcortex_md + fornix_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + forcepsmajor_wmv + forcepsminor_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_WMlat_15_reg_free <- sem(Model_WMlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized white matter metrics model
Model_WMlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*uncinate_fa + a*inferiorfrontooccipitalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*inferiorfrontalsuperiorfrontalcortex_fa + a*fornix_exfimbria_fa + a*uncinate_md + a*inferiorfrontooccipitalfasiculus_md + a*forcepsmajor_md + a*forcepsminor_md + a*striatalinferiorfrontalcortex_md + a*fornix_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*forcepsmajor_wmv + a*forcepsminor_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_WMlat_15_reg_constrained <- sem(Model_WMlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_WMlat_15_reg_free,fit_sem_WMlat_15_reg_constrained)

#Table for the regularized white matter metrics model
BIC(fit_sem_WMlat_15_reg_free)
data_subsample_WMlat_15_regFA_GGSEG <- data.frame(label=c("lh_cingulatecingulum","rh_cingulatecingulum","lh_parahippocampalcingulum","rh_parahippocampalcingulum","lh_uncinate","rh_uncinate","lh_inferiorfrontooccipitalfasiculus","rh_inferiorfrontooccipitalfasiculus","lh_temporalsuperiorlongitudinalfasiculus","rh_temporalsuperiorlongitudinalfasiculus","lh_inferiorfrontalsuperiorfrontalcortex","rh_inferiorfrontalsuperiorfrontalcortex","lh_fornix_exfimbria","rh_fornix_exfimbria"),
                                                    FA_p=rep(NA, times=14),
                                                    FA_std=rep(NA, times=14))
a=1
for (i in seq(from = 1, to = 14, by = 2)) {
  data_subsample_WMlat_15_regFA_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[a+5,7]
  data_subsample_WMlat_15_regFA_GGSEG[i+1,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[a+5,7]
  data_subsample_WMlat_15_regFA_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_WMlat_15_regFA_GGSEG[i+1,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_WMlat_15_regMD_GGSEG <- data.frame(label=c("lh_uncinate","rh_uncinate","lh_inferiorfrontooccipitalfasiculus","rh_inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","lh_striatalinferiorfrontalcortex","rh_striatalinferiorfrontalcortex"),
                                                    MD_p=rep(NA, times=8),
                                                    MD_std=rep(NA, times=8))
a=1
for (i in seq(from = 1, to = 4, by = 2)) {
  data_subsample_WMlat_15_regMD_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[a+12,7]
  data_subsample_WMlat_15_regMD_GGSEG[i+1,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[a+12,7]
  data_subsample_WMlat_15_regMD_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,a+8]
  data_subsample_WMlat_15_regMD_GGSEG[i+1,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,a+8]
  a=a+1
}

data_subsample_WMlat_15_regMD_GGSEG[5,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[15,7]
data_subsample_WMlat_15_regMD_GGSEG[5,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,11]
data_subsample_WMlat_15_regMD_GGSEG[6,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[16,7]
data_subsample_WMlat_15_regMD_GGSEG[6,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,12]
data_subsample_WMlat_15_regMD_GGSEG[7,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[17,7]
data_subsample_WMlat_15_regMD_GGSEG[8,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[17,7]
data_subsample_WMlat_15_regMD_GGSEG[7,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,13]
data_subsample_WMlat_15_regMD_GGSEG[8,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,13]


data_subsample_WMlat_15_regWMV_GGSEG <- data.frame(label=c("lh_fornix","rh_fornix","lh_anteriorthalamicradiations","rh_anteriorthalamicradiations","lh_uncinate","rh_uncinate","lh_inferiorlongitudinalfasiculus","rh_inferiorlongitudinalfasiculus","lh_inferiorfrontooccipitalfasiculus","rh_inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","corpuscallosum","lh_superiorlongitudinalfasiculus","rh_superiorlongitudinalfasiculus","lh_parietalsuperiorlongitudinalfasiculus","rh_parietalsuperiorlongitudinalfasiculus","lh_inferiorfrontalsuperiorfrontalcortex","rh_inferiorfrontalsuperiorfrontalcortex"),
                                                     WMV_p=rep(NA, times=19),
                                                     WMV_std=rep(NA, times=19))
a=1
for (i in seq(from = 1, to = 10, by = 2)) {
  data_subsample_WMlat_15_regWMV_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[a+17,7]
  data_subsample_WMlat_15_regWMV_GGSEG[i+1,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[a+17,7]
  data_subsample_WMlat_15_regWMV_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,a+13]
  data_subsample_WMlat_15_regWMV_GGSEG[i+1,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,a+13]
  a=a+1
}
data_subsample_WMlat_15_regWMV_GGSEG[11,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[23,7]
data_subsample_WMlat_15_regWMV_GGSEG[11,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,19]
data_subsample_WMlat_15_regWMV_GGSEG[12,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[24,7]
data_subsample_WMlat_15_regWMV_GGSEG[12,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,20]
data_subsample_WMlat_15_regWMV_GGSEG[13,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[25,7]
data_subsample_WMlat_15_regWMV_GGSEG[13,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,21]

a=9
for (i in seq(from = 14, to = 19, by = 2)) {
  data_subsample_WMlat_15_regWMV_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[a+17,7]
  data_subsample_WMlat_15_regWMV_GGSEG[i+1,2]<-parameterEstimates(fit_sem_WMlat_15_reg_free)[a+17,7]
  data_subsample_WMlat_15_regWMV_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,a+13]
  data_subsample_WMlat_15_regWMV_GGSEG[i+1,3]<-lavInspect(fit_sem_WMlat_15_reg_free,what = "std.all")$beta[1,a+13]
  a=a+1
}


#Model for white matter metrics with all the regions of interest (60 ROIs)
Model_WMlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                         Cognitive_factor ~ fornix_fa + cingulatecingulum_fa + parahippocampalcingulum_fa + corticospinalpyramidal_fa + anteriorthalamicradiations_fa + uncinate_fa + inferiorlongitudinalfasiculus_fa + inferiorfrontooccipitalfasiculus_fa + forcepsmajor_fa + forcepsminor_fa + corpuscallosum_fa + superiorlongitudinalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + parietalsuperiorlongitudinalfasiculus_fa + superiorcorticostriate_fa + superiorcorticostriatefrontalcortex_fa + superiorcorticostriateparietalcortex_fa + striatalinferiorfrontalcortex_fa + inferiorfrontalsuperiorfrontalcortex_fa + fornix_exfimbria_fa + fornix_md + cingulatecingulum_md + parahippocampalcingulum_md + corticospinalpyramidal_md + anteriorthalamicradiations_md + uncinate_md + inferiorlongitudinalfasiculus_md + inferiorfrontooccipitalfasiculus_md + forcepsmajor_md + forcepsminor_md + corpuscallosum_md + superiorlongitudinalfasiculus_md + temporalsuperiorlongitudinalfasiculus_md + parietalsuperiorlongitudinalfasiculus_md + superiorcorticostriate_md + superiorcorticostriatefrontalcortex_md + superiorcorticostriateparietalcortex_md + striatalinferiorfrontalcortex_md + inferiorfrontalsuperiorfrontalcortex_md + fornix_exfimbria_md + fornix_wmv + cingulatecingulum_wmv + parahippocampalcingulum_wmv + corticospinalpyramidal_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + forcepsmajor_wmv + forcepsminor_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + temporalsuperiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + superiorcorticostriate_wmv + superiorcorticostriatefrontalcortex_wmv + superiorcorticostriateparietalcortex_wmv + striatalinferiorfrontalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv + fornix_exfimbria_wmv'

fit_sem_WMlat_15_all_free <- sem(Model_WMlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the white matter metrics model with all the regions
Model_WMlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ a*fornix_fa + a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*corticospinalpyramidal_fa + a*anteriorthalamicradiations_fa + a*uncinate_fa + a*inferiorlongitudinalfasiculus_fa + a*inferiorfrontooccipitalfasiculus_fa + a*forcepsmajor_fa + a*forcepsminor_fa + a*corpuscallosum_fa + a*superiorlongitudinalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*parietalsuperiorlongitudinalfasiculus_fa + a*superiorcorticostriate_fa + a*superiorcorticostriatefrontalcortex_fa + a*superiorcorticostriateparietalcortex_fa + a*striatalinferiorfrontalcortex_fa + a*inferiorfrontalsuperiorfrontalcortex_fa + a*fornix_exfimbria_fa + a*fornix_md + a*cingulatecingulum_md + a*parahippocampalcingulum_md + a*corticospinalpyramidal_md + a*anteriorthalamicradiations_md + a*uncinate_md + a*inferiorlongitudinalfasiculus_md + a*inferiorfrontooccipitalfasiculus_md + a*forcepsmajor_md + a*forcepsminor_md + a*corpuscallosum_md + a*superiorlongitudinalfasiculus_md + a*temporalsuperiorlongitudinalfasiculus_md + a*parietalsuperiorlongitudinalfasiculus_md + a*superiorcorticostriate_md + a*superiorcorticostriatefrontalcortex_md + a*superiorcorticostriateparietalcortex_md + a*striatalinferiorfrontalcortex_md + a*inferiorfrontalsuperiorfrontalcortex_md + a*fornix_exfimbria_md + a*fornix_wmv + a*cingulatecingulum_wmv + a*parahippocampalcingulum_wmv + a*corticospinalpyramidal_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*forcepsmajor_wmv + a*forcepsminor_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*temporalsuperiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*superiorcorticostriate_wmv + a*superiorcorticostriatefrontalcortex_wmv + a*superiorcorticostriateparietalcortex_wmv + a*striatalinferiorfrontalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv + a*fornix_exfimbria_wmv'

fit_sem_WMlat_15_all_constrained <- sem(Model_WMlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_WMlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_WMlat_15_all_free,fit_sem_WMlat_15_all_constrained)

#Compare the information provide by the regularized regions in the three white matter metrics 
data_subsample_WMlat_15_reg_GGSEG<-merge(data_subsample_WMlat_15_regFA_GGSEG[,c(1,3)],data_subsample_WMlat_15_regMD_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_WMlat_15_reg_GGSEG<-merge(data_subsample_WMlat_15_reg_GGSEG,data_subsample_WMlat_15_regWMV_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_WMlat_15_reg_GGSEG, "STD_WMlat_15_reg.csv", row.names=FALSE)

data_subsample_WMlat_15_reg_GGSEG_long <- data_subsample_WMlat_15_reg_GGSEG %>%
  rename(c("FA_std" = "FA")) %>%
  rename(c("MD_std" = "MD")) %>%
  rename(c("WMV_std" = "WMV")) %>%
  gather(Metric, STD, FA, MD, WMV) 
data_subsample_WMlat_15_reg_GGSEG_long$Metric<-factor(data_subsample_WMlat_15_reg_GGSEG_long$Metric,levels=c("FA","MD","WMV"))

#-------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Models including grey & white matter metrics
#-------------------------------------------------------------------------------------------------------------------------------

#Lasso regularization of the white matter measures from all 60 regions of interest as simultaneously predicting cognitive performance
set.seed(2)
lambdas <- 10^seq(3, -2, by = -.1)
df_gmwm <- data_total_subsample15[,10:172]
df_gmwm <- makeX(df_gmwm, na.impute = TRUE)
x_gmwm <- df_gmwm[,1:162]
y_gmwm <- df_gmwm[,163]
cv_lasso_gmwmlat <- cv.glmnet(x_gmwm,y_gmwm, alpha = 1, lambda = lambdas, intercept = F, grouped=F)
plot(cv_lasso_gmwmlat)
new_coef_gmwm <- coef(cv_lasso_gmwmlat, s = "lambda.min")
new_coef_gmwm <- as.data.frame(summary(new_coef_gmwm))
new_coef_gmwm <- new_coef_gmwm[,c(1,3)]
names(new_coef_gmwm) <- c("index","coef")


ROIs_coef_gmwm <- data.frame(label=c("bankssts_ct","caudalanteriorcingulate_ct","caudalmiddlefrontal_ct","cuneus_ct","entorhinal_ct","fusiform_ct","inferiorparietal_ct","inferiortemporal_ct","isthmuscingulate_ct","lateraloccipital_ct","lateralorbitofrontal_ct","lingual_ct","medialorbitofrontal_ct","middletemporal_ct","parahippocampal_ct","paracentral_ct","parsopercularis_ct","parsorbitalis_ct","parstriangularis_ct","pericalcarine_ct","postcentral_ct","posteriorcingulate_ct","precentral_ct","precuneus_ct","rostralanteriorcingulate_ct","rostralmiddlefrontal_ct","superiorfrontal_ct","superiorparietal_ct","superiortemporal_ct","supramarginal_ct","frontalpole_ct","temporalpole_ct","transversetemporal_ct","insula_ct","bankssts_sa","caudalanteriorcingulate_sa","caudalmiddlefrontal_sa","cuneus_sa","entorhinal_sa","fusiform_sa","inferiorparietal_sa","inferiortemporal_sa","isthmuscingulate_sa","lateraloccipital_sa","lateralorbitofrontal_sa","lingual_sa","medialorbitofrontal_sa","middletemporal_sa","parahippocampal_sa","paracentral_sa","parsopercularis_sa","parsorbitalis_sa","parstriangularis_sa","pericalcarine_sa","postcentral_sa","posteriorcingulate_sa","precentral_sa","precuneus_sa","rostralanteriorcingulate_sa","rostralmiddlefrontal_sa","superiorfrontal_sa","superiorparietal_sa","superiortemporal_sa","supramarginal_sa","frontalpole_sa","temporalpole_sa","transversetemporal_sa","insula_sa","bankssts_gmv","caudalanteriorcingulate_gmv","caudalmiddlefrontal_gmv","cuneus_gmv","entorhinal_gmv","fusiform_gmv","inferiorparietal_gmv","inferiortemporal_gmv","isthmuscingulate_gmv","lateraloccipital_gmv","lateralorbitofrontal_gmv","lingual_gmv","medialorbitofrontal_gmv","middletemporal_gmv","parahippocampal_gmv","paracentral_gmv","parsopercularis_gmv","parsorbitalis_gmv","parstriangularis_gmv","pericalcarine_gmv","postcentral_gmv","posteriorcingulate_gmv","precentral_gmv","precuneus_gmv","rostralanteriorcingulate_gmv","rostralmiddlefrontal_gmv","superiorfrontal_gmv","superiorparietal_gmv","superiortemporal_gmv","supramarginal_gmv","frontalpole_gmv","temporalpole_gmv","transversetemporal_gmv","insula_gmv","fornix_fa","cingulatecingulum_fa","parahippocampalcingulum_fa","corticospinalpyramidal_fa","anteriorthalamicradiations_fa","uncinate_fa","inferiorlongitudinalfasiculus_fa","inferiorfrontooccipitalfasiculus_fa","forcepsmajor_fa","forcepsminor_fa","corpuscallosum_fa","superiorlongitudinalfasiculus_fa","temporalsuperiorlongitudinalfasiculus_fa","parietalsuperiorlongitudinalfasiculus_fa","superiorcorticostriate_fa","superiorcorticostriatefrontalcortex_fa","superiorcorticostriateparietalcortex_fa","striatalinferiorfrontalcortex_fa","inferiorfrontalsuperiorfrontalcortex_fa","fornix_exfimbria_fa","fornix_md","cingulatecingulum_md","parahippocampalcingulum_md","corticospinalpyramidal_md","anteriorthalamicradiations_md","uncinate_md","inferiorlongitudinalfasiculus_md","inferiorfrontooccipitalfasiculus_md","forcepsmajor_md","forcepsminor_md","corpuscallosum_md","superiorlongitudinalfasiculus_md","temporalsuperiorlongitudinalfasiculus_md","parietalsuperiorlongitudinalfasiculus_md","superiorcorticostriate_md","superiorcorticostriatefrontalcortex_md","superiorcorticostriateparietalcortex_md","striatalinferiorfrontalcortex_md","inferiorfrontalsuperiorfrontalcortex_md","fornix_exfimbria_md","fornix_wmv","cingulatecingulum_wmv","parahippocampalcingulum_wmv","corticospinalpyramidal_wmv","anteriorthalamicradiations_wmv","uncinate_wmv","inferiorlongitudinalfasiculus_wmv","inferiorfrontooccipitalfasiculus_wmv","forcepsmajor_wmv","forcepsminor_wmv","corpuscallosum_wmv","superiorlongitudinalfasiculus_wmv","temporalsuperiorlongitudinalfasiculus_wmv","parietalsuperiorlongitudinalfasiculus_wmv","superiorcorticostriate_wmv","superiorcorticostriatefrontalcortex_wmv","superiorcorticostriateparietalcortex_wmv","striatalinferiorfrontalcortex_wmv","inferiorfrontalsuperiorfrontalcortex_wmv","fornix_exfimbria_wmv"),
                             coef=rep("NA",times=162))

ROIs_coef_gmwm[new_coef_gmwm[1,1],2] <- new_coef_gmwm[1,2]

for (i in 1:nrow(new_coef_gmwm)) {
  ROIs_coef_gmwm[new_coef_gmwm[i,1],2] <- new_coef_gmwm[i,2]
}
ROI_gmwm<-c()
for (i in 1:162) {
  if (ROIs_coef_gmwm[i,2] != "NA") ROI_gmwm[i]<- ROIs_coef_gmwm[i,1]
}
ROI_gmwm<-na.omit(ROI_gmwm)
ROI_gmwm_plus<-c(ROI_gmwm[1])     #Compute the regularized regions in a form easily copy/paste to the model
for (i in 2:length(ROI_gmwm)) {
  ROI_gmwm_plus<-paste(ROI_gmwm_plus,"+",ROI_gmwm[i])
}
cat(ROI_gmwm_plus)
  
#Regularized white matter metrics model 
Model_GMWMlat_15_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ caudalmiddlefrontal_ct + fusiform_ct + lateralorbitofrontal_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsorbitalis_ct + pericalcarine_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiortemporal_ct + insula_ct + bankssts_sa + cuneus_sa + entorhinal_sa + lateraloccipital_sa + parsorbitalis_sa + postcentral_sa + precentral_sa + superiortemporal_sa + fusiform_gmv + inferiorparietal_gmv + parahippocampal_gmv + paracentral_gmv + pericalcarine_gmv + precentral_gmv + precuneus_gmv + rostralmiddlefrontal_gmv + supramarginal_gmv + cingulatecingulum_fa + parahippocampalcingulum_fa + uncinate_fa + inferiorfrontooccipitalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + superiorcorticostriateparietalcortex_fa + fornix_exfimbria_fa + parahippocampalcingulum_md + inferiorlongitudinalfasiculus_md + forcepsmajor_md + temporalsuperiorlongitudinalfasiculus_md + striatalinferiorfrontalcortex_md + anteriorthalamicradiations_wmv + inferiorlongitudinalfasiculus_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + superiorcorticostriateparietalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_GMWMlat_15_reg_free <- sem(Model_GMWMlat_15_reg_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMWMlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized grey & white matter metrics model
Model_GMWMlat_15_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*caudalmiddlefrontal_ct + a*fusiform_ct + a*lateralorbitofrontal_ct + a*medialorbitofrontal_ct + a*middletemporal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsorbitalis_ct + a*pericalcarine_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiortemporal_ct + a*insula_ct + a*bankssts_sa + a*cuneus_sa + a*entorhinal_sa + a*lateraloccipital_sa + a*parsorbitalis_sa + a*postcentral_sa + a*precentral_sa + a*superiortemporal_sa + a*fusiform_gmv + a*inferiorparietal_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*pericalcarine_gmv + a*precentral_gmv + a*precuneus_gmv + a*rostralmiddlefrontal_gmv + a*supramarginal_gmv + a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*uncinate_fa + a*inferiorfrontooccipitalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*superiorcorticostriateparietalcortex_fa + a*fornix_exfimbria_fa + a*parahippocampalcingulum_md + a*inferiorlongitudinalfasiculus_md + a*forcepsmajor_md + a*temporalsuperiorlongitudinalfasiculus_md + a*striatalinferiorfrontalcortex_md + a*anteriorthalamicradiations_wmv + a*inferiorlongitudinalfasiculus_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*superiorcorticostriateparietalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_GMWMlat_15_reg_constrained <- sem(Model_GMWMlat_15_reg_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMWMlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_GMWMlat_15_reg_free,fit_sem_GMWMlat_15_reg_constrained)

  
#Model for grey & white matter metrics with all the regions of interest (162 ROIs) 
Model_GMWMlat_15_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ bankssts_ct + caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + cuneus_ct + entorhinal_ct + fusiform_ct + inferiorparietal_ct + inferiortemporal_ct + isthmuscingulate_ct + lateraloccipital_ct + lateralorbitofrontal_ct + lingual_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsopercularis_ct + parsorbitalis_ct + parstriangularis_ct + pericalcarine_ct + postcentral_ct + posteriorcingulate_ct + precentral_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiorparietal_ct + superiortemporal_ct + supramarginal_ct + frontalpole_ct + temporalpole_ct + transversetemporal_ct + insula_ct + 
                                                 bankssts_sa + caudalanteriorcingulate_sa + caudalmiddlefrontal_sa + cuneus_sa + entorhinal_sa + fusiform_sa + inferiorparietal_sa + inferiortemporal_sa + isthmuscingulate_sa + lateraloccipital_sa + lateralorbitofrontal_sa + lingual_sa + medialorbitofrontal_sa + middletemporal_sa + parahippocampal_sa + paracentral_sa + parsopercularis_sa + parsorbitalis_sa + parstriangularis_sa + pericalcarine_sa + postcentral_sa + posteriorcingulate_sa + precentral_sa + precuneus_sa + rostralanteriorcingulate_sa + rostralmiddlefrontal_sa + superiorfrontal_sa + superiorparietal_sa + superiortemporal_sa + supramarginal_sa + frontalpole_sa + temporalpole_sa + transversetemporal_sa + insula_sa + 
                                                 bankssts_gmv + caudalanteriorcingulate_gmv + caudalmiddlefrontal_gmv + cuneus_gmv + entorhinal_gmv + fusiform_gmv + inferiorparietal_gmv + inferiortemporal_gmv + isthmuscingulate_gmv + lateraloccipital_gmv + lateralorbitofrontal_gmv + lingual_gmv + medialorbitofrontal_gmv + middletemporal_gmv + parahippocampal_gmv + paracentral_gmv + parsopercularis_gmv + parsorbitalis_gmv + parstriangularis_gmv + pericalcarine_gmv + postcentral_gmv + posteriorcingulate_gmv + precentral_gmv + precuneus_gmv + rostralanteriorcingulate_gmv + rostralmiddlefrontal_gmv + superiorfrontal_gmv + superiorparietal_gmv + superiortemporal_gmv + supramarginal_gmv + frontalpole_gmv + temporalpole_gmv + transversetemporal_gmv + insula_gmv + 
                                                 fornix_fa + cingulatecingulum_fa + parahippocampalcingulum_fa + corticospinalpyramidal_fa + anteriorthalamicradiations_fa + uncinate_fa + inferiorlongitudinalfasiculus_fa + inferiorfrontooccipitalfasiculus_fa + forcepsmajor_fa + forcepsminor_fa + corpuscallosum_fa + superiorlongitudinalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + parietalsuperiorlongitudinalfasiculus_fa + superiorcorticostriate_fa + superiorcorticostriatefrontalcortex_fa + superiorcorticostriateparietalcortex_fa + striatalinferiorfrontalcortex_fa + inferiorfrontalsuperiorfrontalcortex_fa + fornix_exfimbria_fa + 
                                                 fornix_md + cingulatecingulum_md + parahippocampalcingulum_md + corticospinalpyramidal_md + anteriorthalamicradiations_md + uncinate_md + inferiorlongitudinalfasiculus_md + inferiorfrontooccipitalfasiculus_md + forcepsmajor_md + forcepsminor_md + corpuscallosum_md + superiorlongitudinalfasiculus_md + temporalsuperiorlongitudinalfasiculus_md + parietalsuperiorlongitudinalfasiculus_md + superiorcorticostriate_md + superiorcorticostriatefrontalcortex_md + superiorcorticostriateparietalcortex_md + striatalinferiorfrontalcortex_md + inferiorfrontalsuperiorfrontalcortex_md + fornix_exfimbria_md + 
                                                 fornix_wmv + cingulatecingulum_wmv + parahippocampalcingulum_wmv + corticospinalpyramidal_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + forcepsmajor_wmv + forcepsminor_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + temporalsuperiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + superiorcorticostriate_wmv + superiorcorticostriatefrontalcortex_wmv + superiorcorticostriateparietalcortex_wmv + striatalinferiorfrontalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv + fornix_exfimbria_wmv'

fit_sem_GMWMlat_15_all_free <- sem(Model_GMWMlat_15_all_free, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMWMlat_15_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the grey & white matter metrics model with all the regions
Model_GMWMlat_15_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*bankssts_ct + a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*cuneus_ct + a*entorhinal_ct + a*fusiform_ct + a*inferiorparietal_ct + a*inferiortemporal_ct + a*isthmuscingulate_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*lingual_ct + a*medialorbitofrontal_ct + a*middletemporal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsopercularis_ct + a*parsorbitalis_ct + a*parstriangularis_ct + a*pericalcarine_ct + a*postcentral_ct + a*posteriorcingulate_ct + a*precentral_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*superiortemporal_ct + a*supramarginal_ct + a*frontalpole_ct + a*temporalpole_ct + a*transversetemporal_ct + a*insula_ct + 
                                                        a*bankssts_sa + a*caudalanteriorcingulate_sa + a*caudalmiddlefrontal_sa + a*cuneus_sa + a*entorhinal_sa + a*fusiform_sa + a*inferiorparietal_sa + a*inferiortemporal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*lateralorbitofrontal_sa + a*lingual_sa + a*medialorbitofrontal_sa + a*middletemporal_sa + a*parahippocampal_sa + a*paracentral_sa + a*parsopercularis_sa + a*parsorbitalis_sa + a*parstriangularis_sa + a*pericalcarine_sa + a*postcentral_sa + a*posteriorcingulate_sa + a*precentral_sa + a*precuneus_sa + a*rostralanteriorcingulate_sa + a*rostralmiddlefrontal_sa + a*superiorfrontal_sa + a*superiorparietal_sa + a*superiortemporal_sa + a*supramarginal_sa + a*frontalpole_sa + a*temporalpole_sa + a*transversetemporal_sa + a*insula_sa + 
                                                        a*bankssts_gmv + a*caudalanteriorcingulate_gmv + a*caudalmiddlefrontal_gmv + a*cuneus_gmv + a*entorhinal_gmv + a*fusiform_gmv + a*inferiorparietal_gmv + a*inferiortemporal_gmv + a*isthmuscingulate_gmv + a*lateraloccipital_gmv + a*lateralorbitofrontal_gmv + a*lingual_gmv + a*medialorbitofrontal_gmv + a*middletemporal_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*parsopercularis_gmv + a*parsorbitalis_gmv + a*parstriangularis_gmv + a*pericalcarine_gmv + a*postcentral_gmv + a*posteriorcingulate_gmv + a*precentral_gmv + a*precuneus_gmv + a*rostralanteriorcingulate_gmv + a*rostralmiddlefrontal_gmv + a*superiorfrontal_gmv + a*superiorparietal_gmv + a*superiortemporal_gmv + a*supramarginal_gmv + a*frontalpole_gmv + a*temporalpole_gmv + a*transversetemporal_gmv + a*insula_gmv + 
                                                        a*fornix_fa + a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*corticospinalpyramidal_fa + a*anteriorthalamicradiations_fa + a*uncinate_fa + a*inferiorlongitudinalfasiculus_fa + a*inferiorfrontooccipitalfasiculus_fa + a*forcepsmajor_fa + a*forcepsminor_fa + a*corpuscallosum_fa + a*superiorlongitudinalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*parietalsuperiorlongitudinalfasiculus_fa + a*superiorcorticostriate_fa + a*superiorcorticostriatefrontalcortex_fa + a*superiorcorticostriateparietalcortex_fa + a*striatalinferiorfrontalcortex_fa + a*inferiorfrontalsuperiorfrontalcortex_fa + a*fornix_exfimbria_fa + 
                                                        a*fornix_md + a*cingulatecingulum_md + a*parahippocampalcingulum_md + a*corticospinalpyramidal_md + a*anteriorthalamicradiations_md + a*uncinate_md + a*inferiorlongitudinalfasiculus_md + a*inferiorfrontooccipitalfasiculus_md + a*forcepsmajor_md + a*forcepsminor_md + a*corpuscallosum_md + a*superiorlongitudinalfasiculus_md + a*temporalsuperiorlongitudinalfasiculus_md + a*parietalsuperiorlongitudinalfasiculus_md + a*superiorcorticostriate_md + a*superiorcorticostriatefrontalcortex_md + a*superiorcorticostriateparietalcortex_md + a*striatalinferiorfrontalcortex_md + a*inferiorfrontalsuperiorfrontalcortex_md + a*fornix_exfimbria_md + 
                                                        a*fornix_wmv + a*cingulatecingulum_wmv + a*parahippocampalcingulum_wmv + a*corticospinalpyramidal_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*forcepsmajor_wmv + a*forcepsminor_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*temporalsuperiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*superiorcorticostriate_wmv + a*superiorcorticostriatefrontalcortex_wmv + a*superiorcorticostriateparietalcortex_wmv + a*striatalinferiorfrontalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv + a*fornix_exfimbria_wmv'

fit_sem_GMWMlat_15_all_constrained <- sem(Model_GMWMlat_15_all_constrained, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMWMlat_15_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_GMWMlat_15_all_free,fit_sem_GMWMlat_15_all_constrained)

  
#-------------------------------------------------------------------------------------------------------------------------------
#Sample 15%: Compare the regularized models with grey and white matter
#-------------------------------------------------------------------------------------------------------------------------------
  
fit_models_15<-data.frame(label=c("Model Grey & White Matter","Model Grey Matter","Model White Matter"),
                         AIC=rep(NA,times=3),
                         BIC=rep(NA,times=3))

##Grey & White Matter
summary(fit_sem_GMWMlat_15_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)
fit_models_15[1,2]<-AIC(fit_sem_GMWMlat_15_reg_free)
fit_models_15[1,3]<-BIC(fit_sem_GMWMlat_15_reg_free)

#'I don't need the white matter if I already have the grey matter'
Model_GMWM_WMconstrainedlat_15 <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ caudalmiddlefrontal_ct + fusiform_ct + lateralorbitofrontal_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsorbitalis_ct + pericalcarine_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiortemporal_ct + insula_ct + bankssts_sa + cuneus_sa + entorhinal_sa + lateraloccipital_sa + parsorbitalis_sa + postcentral_sa + precentral_sa + superiortemporal_sa + fusiform_gmv + inferiorparietal_gmv + parahippocampal_gmv + paracentral_gmv + pericalcarine_gmv + precentral_gmv + precuneus_gmv + rostralmiddlefrontal_gmv + supramarginal_gmv + 
                                                      0*cingulatecingulum_fa + 0*parahippocampalcingulum_fa + 0*uncinate_fa + 0*inferiorfrontooccipitalfasiculus_fa + 0*temporalsuperiorlongitudinalfasiculus_fa + 0*superiorcorticostriateparietalcortex_fa + 0*fornix_exfimbria_fa + 0*parahippocampalcingulum_md + 0*inferiorlongitudinalfasiculus_md + 0*forcepsmajor_md + 0*temporalsuperiorlongitudinalfasiculus_md + 0*striatalinferiorfrontalcortex_md + 0*anteriorthalamicradiations_wmv + 0*inferiorlongitudinalfasiculus_wmv + 0*corpuscallosum_wmv + 0*superiorlongitudinalfasiculus_wmv + 0*superiorcorticostriateparietalcortex_wmv + 0*inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_GMWM_WMconstrainedlat_15 <- sem(Model_GMWM_WMconstrainedlat_15, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMWM_WMconstrainedlat_15, fit.measures=TRUE,rsquare=T,standardized=T)
fit_models_T015[2,2]<-AIC(fit_sem_GMWM_WMconstrainedlat_15)
fit_models_T015[2,3]<-BIC(fit_sem_GMWM_WMconstrainedlat_15)

anova(fit_sem_GMWMlat_15_reg_free,fit_sem_GMWM_WMconstrainedlat_15)

#'I don't need the grey matter if I already have the white matter'
Model_GMWM_GMconstrainedlat_15 <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ 0*caudalmiddlefrontal_ct + 0*fusiform_ct + 0*lateralorbitofrontal_ct + 0*medialorbitofrontal_ct + 0*middletemporal_ct + 0*parahippocampal_ct + 0*paracentral_ct + 0*parsorbitalis_ct + 0*pericalcarine_ct + 0*precuneus_ct + 0*rostralanteriorcingulate_ct + 0*rostralmiddlefrontal_ct + 0*superiorfrontal_ct + 0*superiortemporal_ct + 0*insula_ct + 0*bankssts_sa + 0*cuneus_sa + 0*entorhinal_sa + 0*lateraloccipital_sa + 0*parsorbitalis_sa + 0*postcentral_sa + 0*precentral_sa + 0*superiortemporal_sa + 0*fusiform_gmv + 0*inferiorparietal_gmv + 0*parahippocampal_gmv + 0*paracentral_gmv + 0*pericalcarine_gmv + 0*precentral_gmv + 0*precuneus_gmv + 0*rostralmiddlefrontal_gmv + 0*supramarginal_gmv +
                                                        cingulatecingulum_fa + parahippocampalcingulum_fa + uncinate_fa + inferiorfrontooccipitalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + superiorcorticostriateparietalcortex_fa + fornix_exfimbria_fa + parahippocampalcingulum_md + inferiorlongitudinalfasiculus_md + forcepsmajor_md + temporalsuperiorlongitudinalfasiculus_md + striatalinferiorfrontalcortex_md + anteriorthalamicradiations_wmv + inferiorlongitudinalfasiculus_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + superiorcorticostriateparietalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv'
                                  
fit_sem_GMWM_GMconstrainedlat_15 <- sem(Model_GMWM_GMconstrainedlat_15, data=data_total_subsample15,estimator="mlr",missing="fiml")
summary(fit_sem_GMWM_GMconstrainedlat_15, fit.measures=TRUE,rsquare=T,standardized=T)
fit_models_15[3,2]<-AIC(fit_sem_GMWM_GMconstrainedlat_15)
fit_models_15[3,3]<-BIC(fit_sem_GMWM_GMconstrainedlat_15)

anova(fit_sem_GMWMlat_15_reg_free,fit_sem_GMWM_GMconstrainedlat_15)

#'Grey and white matter have equally sized, complementary effects'
summary(fit_sem_GMWMlat_15_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_GMWMlat_15_reg_free,fit_sem_GMWMlat_15_reg_constrained)

#Data Visualization for the comparison between the grey matter model, the white matter model and the model with grey & white matter
#Create a data.frame with model names
AIC_15 <- data.frame(Model = fit_models_15[,1],
                     AIC = factor(rep('AIC', each = 3)),
                     values = fit_models_15[,2])

BIC_15 <- data.frame(Model = fit_models_15[,1],
                     BIC = factor(rep('BIC', each = 3)),
                     values = fit_models_15[,3])

#Plot the 'raw' AIC's
plot_fitraw_aic_15<-ggplot(AIC_15, aes(values,AIC, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1),colour="black") +
  coord_cartesian(xlim = c(58500, 58665)) +
  xlab('Information criterion (AIC)') +
  ylab('') +
  theme_grey(base_size = 30) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
        axis.line = element_line(colour = "black"))+  # adding a black line for x and y axis
  geom_text(aes(label=c("Grey & White\nMatter Model","Grey Matter Model","White Matter Model")), size = 6, color="white",position = position_dodge(width = 1),hjust = 1.1) +
  theme(legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  scale_fill_manual("legend", values = c("Model Grey & White Matter" = "black", "Model Grey Matter" = "coral1", "Model White Matter" = "darkcyan"))
ggsave("plot_fitraw_aic15.png")

#Plot the 'raw' BIC's
plot_fitraw_bic_15<-ggplot(BIC_15, aes(values, BIC, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1),colour="black") +
  coord_cartesian(ylim = c(58840, 59000)) +
  xlab('Information criterion (BIC)') +
  ylab('') +
  theme_grey(base_size = 30) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
        axis.line = element_line(colour = "black"))+  # adding a black line for x and y axis
  geom_text(aes(label=c("Grey & White\nMatter Model","Grey Matter Model","White Matter Model")), size = 6, color="white",position = position_dodge(width = 1),hjust = 1.1) +
  theme(legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  scale_fill_manual("legend", values = c("Model Grey & White Matter" = "black", "Model Grey Matter" = "coral1", "Model White Matter" = "darkcyan"))
ggsave("plot_fitraw_bic15.png")  




#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

  
  
#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Individual models grey matter
#-------------------------------------------------------------------------------------------------------------------------------
  
#INDIVIDUAL MODELS GREY MATTER  
#Create a table to save the regression coefficients or the fit of the ROI model for each ROIs
data_subsample_GMlat_85_GGSEG <- data.frame(label=c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus","entorhinal","fusiform","inferiorparietal","inferiortemporal","isthmuscingulate","lateraloccipital","lateralorbitofrontal","lingual","medialorbitofrontal","middletemporal","parahippocampal","paracentral","parsopercularis","parsorbitalis","parstriangularis","pericalcarine","postcentral","posteriorcingulate","precentral","precuneus","rostralanteriorcingulate","rostralmiddlefrontal","superiorfrontal","superiorparietal","superiortemporal","supramarginal","frontalpole","temporalpole","transversetemporal","insula"),
                                              CT_fit=rep(NA, times=34),
                                              CT_p=rep(NA, times=34),
                                              CT_std=rep(NA, times=34),
                                              SA_fit=rep(NA, times=34),
                                              SA_p=rep(NA, times=34),
                                              SA_std=rep(NA, times=34),
                                              GMV_fit=rep(NA, times=34),
                                              GMV_p=rep(NA, times=34),
                                              GMV_std=rep(NA, times=34))

#Create tables to save the parameters of the models for each ROI in each metric
paramest_GMlat_CT_85 <- data.frame(lhs=rep(NA, times=34),
                                     op=rep(NA, times=34),
                                     rhs=rep(NA, times=34),
                                     est=rep(NA, times=34),
                                     se=rep(NA, times=34),
                                     z=rep(NA, times=34),
                                     pvalue=rep(NA, times=34),
                                     ci.lower=rep(NA, times=34),
                                     ci.upper=rep(NA, times=34),
                                     std.lv=rep(NA, times=34),
                                     std.all=rep(NA, times=34),
                                     std.nox=rep(NA, times=34))
paramest_GMlat_SA_85 <- data.frame(lhs=rep(NA, times=34),
                                     op=rep(NA, times=34),
                                     rhs=rep(NA, times=34),
                                     est=rep(NA, times=34),
                                     se=rep(NA, times=34),
                                     z=rep(NA, times=34),
                                     pvalue=rep(NA, times=34),
                                     ci.lower=rep(NA, times=34),
                                     ci.upper=rep(NA, times=34),
                                     std.lv=rep(NA, times=34),
                                     std.all=rep(NA, times=34),
                                     std.nox=rep(NA, times=34))
paramest_GMlat_GMV_85 <- data.frame(lhs=rep(NA, times=34),
                                     op=rep(NA, times=34),
                                     rhs=rep(NA, times=34),
                                     est=rep(NA, times=34),
                                     se=rep(NA, times=34),
                                     z=rep(NA, times=34),
                                     pvalue=rep(NA, times=34),
                                     ci.lower=rep(NA, times=34),
                                     ci.upper=rep(NA, times=34),
                                     std.lv=rep(NA, times=34),
                                     std.all=rep(NA, times=34),
                                     std.nox=rep(NA, times=34))

#Cortical Thickness (34 regions)
##bankssts
Model_CT_bankssts <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                      Cognitive_factor ~ bankssts_ct'
fit_sem_CT_bankssts <- sem(Model_CT_bankssts, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_bankssts, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[1,2]<-BIC(fit_sem_CT_bankssts)
data_subsample_GMlat_85_GGSEG[1,3]<-parameterEstimates(fit_sem_CT_bankssts)[6,7]
data_subsample_GMlat_85_GGSEG[1,4]<-lavInspect(fit_sem_CT_bankssts,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[1,] <-parameterEstimates(fit_sem_CT_bankssts,standardized = TRUE, rsquare=TRUE)[6,]

##caudalanteriorcingulate
Model_CT_caudalanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ caudalanteriorcingulate_ct'
fit_sem_CT_caudalanteriorcingulate <- sem(Model_CT_caudalanteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_caudalanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[2,2]<-BIC(fit_sem_CT_caudalanteriorcingulate)
data_subsample_GMlat_85_GGSEG[2,3]<-parameterEstimates(fit_sem_CT_caudalanteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[2,4]<-lavInspect(fit_sem_CT_caudalanteriorcingulate,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[2,] <-parameterEstimates(fit_sem_CT_caudalanteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##caudalmiddlefrontal
Model_CT_caudalmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ caudalmiddlefrontal_ct'
fit_sem_CT_caudalmiddlefrontal <- sem(Model_CT_caudalmiddlefrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_caudalmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[3,2]<-BIC(fit_sem_CT_caudalmiddlefrontal)
data_subsample_GMlat_85_GGSEG[3,3]<-parameterEstimates(fit_sem_CT_caudalmiddlefrontal)[6,7]
data_subsample_GMlat_85_GGSEG[3,4]<-lavInspect(fit_sem_CT_caudalmiddlefrontal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[3,] <-parameterEstimates(fit_sem_CT_caudalmiddlefrontal,standardized = TRUE, rsquare=TRUE)[6,]

##cuneus
Model_CT_cuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ cuneus_ct'
fit_sem_CT_cuneus <- sem(Model_CT_cuneus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_cuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[4,2]<-BIC(fit_sem_CT_cuneus)
data_subsample_GMlat_85_GGSEG[4,3]<-parameterEstimates(fit_sem_CT_cuneus)[6,7]
data_subsample_GMlat_85_GGSEG[4,4]<-lavInspect(fit_sem_CT_cuneus,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[4,] <-parameterEstimates(fit_sem_CT_cuneus,standardized = TRUE, rsquare=TRUE)[6,]

##entorhinal
Model_CT_entorhinal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ entorhinal_ct'
fit_sem_CT_entorhinal <- sem(Model_CT_entorhinal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_entorhinal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[5,2]<-BIC(fit_sem_CT_entorhinal)
data_subsample_GMlat_85_GGSEG[5,3]<-parameterEstimates(fit_sem_CT_entorhinal)[6,7]
data_subsample_GMlat_85_GGSEG[5,4]<-lavInspect(fit_sem_CT_entorhinal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[5,] <-parameterEstimates(fit_sem_CT_entorhinal,standardized = TRUE, rsquare=TRUE)[6,]

##fusiform
Model_CT_fusiform <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fusiform_ct'
fit_sem_CT_fusiform <- sem(Model_CT_fusiform, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_fusiform, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[6,2]<-BIC(fit_sem_CT_fusiform)
data_subsample_GMlat_85_GGSEG[6,3]<-parameterEstimates(fit_sem_CT_fusiform)[6,7]
data_subsample_GMlat_85_GGSEG[6,4]<-lavInspect(fit_sem_CT_fusiform,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[6,] <-parameterEstimates(fit_sem_CT_fusiform,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorparietal
Model_CT_inferiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorparietal_ct'
fit_sem_CT_inferiorparietal <- sem(Model_CT_inferiorparietal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_inferiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[7,2]<-BIC(fit_sem_CT_inferiorparietal)
data_subsample_GMlat_85_GGSEG[7,3]<-parameterEstimates(fit_sem_CT_inferiorparietal)[6,7]
data_subsample_GMlat_85_GGSEG[7,4]<-lavInspect(fit_sem_CT_inferiorparietal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[7,] <-parameterEstimates(fit_sem_CT_inferiorparietal,standardized = TRUE, rsquare=TRUE)[6,]

##inferiortemporal
Model_CT_inferiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiortemporal_ct'
fit_sem_CT_inferiortemporal <- sem(Model_CT_inferiortemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_inferiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[8,2]<-BIC(fit_sem_CT_inferiortemporal)
data_subsample_GMlat_85_GGSEG[8,3]<-parameterEstimates(fit_sem_CT_inferiortemporal)[6,7]
data_subsample_GMlat_85_GGSEG[8,4]<-lavInspect(fit_sem_CT_inferiortemporal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[8,] <-parameterEstimates(fit_sem_CT_inferiortemporal,standardized = TRUE, rsquare=TRUE)[6,]

##isthmuscingulate
Model_CT_isthmuscingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ isthmuscingulate_ct'
fit_sem_CT_isthmuscingulate <- sem(Model_CT_isthmuscingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_isthmuscingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[9,2]<-BIC(fit_sem_CT_isthmuscingulate)
data_subsample_GMlat_85_GGSEG[9,3]<-parameterEstimates(fit_sem_CT_isthmuscingulate)[6,7]
data_subsample_GMlat_85_GGSEG[9,4]<-lavInspect(fit_sem_CT_isthmuscingulate,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[9,] <-parameterEstimates(fit_sem_CT_isthmuscingulate,standardized = TRUE, rsquare=TRUE)[6,]

##lateraloccipital
Model_CT_lateraloccipital <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lateraloccipital_ct'
fit_sem_CT_lateraloccipital <- sem(Model_CT_lateraloccipital, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_lateraloccipital, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[10,2]<-BIC(fit_sem_CT_lateraloccipital)
data_subsample_GMlat_85_GGSEG[10,3]<-parameterEstimates(fit_sem_CT_lateraloccipital)[6,7]
data_subsample_GMlat_85_GGSEG[10,4]<-lavInspect(fit_sem_CT_lateraloccipital,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[10,] <-parameterEstimates(fit_sem_CT_lateraloccipital,standardized = TRUE, rsquare=TRUE)[6,]

##lateralorbitofrontal
Model_CT_lateralorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lateralorbitofrontal_ct'
fit_sem_CT_lateralorbitofrontal <- sem(Model_CT_lateralorbitofrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_lateralorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[11,2]<-BIC(fit_sem_CT_lateralorbitofrontal)
data_subsample_GMlat_85_GGSEG[11,3]<-parameterEstimates(fit_sem_CT_lateralorbitofrontal)[6,7]
data_subsample_GMlat_85_GGSEG[11,4]<-lavInspect(fit_sem_CT_lateralorbitofrontal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[11,] <-parameterEstimates(fit_sem_CT_lateralorbitofrontal,standardized = TRUE, rsquare=TRUE)[6,]

##lingual
Model_CT_lingual <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lingual_ct'
fit_sem_CT_lingual <- sem(Model_CT_lingual, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_lingual, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[12,2]<-BIC(fit_sem_CT_lingual)
data_subsample_GMlat_85_GGSEG[12,3]<-parameterEstimates(fit_sem_CT_lingual)[6,7]
data_subsample_GMlat_85_GGSEG[12,4]<-lavInspect(fit_sem_CT_lingual,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[12,] <-parameterEstimates(fit_sem_CT_lingual,standardized = TRUE, rsquare=TRUE)[6,]

##medialorbitofrontal
Model_CT_medialorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ medialorbitofrontal_ct'
fit_sem_CT_medialorbitofrontal <- sem(Model_CT_medialorbitofrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_medialorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[13,2]<-BIC(fit_sem_CT_medialorbitofrontal)
data_subsample_GMlat_85_GGSEG[13,3]<-parameterEstimates(fit_sem_CT_medialorbitofrontal)[6,7]
data_subsample_GMlat_85_GGSEG[13,4]<-lavInspect(fit_sem_CT_medialorbitofrontal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[13,] <-parameterEstimates(fit_sem_CT_medialorbitofrontal,standardized = TRUE, rsquare=TRUE)[6,]

##middletemporal
Model_CT_middletemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ middletemporal_ct'
fit_sem_CT_middletemporal <- sem(Model_CT_middletemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_middletemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[14,2]<-BIC(fit_sem_CT_middletemporal)
data_subsample_GMlat_85_GGSEG[14,3]<-parameterEstimates(fit_sem_CT_middletemporal)[6,7]
data_subsample_GMlat_85_GGSEG[14,4]<-lavInspect(fit_sem_CT_middletemporal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[14,] <-parameterEstimates(fit_sem_CT_middletemporal,standardized = TRUE, rsquare=TRUE)[6,]

##parahippocampal
Model_CT_parahippocampal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parahippocampal_ct'
fit_sem_CT_parahippocampal <- sem(Model_CT_parahippocampal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_parahippocampal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[15,2]<-BIC(fit_sem_CT_parahippocampal)
data_subsample_GMlat_85_GGSEG[15,3]<-parameterEstimates(fit_sem_CT_parahippocampal)[6,7]
data_subsample_GMlat_85_GGSEG[15,4]<-lavInspect(fit_sem_CT_parahippocampal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[15,] <-parameterEstimates(fit_sem_CT_parahippocampal,standardized = TRUE, rsquare=TRUE)[6,]

##paracentral
Model_CT_paracentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ paracentral_ct'
fit_sem_CT_paracentral <- sem(Model_CT_paracentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_paracentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[16,2]<-BIC(fit_sem_CT_paracentral)
data_subsample_GMlat_85_GGSEG[16,3]<-parameterEstimates(fit_sem_CT_paracentral)[6,7]
data_subsample_GMlat_85_GGSEG[16,4]<-lavInspect(fit_sem_CT_paracentral,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[16,] <-parameterEstimates(fit_sem_CT_paracentral,standardized = TRUE, rsquare=TRUE)[6,]

##parsopercularis
Model_CT_parsopercularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parsopercularis_ct'
fit_sem_CT_parsopercularis <- sem(Model_CT_parsopercularis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_parsopercularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[17,2]<-BIC(fit_sem_CT_parsopercularis)
data_subsample_GMlat_85_GGSEG[17,3]<-parameterEstimates(fit_sem_CT_parsopercularis)[6,7]
data_subsample_GMlat_85_GGSEG[17,4]<-lavInspect(fit_sem_CT_parsopercularis,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[17,] <-parameterEstimates(fit_sem_CT_parsopercularis,standardized = TRUE, rsquare=TRUE)[6,]

##parsorbitalis
Model_CT_parsorbitalis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parsorbitalis_ct'
fit_sem_CT_parsorbitalis <- sem(Model_CT_parsorbitalis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_parsorbitalis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[18,2]<-BIC(fit_sem_CT_parsorbitalis)
data_subsample_GMlat_85_GGSEG[18,3]<-parameterEstimates(fit_sem_CT_parsorbitalis)[6,7]
data_subsample_GMlat_85_GGSEG[18,4]<-lavInspect(fit_sem_CT_parsorbitalis,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[18,] <-parameterEstimates(fit_sem_CT_parsorbitalis,standardized = TRUE, rsquare=TRUE)[6,]

##parstriangularis
Model_CT_parstriangularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parstriangularis_ct'
fit_sem_CT_parstriangularis <- sem(Model_CT_parstriangularis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_parstriangularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[19,2]<-BIC(fit_sem_CT_parstriangularis)
data_subsample_GMlat_85_GGSEG[19,3]<-parameterEstimates(fit_sem_CT_parstriangularis)[6,7]
data_subsample_GMlat_85_GGSEG[19,4]<-lavInspect(fit_sem_CT_parstriangularis,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[19,] <-parameterEstimates(fit_sem_CT_parstriangularis,standardized = TRUE, rsquare=TRUE)[6,]

##pericalcarine
Model_CT_pericalcarine <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ pericalcarine_ct'
fit_sem_CT_pericalcarine <- sem(Model_CT_pericalcarine, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_pericalcarine, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[20,2]<-BIC(fit_sem_CT_pericalcarine)
data_subsample_GMlat_85_GGSEG[20,3]<-parameterEstimates(fit_sem_CT_pericalcarine)[6,7]
data_subsample_GMlat_85_GGSEG[20,4]<-lavInspect(fit_sem_CT_pericalcarine,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[20,] <-parameterEstimates(fit_sem_CT_pericalcarine,standardized = TRUE, rsquare=TRUE)[6,]

##postcentral
Model_CT_postcentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ postcentral_ct'
fit_sem_CT_postcentral <- sem(Model_CT_postcentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_postcentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[21,2]<-BIC(fit_sem_CT_postcentral)
data_subsample_GMlat_85_GGSEG[21,3]<-parameterEstimates(fit_sem_CT_postcentral)[6,7]
data_subsample_GMlat_85_GGSEG[21,4]<-lavInspect(fit_sem_CT_postcentral,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[21,] <-parameterEstimates(fit_sem_CT_postcentral,standardized = TRUE, rsquare=TRUE)[6,]

##posteriorcingulate
Model_CT_posteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ posteriorcingulate_ct'
fit_sem_CT_posteriorcingulate <- sem(Model_CT_posteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_posteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[22,2]<-BIC(fit_sem_CT_posteriorcingulate)
data_subsample_GMlat_85_GGSEG[22,3]<-parameterEstimates(fit_sem_CT_posteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[22,4]<-lavInspect(fit_sem_CT_posteriorcingulate,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[22,] <-parameterEstimates(fit_sem_CT_posteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##precentral
Model_CT_precentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ precentral_ct'
fit_sem_CT_precentral <- sem(Model_CT_precentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_precentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[23,2]<-BIC(fit_sem_CT_precentral)
data_subsample_GMlat_85_GGSEG[23,3]<-parameterEstimates(fit_sem_CT_precentral)[6,7]
data_subsample_GMlat_85_GGSEG[23,4]<-lavInspect(fit_sem_CT_precentral,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[23,] <-parameterEstimates(fit_sem_CT_precentral,standardized = TRUE, rsquare=TRUE)[6,]

##precuneus
Model_CT_precuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ precuneus_ct'
fit_sem_CT_precuneus <- sem(Model_CT_precuneus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_precuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[24,2]<-BIC(fit_sem_CT_precuneus)
data_subsample_GMlat_85_GGSEG[24,3]<-parameterEstimates(fit_sem_CT_precuneus)[6,7]
data_subsample_GMlat_85_GGSEG[24,4]<-lavInspect(fit_sem_CT_precuneus,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[24,] <-parameterEstimates(fit_sem_CT_precuneus,standardized = TRUE, rsquare=TRUE)[6,]

##rostralanteriorcingulate
Model_CT_rostralanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ rostralanteriorcingulate_ct'
fit_sem_CT_rostralanteriorcingulate <- sem(Model_CT_rostralanteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_rostralanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[25,2]<-BIC(fit_sem_CT_rostralanteriorcingulate)
data_subsample_GMlat_85_GGSEG[25,3]<-parameterEstimates(fit_sem_CT_rostralanteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[25,4]<-lavInspect(fit_sem_CT_rostralanteriorcingulate,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[25,] <-parameterEstimates(fit_sem_CT_rostralanteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##rostralmiddlefrontal
Model_CT_rostralmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ rostralmiddlefrontal_ct'
fit_sem_CT_rostralmiddlefrontal <- sem(Model_CT_rostralmiddlefrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_rostralmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[26,2]<-BIC(fit_sem_CT_rostralmiddlefrontal)
data_subsample_GMlat_85_GGSEG[26,3]<-parameterEstimates(fit_sem_CT_rostralmiddlefrontal)[6,7]
data_subsample_GMlat_85_GGSEG[26,4]<-lavInspect(fit_sem_CT_rostralmiddlefrontal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[26,] <-parameterEstimates(fit_sem_CT_rostralmiddlefrontal,standardized = TRUE, rsquare=TRUE)[6,]

##superiorfrontal
Model_CT_superiorfrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorfrontal_ct'
fit_sem_CT_superiorfrontal <- sem(Model_CT_superiorfrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_superiorfrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[27,2]<-BIC(fit_sem_CT_superiorfrontal)
data_subsample_GMlat_85_GGSEG[27,3]<-parameterEstimates(fit_sem_CT_superiorfrontal)[6,7]
data_subsample_GMlat_85_GGSEG[27,4]<-lavInspect(fit_sem_CT_superiorfrontal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[27,] <-parameterEstimates(fit_sem_CT_superiorfrontal,standardized = TRUE, rsquare=TRUE)[6,]

##superiorparietal
Model_CT_superiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorparietal_ct'
fit_sem_CT_superiorparietal <- sem(Model_CT_superiorparietal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_superiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[28,2]<-BIC(fit_sem_CT_superiorparietal)
data_subsample_GMlat_85_GGSEG[28,3]<-parameterEstimates(fit_sem_CT_superiorparietal)[6,7]
data_subsample_GMlat_85_GGSEG[28,4]<-lavInspect(fit_sem_CT_superiorparietal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[28,] <-parameterEstimates(fit_sem_CT_superiorparietal,standardized = TRUE, rsquare=TRUE)[6,]

##superiortemporal
Model_CT_superiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiortemporal_ct'
fit_sem_CT_superiortemporal <- sem(Model_CT_superiortemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_superiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[29,2]<-BIC(fit_sem_CT_superiortemporal)
data_subsample_GMlat_85_GGSEG[29,3]<-parameterEstimates(fit_sem_CT_superiortemporal)[6,7]
data_subsample_GMlat_85_GGSEG[29,4]<-lavInspect(fit_sem_CT_superiortemporal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[29,] <-parameterEstimates(fit_sem_CT_superiortemporal,standardized = TRUE, rsquare=TRUE)[6,]

##supramarginal
Model_CT_supramarginal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ supramarginal_ct'
fit_sem_CT_supramarginal <- sem(Model_CT_supramarginal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_supramarginal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[30,2]<-BIC(fit_sem_CT_supramarginal)
data_subsample_GMlat_85_GGSEG[30,3]<-parameterEstimates(fit_sem_CT_supramarginal)[6,7]
data_subsample_GMlat_85_GGSEG[30,4]<-lavInspect(fit_sem_CT_supramarginal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[30,] <-parameterEstimates(fit_sem_CT_supramarginal,standardized = TRUE, rsquare=TRUE)[6,]

##frontalpole
Model_CT_frontalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ frontalpole_ct'
fit_sem_CT_frontalpole <- sem(Model_CT_frontalpole, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_frontalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[31,2]<-BIC(fit_sem_CT_frontalpole)
data_subsample_GMlat_85_GGSEG[31,3]<-parameterEstimates(fit_sem_CT_frontalpole)[6,7]
data_subsample_GMlat_85_GGSEG[31,4]<-lavInspect(fit_sem_CT_frontalpole,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[31,] <-parameterEstimates(fit_sem_CT_frontalpole,standardized = TRUE, rsquare=TRUE)[6,]

##temporalpole
Model_CT_temporalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalpole_ct'
fit_sem_CT_temporalpole <- sem(Model_CT_temporalpole, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_temporalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[32,2]<-BIC(fit_sem_CT_temporalpole)
data_subsample_GMlat_85_GGSEG[32,3]<-parameterEstimates(fit_sem_CT_temporalpole)[6,7]
data_subsample_GMlat_85_GGSEG[32,4]<-lavInspect(fit_sem_CT_temporalpole,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[32,] <-parameterEstimates(fit_sem_CT_temporalpole,standardized = TRUE, rsquare=TRUE)[6,]

##transversetemporal
Model_CT_transversetemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ transversetemporal_ct'
fit_sem_CT_transversetemporal <- sem(Model_CT_transversetemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_transversetemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[33,2]<-BIC(fit_sem_CT_transversetemporal)
data_subsample_GMlat_85_GGSEG[33,3]<-parameterEstimates(fit_sem_CT_transversetemporal)[6,7]
data_subsample_GMlat_85_GGSEG[33,4]<-lavInspect(fit_sem_CT_transversetemporal,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[33,] <-parameterEstimates(fit_sem_CT_transversetemporal,standardized = TRUE, rsquare=TRUE)[6,]

##insula
Model_CT_insula <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ insula_ct'
fit_sem_CT_insula <- sem(Model_CT_insula, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CT_insula, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[34,2]<-BIC(fit_sem_CT_insula)
data_subsample_GMlat_85_GGSEG[34,3]<-parameterEstimates(fit_sem_CT_insula)[6,7]
data_subsample_GMlat_85_GGSEG[34,4]<-lavInspect(fit_sem_CT_insula,what = "std.all")$beta[1,2]
paramest_GMlat_CT_85[34,] <-parameterEstimates(fit_sem_CT_insula,standardized = TRUE, rsquare=TRUE)[6,]

write.csv(paramest_GMlat_CT_85, "paramest_GMlat_CT_85.csv", row.names=FALSE) 

#Surface Area (34 regions)
##bankssts
Model_SA_bankssts <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ bankssts_sa'
fit_sem_SA_bankssts <- sem(Model_SA_bankssts, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_bankssts, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[1,5]<-BIC(fit_sem_SA_bankssts)
data_subsample_GMlat_85_GGSEG[1,6]<-parameterEstimates(fit_sem_SA_bankssts)[6,7]
data_subsample_GMlat_85_GGSEG[1,7]<-lavInspect(fit_sem_SA_bankssts,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[1,] <-parameterEstimates(fit_sem_SA_bankssts,standardized = TRUE, rsquare=TRUE)[6,]

##caudalanteriorcingulate
Model_SA_caudalanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ caudalanteriorcingulate_sa'
fit_sem_SA_caudalanteriorcingulate <- sem(Model_SA_caudalanteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_caudalanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[2,5]<-BIC(fit_sem_SA_caudalanteriorcingulate)
data_subsample_GMlat_85_GGSEG[2,6]<-parameterEstimates(fit_sem_SA_caudalanteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[2,7]<-lavInspect(fit_sem_SA_caudalanteriorcingulate,what = "std.all")$beta[1,2]
paramtable_SA_caudalanteriorcingulate<-parameterEstimates(fit_sem_SA_caudalanteriorcingulate,standardized = TRUE, rsquare=TRUE)
paramest_GMlat_SA_85[2,] <-parameterEstimates(fit_sem_SA_caudalanteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##caudalmiddlefrontal
Model_SA_caudalmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ caudalmiddlefrontal_sa'
fit_sem_SA_caudalmiddlefrontal <- sem(Model_SA_caudalmiddlefrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_caudalmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[3,5]<-BIC(fit_sem_SA_caudalmiddlefrontal)
data_subsample_GMlat_85_GGSEG[3,6]<-parameterEstimates(fit_sem_SA_caudalmiddlefrontal)[6,7]
data_subsample_GMlat_85_GGSEG[3,7]<-lavInspect(fit_sem_SA_caudalmiddlefrontal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[3,] <-parameterEstimates(fit_sem_SA_caudalmiddlefrontal,standardized = TRUE, rsquare=TRUE)[6,]

##cuneus
Model_SA_cuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ cuneus_sa'
fit_sem_SA_cuneus <- sem(Model_SA_cuneus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_cuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[4,5]<-BIC(fit_sem_SA_cuneus)
data_subsample_GMlat_85_GGSEG[4,6]<-parameterEstimates(fit_sem_SA_cuneus)[6,7]
data_subsample_GMlat_85_GGSEG[4,7]<-lavInspect(fit_sem_SA_cuneus,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[4,] <-parameterEstimates(fit_sem_SA_cuneus,standardized = TRUE, rsquare=TRUE)[6,]

##entorhinal
Model_SA_entorhinal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ entorhinal_sa'
fit_sem_SA_entorhinal <- sem(Model_SA_entorhinal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_entorhinal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[5,5]<-BIC(fit_sem_SA_entorhinal)
data_subsample_GMlat_85_GGSEG[5,6]<-parameterEstimates(fit_sem_SA_entorhinal)[6,7]
data_subsample_GMlat_85_GGSEG[5,7]<-lavInspect(fit_sem_SA_entorhinal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[5,] <-parameterEstimates(fit_sem_SA_entorhinal,standardized = TRUE, rsquare=TRUE)[6,]

##fusiform
Model_SA_fusiform <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fusiform_sa'
fit_sem_SA_fusiform <- sem(Model_SA_fusiform, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_fusiform, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[6,5]<-BIC(fit_sem_SA_fusiform)
data_subsample_GMlat_85_GGSEG[6,6]<-parameterEstimates(fit_sem_SA_fusiform)[6,7]
data_subsample_GMlat_85_GGSEG[6,7]<-lavInspect(fit_sem_SA_fusiform,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[6,] <-parameterEstimates(fit_sem_SA_fusiform,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorparietal
Model_SA_inferiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorparietal_sa'
fit_sem_SA_inferiorparietal <- sem(Model_SA_inferiorparietal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_inferiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[7,5]<-BIC(fit_sem_SA_inferiorparietal)
data_subsample_GMlat_85_GGSEG[7,6]<-parameterEstimates(fit_sem_SA_inferiorparietal)[6,7]
data_subsample_GMlat_85_GGSEG[7,7]<-lavInspect(fit_sem_SA_inferiorparietal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[7,] <-parameterEstimates(fit_sem_SA_inferiorparietal,standardized = TRUE, rsquare=TRUE)[6,]

##inferiortemporal
Model_SA_inferiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiortemporal_sa'
fit_sem_SA_inferiortemporal <- sem(Model_SA_inferiortemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_inferiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[8,5]<-BIC(fit_sem_SA_inferiortemporal)
data_subsample_GMlat_85_GGSEG[8,6]<-parameterEstimates(fit_sem_SA_inferiortemporal)[6,7]
data_subsample_GMlat_85_GGSEG[8,7]<-lavInspect(fit_sem_SA_inferiortemporal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[8,] <-parameterEstimates(fit_sem_SA_inferiortemporal,standardized = TRUE, rsquare=TRUE)[6,]

##isthmuscingulate
Model_SA_isthmuscingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ isthmuscingulate_sa'
fit_sem_SA_isthmuscingulate <- sem(Model_SA_isthmuscingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_isthmuscingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[9,5]<-BIC(fit_sem_SA_isthmuscingulate)
data_subsample_GMlat_85_GGSEG[9,6]<-parameterEstimates(fit_sem_SA_isthmuscingulate)[6,7]
data_subsample_GMlat_85_GGSEG[9,7]<-lavInspect(fit_sem_SA_isthmuscingulate,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[9,] <-parameterEstimates(fit_sem_SA_isthmuscingulate,standardized = TRUE, rsquare=TRUE)[6,]

##lateraloccipital
Model_SA_lateraloccipital <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lateraloccipital_sa'
fit_sem_SA_lateraloccipital <- sem(Model_SA_lateraloccipital, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_lateraloccipital, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[10,5]<-BIC(fit_sem_SA_lateraloccipital)
data_subsample_GMlat_85_GGSEG[10,6]<-parameterEstimates(fit_sem_SA_lateraloccipital)[6,7]
data_subsample_GMlat_85_GGSEG[10,7]<-lavInspect(fit_sem_SA_lateraloccipital,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[10,] <-parameterEstimates(fit_sem_SA_lateraloccipital,standardized = TRUE, rsquare=TRUE)[6,]

##lateralorbitofrontal
Model_SA_lateralorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lateralorbitofrontal_sa'
fit_sem_SA_lateralorbitofrontal <- sem(Model_SA_lateralorbitofrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_lateralorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[11,5]<-BIC(fit_sem_SA_lateralorbitofrontal)
data_subsample_GMlat_85_GGSEG[11,6]<-parameterEstimates(fit_sem_SA_lateralorbitofrontal)[6,7]
data_subsample_GMlat_85_GGSEG[11,7]<-lavInspect(fit_sem_SA_lateralorbitofrontal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[11,] <-parameterEstimates(fit_sem_SA_lateralorbitofrontal,standardized = TRUE, rsquare=TRUE)[6,]

##lingual
Model_SA_lingual <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lingual_sa'
fit_sem_SA_lingual <- sem(Model_SA_lingual, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_lingual, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[12,5]<-BIC(fit_sem_SA_lingual)
data_subsample_GMlat_85_GGSEG[12,6]<-parameterEstimates(fit_sem_SA_lingual)[6,7]
data_subsample_GMlat_85_GGSEG[12,7]<-lavInspect(fit_sem_SA_lingual,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[12,] <-parameterEstimates(fit_sem_SA_lingual,standardized = TRUE, rsquare=TRUE)[6,]

##medialorbitofrontal
Model_SA_medialorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ medialorbitofrontal_sa'
fit_sem_SA_medialorbitofrontal <- sem(Model_SA_medialorbitofrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_medialorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[13,5]<-BIC(fit_sem_SA_medialorbitofrontal)
data_subsample_GMlat_85_GGSEG[13,6]<-parameterEstimates(fit_sem_SA_medialorbitofrontal)[6,7]
data_subsample_GMlat_85_GGSEG[13,7]<-lavInspect(fit_sem_SA_medialorbitofrontal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[13,] <-parameterEstimates(fit_sem_SA_medialorbitofrontal,standardized = TRUE, rsquare=TRUE)[6,]

##middletemporal
Model_SA_middletemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ middletemporal_sa'
fit_sem_SA_middletemporal <- sem(Model_SA_middletemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_middletemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[14,5]<-BIC(fit_sem_SA_middletemporal)
data_subsample_GMlat_85_GGSEG[14,6]<-parameterEstimates(fit_sem_SA_middletemporal)[6,7]
data_subsample_GMlat_85_GGSEG[14,7]<-lavInspect(fit_sem_SA_middletemporal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[14,] <-parameterEstimates(fit_sem_SA_middletemporal,standardized = TRUE, rsquare=TRUE)[6,]

##parahippocampal
Model_SA_parahippocampal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parahippocampal_sa'
fit_sem_SA_parahippocampal <- sem(Model_SA_parahippocampal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_parahippocampal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[15,5]<-BIC(fit_sem_SA_parahippocampal)
data_subsample_GMlat_85_GGSEG[15,6]<-parameterEstimates(fit_sem_SA_parahippocampal)[6,7]
data_subsample_GMlat_85_GGSEG[15,7]<-lavInspect(fit_sem_SA_parahippocampal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[15,] <-parameterEstimates(fit_sem_SA_parahippocampal,standardized = TRUE, rsquare=TRUE)[6,]

##paracentral
Model_SA_paracentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ paracentral_sa'
fit_sem_SA_paracentral <- sem(Model_SA_paracentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_paracentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[16,5]<-BIC(fit_sem_SA_paracentral)
data_subsample_GMlat_85_GGSEG[16,6]<-parameterEstimates(fit_sem_SA_paracentral)[6,7]
data_subsample_GMlat_85_GGSEG[16,7]<-lavInspect(fit_sem_SA_paracentral,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[16,] <-parameterEstimates(fit_sem_SA_paracentral,standardized = TRUE, rsquare=TRUE)[6,]

##parsopercularis
Model_SA_parsopercularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parsopercularis_sa'
fit_sem_SA_parsopercularis <- sem(Model_SA_parsopercularis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_parsopercularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[17,5]<-BIC(fit_sem_SA_parsopercularis)
data_subsample_GMlat_85_GGSEG[17,6]<-parameterEstimates(fit_sem_SA_parsopercularis)[6,7]
data_subsample_GMlat_85_GGSEG[17,7]<-lavInspect(fit_sem_SA_parsopercularis,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[17,] <-parameterEstimates(fit_sem_SA_parsopercularis,standardized = TRUE, rsquare=TRUE)[6,]

##parsorbitalis
Model_SA_parsorbitalis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parsorbitalis_sa'
fit_sem_SA_parsorbitalis <- sem(Model_SA_parsorbitalis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_parsorbitalis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[18,5]<-BIC(fit_sem_SA_parsorbitalis)
data_subsample_GMlat_85_GGSEG[18,6]<-parameterEstimates(fit_sem_SA_parsorbitalis)[6,7]
data_subsample_GMlat_85_GGSEG[18,7]<-lavInspect(fit_sem_SA_parsorbitalis,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[18,] <-parameterEstimates(fit_sem_SA_parsorbitalis,standardized = TRUE, rsquare=TRUE)[6,]

##parstriangularis
Model_SA_parstriangularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parstriangularis_sa'
fit_sem_SA_parstriangularis <- sem(Model_SA_parstriangularis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_parstriangularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[19,5]<-BIC(fit_sem_SA_parstriangularis)
data_subsample_GMlat_85_GGSEG[19,6]<-parameterEstimates(fit_sem_SA_parstriangularis)[6,7]
data_subsample_GMlat_85_GGSEG[19,7]<-lavInspect(fit_sem_SA_parstriangularis,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[19,] <-parameterEstimates(fit_sem_SA_parstriangularis,standardized = TRUE, rsquare=TRUE)[6,]

##pericalcarine
Model_SA_pericalcarine <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ pericalcarine_sa'
fit_sem_SA_pericalcarine <- sem(Model_SA_pericalcarine, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_pericalcarine, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[20,5]<-BIC(fit_sem_SA_pericalcarine)
data_subsample_GMlat_85_GGSEG[20,6]<-parameterEstimates(fit_sem_SA_pericalcarine)[6,7]
data_subsample_GMlat_85_GGSEG[20,7]<-lavInspect(fit_sem_SA_pericalcarine,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[20,] <-parameterEstimates(fit_sem_SA_pericalcarine,standardized = TRUE, rsquare=TRUE)[6,]

##postcentral
Model_SA_postcentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ postcentral_sa'
fit_sem_SA_postcentral <- sem(Model_SA_postcentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_postcentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[21,5]<-BIC(fit_sem_SA_postcentral)
data_subsample_GMlat_85_GGSEG[21,6]<-parameterEstimates(fit_sem_SA_postcentral)[6,7]
data_subsample_GMlat_85_GGSEG[21,7]<-lavInspect(fit_sem_SA_postcentral,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[21,] <-parameterEstimates(fit_sem_SA_postcentral,standardized = TRUE, rsquare=TRUE)[6,]

##posteriorcingulate
Model_SA_posteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ posteriorcingulate_sa'
fit_sem_SA_posteriorcingulate <- sem(Model_SA_posteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_posteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[22,5]<-BIC(fit_sem_SA_posteriorcingulate)
data_subsample_GMlat_85_GGSEG[22,6]<-parameterEstimates(fit_sem_SA_posteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[22,7]<-lavInspect(fit_sem_SA_posteriorcingulate,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[22,] <-parameterEstimates(fit_sem_SA_posteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##precentral
Model_SA_precentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ precentral_sa'
fit_sem_SA_precentral <- sem(Model_SA_precentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_precentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[23,5]<-BIC(fit_sem_SA_precentral)
data_subsample_GMlat_85_GGSEG[23,6]<-parameterEstimates(fit_sem_SA_precentral)[6,7]
data_subsample_GMlat_85_GGSEG[23,7]<-lavInspect(fit_sem_SA_precentral,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[23,] <-parameterEstimates(fit_sem_SA_precentral,standardized = TRUE, rsquare=TRUE)[6,]

##precuneus
Model_SA_precuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ precuneus_sa'
fit_sem_SA_precuneus <- sem(Model_SA_precuneus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_precuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[24,5]<-BIC(fit_sem_SA_precuneus)
data_subsample_GMlat_85_GGSEG[24,6]<-parameterEstimates(fit_sem_SA_precuneus)[6,7]
data_subsample_GMlat_85_GGSEG[24,7]<-lavInspect(fit_sem_SA_precuneus,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[24,] <-parameterEstimates(fit_sem_SA_precuneus,standardized = TRUE, rsquare=TRUE)[6,]

##rostralanteriorcingulate
Model_SA_rostralanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ rostralanteriorcingulate_sa'
fit_sem_SA_rostralanteriorcingulate <- sem(Model_SA_rostralanteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_rostralanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[25,5]<-BIC(fit_sem_SA_rostralanteriorcingulate)
data_subsample_GMlat_85_GGSEG[25,6]<-parameterEstimates(fit_sem_SA_rostralanteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[25,7]<-lavInspect(fit_sem_SA_rostralanteriorcingulate,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[25,] <-parameterEstimates(fit_sem_SA_rostralanteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##rostralmiddlefrontal
Model_SA_rostralmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ rostralmiddlefrontal_sa'
fit_sem_SA_rostralmiddlefrontal <- sem(Model_SA_rostralmiddlefrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_rostralmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[26,5]<-BIC(fit_sem_SA_rostralmiddlefrontal)
data_subsample_GMlat_85_GGSEG[26,6]<-parameterEstimates(fit_sem_SA_rostralmiddlefrontal)[6,7]
data_subsample_GMlat_85_GGSEG[26,7]<-lavInspect(fit_sem_SA_rostralmiddlefrontal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[26,] <-parameterEstimates(fit_sem_SA_rostralmiddlefrontal,standardized = TRUE, rsquare=TRUE)[6,]

##superiorfrontal
Model_SA_superiorfrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorfrontal_sa'
fit_sem_SA_superiorfrontal <- sem(Model_SA_superiorfrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_superiorfrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[27,5]<-BIC(fit_sem_SA_superiorfrontal)
data_subsample_GMlat_85_GGSEG[27,6]<-parameterEstimates(fit_sem_SA_superiorfrontal)[6,7]
data_subsample_GMlat_85_GGSEG[27,7]<-lavInspect(fit_sem_SA_superiorfrontal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[27,] <-parameterEstimates(fit_sem_SA_superiorfrontal,standardized = TRUE, rsquare=TRUE)[6,]

##superiorparietal
Model_SA_superiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorparietal_sa'
fit_sem_SA_superiorparietal <- sem(Model_SA_superiorparietal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_superiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[28,5]<-BIC(fit_sem_SA_superiorparietal)
data_subsample_GMlat_85_GGSEG[28,6]<-parameterEstimates(fit_sem_SA_superiorparietal)[6,7]
data_subsample_GMlat_85_GGSEG[28,7]<-lavInspect(fit_sem_SA_superiorparietal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[28,] <-parameterEstimates(fit_sem_SA_superiorparietal,standardized = TRUE, rsquare=TRUE)[6,]

##superiortemporal
Model_SA_superiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiortemporal_sa'
fit_sem_SA_superiortemporal <- sem(Model_SA_superiortemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_superiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[29,5]<-BIC(fit_sem_SA_superiortemporal)
data_subsample_GMlat_85_GGSEG[29,6]<-parameterEstimates(fit_sem_SA_superiortemporal)[6,7]
data_subsample_GMlat_85_GGSEG[29,7]<-lavInspect(fit_sem_SA_superiortemporal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[29,] <-parameterEstimates(fit_sem_SA_superiortemporal,standardized = TRUE, rsquare=TRUE)[6,]

##supramarginal
Model_SA_supramarginal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ supramarginal_sa'
fit_sem_SA_supramarginal <- sem(Model_SA_supramarginal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_supramarginal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[30,5]<-BIC(fit_sem_SA_supramarginal)
data_subsample_GMlat_85_GGSEG[30,6]<-parameterEstimates(fit_sem_SA_supramarginal)[6,7]
data_subsample_GMlat_85_GGSEG[30,7]<-lavInspect(fit_sem_SA_supramarginal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[30,] <-parameterEstimates(fit_sem_SA_supramarginal,standardized = TRUE, rsquare=TRUE)[6,]

##frontalpole
Model_SA_frontalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ frontalpole_sa'
fit_sem_SA_frontalpole <- sem(Model_SA_frontalpole, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_frontalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[31,5]<-BIC(fit_sem_SA_frontalpole)
data_subsample_GMlat_85_GGSEG[31,6]<-parameterEstimates(fit_sem_SA_frontalpole)[6,7]
data_subsample_GMlat_85_GGSEG[31,7]<-lavInspect(fit_sem_SA_frontalpole,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[31,] <-parameterEstimates(fit_sem_SA_frontalpole,standardized = TRUE, rsquare=TRUE)[6,]

##temporalpole
Model_SA_temporalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalpole_sa'
fit_sem_SA_temporalpole <- sem(Model_SA_temporalpole, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_temporalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[32,5]<-BIC(fit_sem_SA_temporalpole)
data_subsample_GMlat_85_GGSEG[32,6]<-parameterEstimates(fit_sem_SA_temporalpole)[6,7]
data_subsample_GMlat_85_GGSEG[32,7]<-lavInspect(fit_sem_SA_temporalpole,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[32,] <-parameterEstimates(fit_sem_SA_temporalpole,standardized = TRUE, rsquare=TRUE)[6,]

##transversetemporal
Model_SA_transversetemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ transversetemporal_sa'
fit_sem_SA_transversetemporal <- sem(Model_SA_transversetemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_transversetemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[33,5]<-BIC(fit_sem_SA_transversetemporal)
data_subsample_GMlat_85_GGSEG[33,6]<-parameterEstimates(fit_sem_SA_transversetemporal)[6,7]
data_subsample_GMlat_85_GGSEG[33,7]<-lavInspect(fit_sem_SA_transversetemporal,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[33,] <-parameterEstimates(fit_sem_SA_transversetemporal,standardized = TRUE, rsquare=TRUE)[6,]

##insula
Model_SA_insula <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ insula_sa'
fit_sem_SA_insula <- sem(Model_SA_insula, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SA_insula, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[34,5]<-BIC(fit_sem_SA_insula)
data_subsample_GMlat_85_GGSEG[34,6]<-parameterEstimates(fit_sem_SA_insula)[6,7]
data_subsample_GMlat_85_GGSEG[34,7]<-lavInspect(fit_sem_SA_insula,what = "std.all")$beta[1,2]
paramest_GMlat_SA_85[34,] <-parameterEstimates(fit_sem_SA_insula,standardized = TRUE, rsquare=TRUE)[6,]

write.csv(paramest_GMlat_SA_85, "paramest_GMlat_SA_85.csv", row.names=FALSE)

#Grey Matter Volume (34 regions)
##bankssts
Model_GMV_bankssts <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ bankssts_gmv'
fit_sem_GMV_bankssts <- sem(Model_GMV_bankssts, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_bankssts, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[1,8]<-BIC(fit_sem_GMV_bankssts)
data_subsample_GMlat_85_GGSEG[1,9]<-parameterEstimates(fit_sem_GMV_bankssts)[6,7]
data_subsample_GMlat_85_GGSEG[1,10]<-lavInspect(fit_sem_GMV_bankssts,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[1,] <-parameterEstimates(fit_sem_GMV_bankssts,standardized = TRUE, rsquare=TRUE)[6,]

##caudalanteriorcingulate
Model_GMV_caudalanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ caudalanteriorcingulate_gmv'
fit_sem_GMV_caudalanteriorcingulate <- sem(Model_GMV_caudalanteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_caudalanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[2,8]<-BIC(fit_sem_GMV_caudalanteriorcingulate)
data_subsample_GMlat_85_GGSEG[2,9]<-parameterEstimates(fit_sem_GMV_caudalanteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[2,10]<-lavInspect(fit_sem_GMV_caudalanteriorcingulate,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[2,] <-parameterEstimates(fit_sem_GMV_caudalanteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##caudalmiddlefrontal
Model_GMV_caudalmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ caudalmiddlefrontal_gmv'
fit_sem_GMV_caudalmiddlefrontal <- sem(Model_GMV_caudalmiddlefrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_caudalmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[3,8]<-BIC(fit_sem_GMV_caudalmiddlefrontal)
data_subsample_GMlat_85_GGSEG[3,9]<-parameterEstimates(fit_sem_GMV_caudalmiddlefrontal)[6,7]
data_subsample_GMlat_85_GGSEG[3,10]<-lavInspect(fit_sem_GMV_caudalmiddlefrontal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[3,] <-parameterEstimates(fit_sem_GMV_caudalmiddlefrontal,standardized = TRUE, rsquare=TRUE)[6,]

##cuneus
Model_GMV_cuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ cuneus_gmv'
fit_sem_GMV_cuneus <- sem(Model_GMV_cuneus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_cuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[4,8]<-BIC(fit_sem_GMV_cuneus)
data_subsample_GMlat_85_GGSEG[4,9]<-parameterEstimates(fit_sem_GMV_cuneus)[6,7]
data_subsample_GMlat_85_GGSEG[4,10]<-lavInspect(fit_sem_GMV_cuneus,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[4,] <-parameterEstimates(fit_sem_GMV_cuneus,standardized = TRUE, rsquare=TRUE)[6,]

##entorhinal
Model_GMV_entorhinal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ entorhinal_gmv'
fit_sem_GMV_entorhinal <- sem(Model_GMV_entorhinal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_entorhinal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[5,8]<-BIC(fit_sem_GMV_entorhinal)
data_subsample_GMlat_85_GGSEG[5,9]<-parameterEstimates(fit_sem_GMV_entorhinal)[6,7]
data_subsample_GMlat_85_GGSEG[5,10]<-lavInspect(fit_sem_GMV_entorhinal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[5,] <-parameterEstimates(fit_sem_GMV_entorhinal,standardized = TRUE, rsquare=TRUE)[6,]

##fusiform
Model_GMV_fusiform <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fusiform_gmv'
fit_sem_GMV_fusiform <- sem(Model_GMV_fusiform, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_fusiform, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[6,8]<-BIC(fit_sem_GMV_fusiform)
data_subsample_GMlat_85_GGSEG[6,9]<-parameterEstimates(fit_sem_GMV_fusiform)[6,7]
data_subsample_GMlat_85_GGSEG[6,10]<-lavInspect(fit_sem_GMV_fusiform,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[6,] <-parameterEstimates(fit_sem_GMV_fusiform,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorparietal
Model_GMV_inferiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorparietal_gmv'
fit_sem_GMV_inferiorparietal <- sem(Model_GMV_inferiorparietal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_inferiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[7,8]<-BIC(fit_sem_GMV_inferiorparietal)
data_subsample_GMlat_85_GGSEG[7,9]<-parameterEstimates(fit_sem_GMV_inferiorparietal)[6,7]
data_subsample_GMlat_85_GGSEG[7,10]<-lavInspect(fit_sem_GMV_inferiorparietal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[7,] <-parameterEstimates(fit_sem_GMV_inferiorparietal,standardized = TRUE, rsquare=TRUE)[6,]

##inferiortemporal
Model_GMV_inferiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiortemporal_gmv'
fit_sem_GMV_inferiortemporal <- sem(Model_GMV_inferiortemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_inferiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[8,8]<-BIC(fit_sem_GMV_inferiortemporal)
data_subsample_GMlat_85_GGSEG[8,9]<-parameterEstimates(fit_sem_GMV_inferiortemporal)[6,7]
data_subsample_GMlat_85_GGSEG[8,10]<-lavInspect(fit_sem_GMV_inferiortemporal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[8,] <-parameterEstimates(fit_sem_GMV_inferiortemporal,standardized = TRUE, rsquare=TRUE)[6,]

##isthmuscingulate
Model_GMV_isthmuscingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ isthmuscingulate_gmv'
fit_sem_GMV_isthmuscingulate <- sem(Model_GMV_isthmuscingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_isthmuscingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[9,8]<-BIC(fit_sem_GMV_isthmuscingulate)
data_subsample_GMlat_85_GGSEG[9,9]<-parameterEstimates(fit_sem_GMV_isthmuscingulate)[6,7]
data_subsample_GMlat_85_GGSEG[9,10]<-lavInspect(fit_sem_GMV_isthmuscingulate,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[9,] <-parameterEstimates(fit_sem_GMV_isthmuscingulate,standardized = TRUE, rsquare=TRUE)[6,]

##lateraloccipital
Model_GMV_lateraloccipital <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lateraloccipital_gmv'
fit_sem_GMV_lateraloccipital <- sem(Model_GMV_lateraloccipital, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_lateraloccipital, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[10,8]<-BIC(fit_sem_GMV_lateraloccipital)
data_subsample_GMlat_85_GGSEG[10,9]<-parameterEstimates(fit_sem_GMV_lateraloccipital)[6,7]
data_subsample_GMlat_85_GGSEG[10,10]<-lavInspect(fit_sem_GMV_lateraloccipital,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[10,] <-parameterEstimates(fit_sem_GMV_lateraloccipital,standardized = TRUE, rsquare=TRUE)[6,]

##lateralorbitofrontal
Model_GMV_lateralorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lateralorbitofrontal_gmv'
fit_sem_GMV_lateralorbitofrontal <- sem(Model_GMV_lateralorbitofrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_lateralorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[11,8]<-BIC(fit_sem_GMV_lateralorbitofrontal)
data_subsample_GMlat_85_GGSEG[11,9]<-parameterEstimates(fit_sem_GMV_lateralorbitofrontal)[6,7]
data_subsample_GMlat_85_GGSEG[11,10]<-lavInspect(fit_sem_GMV_lateralorbitofrontal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[11,] <-parameterEstimates(fit_sem_GMV_lateralorbitofrontal,standardized = TRUE, rsquare=TRUE)[6,]

##lingual
Model_GMV_lingual <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ lingual_gmv'
fit_sem_GMV_lingual <- sem(Model_GMV_lingual, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_lingual, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[12,8]<-BIC(fit_sem_GMV_lingual)
data_subsample_GMlat_85_GGSEG[12,9]<-parameterEstimates(fit_sem_GMV_lingual)[6,7]
data_subsample_GMlat_85_GGSEG[12,10]<-lavInspect(fit_sem_GMV_lingual,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[12,] <-parameterEstimates(fit_sem_GMV_lingual,standardized = TRUE, rsquare=TRUE)[6,]

##medialorbitofrontal
Model_GMV_medialorbitofrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ medialorbitofrontal_gmv'
fit_sem_GMV_medialorbitofrontal <- sem(Model_GMV_medialorbitofrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_medialorbitofrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[13,8]<-BIC(fit_sem_GMV_medialorbitofrontal)
data_subsample_GMlat_85_GGSEG[13,9]<-parameterEstimates(fit_sem_GMV_medialorbitofrontal)[6,7]
data_subsample_GMlat_85_GGSEG[13,10]<-lavInspect(fit_sem_GMV_medialorbitofrontal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[13,] <-parameterEstimates(fit_sem_GMV_medialorbitofrontal,standardized = TRUE, rsquare=TRUE)[6,]

##middletemporal
Model_GMV_middletemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ middletemporal_gmv'
fit_sem_GMV_middletemporal <- sem(Model_GMV_middletemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_middletemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[14,8]<-BIC(fit_sem_GMV_middletemporal)
data_subsample_GMlat_85_GGSEG[14,9]<-parameterEstimates(fit_sem_GMV_middletemporal)[6,7]
data_subsample_GMlat_85_GGSEG[14,10]<-lavInspect(fit_sem_GMV_middletemporal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[14,] <-parameterEstimates(fit_sem_GMV_middletemporal,standardized = TRUE, rsquare=TRUE)[6,]

##parahippocampal
Model_GMV_parahippocampal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parahippocampal_gmv'
fit_sem_GMV_parahippocampal <- sem(Model_GMV_parahippocampal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_parahippocampal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[15,8]<-BIC(fit_sem_GMV_parahippocampal)
data_subsample_GMlat_85_GGSEG[15,9]<-parameterEstimates(fit_sem_GMV_parahippocampal)[6,7]
data_subsample_GMlat_85_GGSEG[15,10]<-lavInspect(fit_sem_GMV_parahippocampal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[15,] <-parameterEstimates(fit_sem_GMV_parahippocampal,standardized = TRUE, rsquare=TRUE)[6,]

##paracentral
Model_GMV_paracentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ paracentral_gmv'
fit_sem_GMV_paracentral <- sem(Model_GMV_paracentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_paracentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[16,8]<-BIC(fit_sem_GMV_paracentral)
data_subsample_GMlat_85_GGSEG[16,9]<-parameterEstimates(fit_sem_GMV_paracentral)[6,7]
data_subsample_GMlat_85_GGSEG[16,10]<-lavInspect(fit_sem_GMV_paracentral,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[16,] <-parameterEstimates(fit_sem_GMV_paracentral,standardized = TRUE, rsquare=TRUE)[6,]

##parsopercularis
Model_GMV_parsopercularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parsopercularis_gmv'
fit_sem_GMV_parsopercularis <- sem(Model_GMV_parsopercularis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_parsopercularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[17,8]<-BIC(fit_sem_GMV_parsopercularis)
data_subsample_GMlat_85_GGSEG[17,9]<-parameterEstimates(fit_sem_GMV_parsopercularis)[6,7]
data_subsample_GMlat_85_GGSEG[17,10]<-lavInspect(fit_sem_GMV_parsopercularis,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[17,] <-parameterEstimates(fit_sem_GMV_parsopercularis,standardized = TRUE, rsquare=TRUE)[6,]

##parsorbitalis
Model_GMV_parsorbitalis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parsorbitalis_gmv'
fit_sem_GMV_parsorbitalis <- sem(Model_GMV_parsorbitalis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_parsorbitalis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[18,8]<-BIC(fit_sem_GMV_parsorbitalis)
data_subsample_GMlat_85_GGSEG[18,9]<-parameterEstimates(fit_sem_GMV_parsorbitalis)[6,7]
data_subsample_GMlat_85_GGSEG[18,10]<-lavInspect(fit_sem_GMV_parsorbitalis,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[18,] <-parameterEstimates(fit_sem_GMV_parsorbitalis,standardized = TRUE, rsquare=TRUE)[6,]

##parstriangularis
Model_GMV_parstriangularis <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parstriangularis_gmv'
fit_sem_GMV_parstriangularis <- sem(Model_GMV_parstriangularis, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_parstriangularis, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[19,8]<-BIC(fit_sem_GMV_parstriangularis)
data_subsample_GMlat_85_GGSEG[19,9]<-parameterEstimates(fit_sem_GMV_parstriangularis)[6,7]
data_subsample_GMlat_85_GGSEG[19,10]<-lavInspect(fit_sem_GMV_parstriangularis,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[19,] <-parameterEstimates(fit_sem_GMV_parstriangularis,standardized = TRUE, rsquare=TRUE)[6,]

##pericalcarine
Model_GMV_pericalcarine <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ pericalcarine_gmv'
fit_sem_GMV_pericalcarine <- sem(Model_GMV_pericalcarine, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_pericalcarine, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[20,8]<-BIC(fit_sem_GMV_pericalcarine)
data_subsample_GMlat_85_GGSEG[20,9]<-parameterEstimates(fit_sem_GMV_pericalcarine)[6,7]
data_subsample_GMlat_85_GGSEG[20,10]<-lavInspect(fit_sem_GMV_pericalcarine,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[20,] <-parameterEstimates(fit_sem_GMV_pericalcarine,standardized = TRUE, rsquare=TRUE)[6,]

##postcentral
Model_GMV_postcentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ postcentral_gmv'
fit_sem_GMV_postcentral <- sem(Model_GMV_postcentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_postcentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[21,8]<-BIC(fit_sem_GMV_postcentral)
data_subsample_GMlat_85_GGSEG[21,9]<-parameterEstimates(fit_sem_GMV_postcentral)[6,7]
data_subsample_GMlat_85_GGSEG[21,10]<-lavInspect(fit_sem_GMV_postcentral,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[21,] <-parameterEstimates(fit_sem_GMV_postcentral,standardized = TRUE, rsquare=TRUE)[6,]

##posteriorcingulate
Model_GMV_posteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ posteriorcingulate_gmv'
fit_sem_GMV_posteriorcingulate <- sem(Model_GMV_posteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_posteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[22,8]<-BIC(fit_sem_GMV_posteriorcingulate)
data_subsample_GMlat_85_GGSEG[22,9]<-parameterEstimates(fit_sem_GMV_posteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[22,10]<-lavInspect(fit_sem_GMV_posteriorcingulate,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[22,] <-parameterEstimates(fit_sem_GMV_posteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##precentral
Model_GMV_precentral <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ precentral_gmv'
fit_sem_GMV_precentral <- sem(Model_GMV_precentral, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_precentral, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[23,8]<-BIC(fit_sem_GMV_precentral)
data_subsample_GMlat_85_GGSEG[23,9]<-parameterEstimates(fit_sem_GMV_precentral)[6,7]
data_subsample_GMlat_85_GGSEG[23,10]<-lavInspect(fit_sem_GMV_precentral,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[23,] <-parameterEstimates(fit_sem_GMV_precentral,standardized = TRUE, rsquare=TRUE)[6,]

##precuneus
Model_GMV_precuneus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ precuneus_gmv'
fit_sem_GMV_precuneus <- sem(Model_GMV_precuneus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_precuneus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[24,8]<-BIC(fit_sem_GMV_precuneus)
data_subsample_GMlat_85_GGSEG[24,9]<-parameterEstimates(fit_sem_GMV_precuneus)[6,7]
data_subsample_GMlat_85_GGSEG[24,10]<-lavInspect(fit_sem_GMV_precuneus,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[24,] <-parameterEstimates(fit_sem_GMV_precuneus,standardized = TRUE, rsquare=TRUE)[6,]

##rostralanteriorcingulate
Model_GMV_rostralanteriorcingulate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ rostralanteriorcingulate_gmv'
fit_sem_GMV_rostralanteriorcingulate <- sem(Model_GMV_rostralanteriorcingulate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_rostralanteriorcingulate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[25,8]<-BIC(fit_sem_GMV_rostralanteriorcingulate)
data_subsample_GMlat_85_GGSEG[25,9]<-parameterEstimates(fit_sem_GMV_rostralanteriorcingulate)[6,7]
data_subsample_GMlat_85_GGSEG[25,10]<-lavInspect(fit_sem_GMV_rostralanteriorcingulate,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[25,] <-parameterEstimates(fit_sem_GMV_rostralanteriorcingulate,standardized = TRUE, rsquare=TRUE)[6,]

##rostralmiddlefrontal
Model_GMV_rostralmiddlefrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ rostralmiddlefrontal_gmv'
fit_sem_GMV_rostralmiddlefrontal <- sem(Model_GMV_rostralmiddlefrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_rostralmiddlefrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[26,8]<-BIC(fit_sem_GMV_rostralmiddlefrontal)
data_subsample_GMlat_85_GGSEG[26,9]<-parameterEstimates(fit_sem_GMV_rostralmiddlefrontal)[6,7]
data_subsample_GMlat_85_GGSEG[26,10]<-lavInspect(fit_sem_GMV_rostralmiddlefrontal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[26,] <-parameterEstimates(fit_sem_GMV_rostralmiddlefrontal,standardized = TRUE, rsquare=TRUE)[6,]

##superiorfrontal
Model_GMV_superiorfrontal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorfrontal_gmv'
fit_sem_GMV_superiorfrontal <- sem(Model_GMV_superiorfrontal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_superiorfrontal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[27,8]<-BIC(fit_sem_GMV_superiorfrontal)
data_subsample_GMlat_85_GGSEG[27,9]<-parameterEstimates(fit_sem_GMV_superiorfrontal)[6,7]
data_subsample_GMlat_85_GGSEG[27,10]<-lavInspect(fit_sem_GMV_superiorfrontal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[27,] <-parameterEstimates(fit_sem_GMV_superiorfrontal,standardized = TRUE, rsquare=TRUE)[6,]

##superiorparietal
Model_GMV_superiorparietal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorparietal_gmv'
fit_sem_GMV_superiorparietal <- sem(Model_GMV_superiorparietal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_superiorparietal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[28,8]<-BIC(fit_sem_GMV_superiorparietal)
data_subsample_GMlat_85_GGSEG[28,9]<-parameterEstimates(fit_sem_GMV_superiorparietal)[6,7]
data_subsample_GMlat_85_GGSEG[28,10]<-lavInspect(fit_sem_GMV_superiorparietal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[28,] <-parameterEstimates(fit_sem_GMV_superiorparietal,standardized = TRUE, rsquare=TRUE)[6,]

##superiortemporal
Model_GMV_superiortemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiortemporal_gmv'
fit_sem_GMV_superiortemporal <- sem(Model_GMV_superiortemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_superiortemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[29,8]<-BIC(fit_sem_GMV_superiortemporal)
data_subsample_GMlat_85_GGSEG[29,9]<-parameterEstimates(fit_sem_GMV_superiortemporal)[6,7]
data_subsample_GMlat_85_GGSEG[29,10]<-lavInspect(fit_sem_GMV_superiortemporal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[29,] <-parameterEstimates(fit_sem_GMV_superiortemporal,standardized = TRUE, rsquare=TRUE)[6,]

##supramarginal
Model_GMV_supramarginal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ supramarginal_gmv'
fit_sem_GMV_supramarginal <- sem(Model_GMV_supramarginal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_supramarginal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[30,8]<-BIC(fit_sem_GMV_supramarginal)
data_subsample_GMlat_85_GGSEG[30,9]<-parameterEstimates(fit_sem_GMV_supramarginal)[6,7]
data_subsample_GMlat_85_GGSEG[30,10]<-lavInspect(fit_sem_GMV_supramarginal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[30,] <-parameterEstimates(fit_sem_GMV_supramarginal,standardized = TRUE, rsquare=TRUE)[6,]

##frontalpole
Model_GMV_frontalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ frontalpole_gmv'
fit_sem_GMV_frontalpole <- sem(Model_GMV_frontalpole, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_frontalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[31,8]<-BIC(fit_sem_GMV_frontalpole)
data_subsample_GMlat_85_GGSEG[31,9]<-parameterEstimates(fit_sem_GMV_frontalpole)[6,7]
data_subsample_GMlat_85_GGSEG[31,10]<-lavInspect(fit_sem_GMV_frontalpole,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[31,] <-parameterEstimates(fit_sem_GMV_frontalpole,standardized = TRUE, rsquare=TRUE)[6,]

##temporalpole
Model_GMV_temporalpole <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalpole_gmv'
fit_sem_GMV_temporalpole <- sem(Model_GMV_temporalpole, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_temporalpole, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[32,8]<-BIC(fit_sem_GMV_temporalpole)
data_subsample_GMlat_85_GGSEG[32,9]<-parameterEstimates(fit_sem_GMV_temporalpole)[6,7]
data_subsample_GMlat_85_GGSEG[32,10]<-lavInspect(fit_sem_GMV_temporalpole,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[32,] <-parameterEstimates(fit_sem_GMV_temporalpole,standardized = TRUE, rsquare=TRUE)[6,]

##transversetemporal
Model_GMV_transversetemporal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ transversetemporal_gmv'
fit_sem_GMV_transversetemporal <- sem(Model_GMV_transversetemporal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_transversetemporal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[33,8]<-BIC(fit_sem_GMV_transversetemporal)
data_subsample_GMlat_85_GGSEG[33,9]<-parameterEstimates(fit_sem_GMV_transversetemporal)[6,7]
data_subsample_GMlat_85_GGSEG[33,10]<-lavInspect(fit_sem_GMV_transversetemporal,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[33,] <-parameterEstimates(fit_sem_GMV_transversetemporal,standardized = TRUE, rsquare=TRUE)[6,]

##insula
Model_GMV_insula <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ insula_gmv'
fit_sem_GMV_insula <- sem(Model_GMV_insula, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMV_insula, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_GMlat_85_GGSEG[34,8]<-BIC(fit_sem_GMV_insula)
data_subsample_GMlat_85_GGSEG[34,9]<-parameterEstimates(fit_sem_GMV_insula)[6,7]
data_subsample_GMlat_85_GGSEG[34,10]<-lavInspect(fit_sem_GMV_insula,what = "std.all")$beta[1,2]
paramest_GMlat_GMV_85[34,] <-parameterEstimates(fit_sem_GMV_insula,standardized = TRUE, rsquare=TRUE)[6,]

write.csv(paramest_GMlat_GMV_85, "paramest_GMlat_GMV_85.csv", row.names=FALSE)

#Data Visualization for the individual models in the three grey matter metrics
data_subsample_GMlat_85_rl_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                 CT_fit=rep(NA, times=68),
                                                 CT_p=rep(NA, times=68),
                                                 CT_std=rep(NA, times=68),
                                                 SA_fit=rep(NA, times=68),
                                                 SA_p=rep(NA, times=68),
                                                 SA_std=rep(NA, times=68),
                                                 GMV_fit=rep(NA, times=68),
                                                 GMV_p=rep(NA, times=68),
                                                 GMV_std=rep(NA, times=68))

for (j in 2:10) {
  a=1
  for (i in seq(from = 1, to = 68, by = 2)) {
    data_subsample_GMlat_85_rl_GGSEG[i,j]<-data_subsample_GMlat_85_GGSEG[a,j]
    data_subsample_GMlat_85_rl_GGSEG[i+1,j]<-data_subsample_GMlat_85_GGSEG[a,j]
    a=a+1
  }
}

data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_fit)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="Fit (BIC)")
data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_fit)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="Fit (BIC)")
data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_fit)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="Fit (BIC)")
data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_rl_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="Standardized\nparameter\nestimates")

#Select all the std for the grey matter metrics
data_subsample_GMlat_85_GGSEG_CT <- data_subsample_GMlat_85_GGSEG %>%
  dplyr::select(-c("CT_fit","CT_p","SA_fit","SA_p","SA_std","GMV_fit","GMV_p","GMV_std")) %>%
  rename(c("CT_std" = "Std")) %>%
  add_column(Metric = "CT")

data_subsample_GMlat_85_GGSEG_SA <- data_subsample_GMlat_85_GGSEG %>%
  dplyr::select(-c("CT_fit","CT_p","CT_std","SA_fit","SA_p","GMV_fit","GMV_p","GMV_std")) %>%
  rename(c("SA_std" = "Std")) %>%
  add_column(Metric = "SA")

data_subsample_GMlat_85_GGSEG_GMV <- data_subsample_GMlat_85_GGSEG %>%
  dplyr::select(-c("CT_fit","CT_p","CT_std","SA_fit","SA_p","SA_std","GMV_fit","GMV_p")) %>%
  rename(c("GMV_std" = "Std")) %>%
  add_column(Metric = "GMV")

#Compute the maximum and the minimum
min(data_subsample_GMlat_85_GGSEG_CT[,2])
max(data_subsample_GMlat_85_GGSEG_CT[,2])
min(data_subsample_GMlat_85_GGSEG_SA[,2])
max(data_subsample_GMlat_85_GGSEG_SA[,2])
min(data_subsample_GMlat_85_GGSEG_GMV[,2])
max(data_subsample_GMlat_85_GGSEG_GMV[,2])

data_subsample_STD_GMlat_85 <- rbind(data_subsample_GMlat_85_GGSEG_CT,data_subsample_GMlat_85_GGSEG_SA,by=c("label"))
data_subsample_STD_GMlat_85 <- rbind(data_subsample_STD_GMlat_85,data_subsample_GMlat_85_GGSEG_GMV,by=c("label"))
data_subsample_STD_GMlat_85 <- data_subsample_STD_GMlat_85[-c(69,104),] 

data_subsample_STD_GMlat_85$mean <- ave(as.numeric(data_subsample_STD_GMlat_85$Std), data_subsample_STD_GMlat_85$Metric)

ggplot(data_subsample_STD_GMlat_85,aes(as.numeric(Std),fill=Metric)) + 
  geom_histogram(binwidth=0.005)+
  geom_density(adjust=2,alpha=0.2) +
  facet_grid(Metric~.) +
  geom_vline(aes(xintercept = mean, group = Metric), linetype="dashed", size=1,alpha=0.4) +
  theme_classic(base_size = 30)+
  xlab('Standardized estimated model parameters') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_fill_viridis(discrete=TRUE, option="rocket",begin=0.3,end=0.9,name="Metrics",labels=c("Cortical Thickness","Surface Area","Grey Matter Volume"))
ggsave("plot_GMlat_85_std.png")

#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Model per metric - cortical thickness
#-------------------------------------------------------------------------------------------------------------------------------
  
#Regularized cortical thickness model
Model_CTlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + fusiform_ct + isthmuscingulate_ct + lateraloccipital_ct + lateralorbitofrontal_ct + medialorbitofrontal_ct + parahippocampal_ct + paracentral_ct + parsopercularis_ct + parsorbitalis_ct + pericalcarine_ct + postcentral_ct + posteriorcingulate_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiorparietal_ct + superiortemporal_ct + transversetemporal_ct + insula_ct'

fit_sem_CTlat_85_reg_free <- sem(Model_CTlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CTlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized cortical thickness model
Model_CTlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*fusiform_ct + a*isthmuscingulate_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*medialorbitofrontal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsopercularis_ct + a*parsorbitalis_ct + a*pericalcarine_ct + a*postcentral_ct + a*posteriorcingulate_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*superiortemporal_ct + a*transversetemporal_ct + a*insula_ct'

fit_sem_CTlat_85_reg_constrained <- sem(Model_CTlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CTlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
kable(anova(fit_sem_CTlat_85_reg_free,fit_sem_CTlat_85_reg_constrained))

#Data Visualization for the regularized cortical thickness model
BIC(fit_sem_CTlat_85_reg_free)
data_subsample_CTlat_85_rl_reg_GGSEG <- data.frame(label=c("lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_fusiform","rh_fusiform","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                   CT_p=rep(NA, times=44),
                                                   CT_std=rep(NA, times=44))
a=1
for (i in seq(from = 1, to = 44, by = 2)) {
  data_subsample_CTlat_85_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_CTlat_85_reg_free)[a+5,7]
  data_subsample_CTlat_85_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_CTlat_85_reg_free)[a+5,7]
  data_subsample_CTlat_85_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_CTlat_85_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_CTlat_85_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_CTlat_85_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_CTlat_85_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_CTlat_85_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#Model for cortical thickness with all the regions of interest (34 ROIs)
Model_CTlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ bankssts_ct + caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + cuneus_ct + entorhinal_ct + fusiform_ct + inferiorparietal_ct + inferiortemporal_ct + isthmuscingulate_ct + lateraloccipital_ct + lateralorbitofrontal_ct + lingual_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsopercularis_ct + parsorbitalis_ct + parstriangularis_ct + pericalcarine_ct + postcentral_ct + posteriorcingulate_ct + precentral_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiorparietal_ct + superiortemporal_ct + supramarginal_ct + frontalpole_ct + temporalpole_ct + transversetemporal_ct + insula_ct'

fit_sem_CTlat_85_all_free <- sem(Model_CTlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CTlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the cortical thickness model with all the regions
Model_CTlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*bankssts_ct + a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*cuneus_ct + a*entorhinal_ct + a*fusiform_ct + a*inferiorparietal_ct + a*inferiortemporal_ct + a*isthmuscingulate_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*lingual_ct + a*medialorbitofrontal_ct + a*middletemporal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsopercularis_ct + a*parsorbitalis_ct + a*parstriangularis_ct + a*pericalcarine_ct + a*postcentral_ct + a*posteriorcingulate_ct + a*precentral_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*superiortemporal_ct + a*supramarginal_ct + a*frontalpole_ct + a*temporalpole_ct + a*transversetemporal_ct + a*insula_ct'

fit_sem_CTlat_85_all_constrained <- sem(Model_CTlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_CTlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)
kable(anova(fit_sem_CTlat_85_all_free,fit_sem_CTlat_85_all_constrained))

#Data Visualization for the cortical thickness model with all the regions
BIC(fit_sem_CTlat_85_all_free)
data_subsample_CTlat_85_rl_all_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                   CT_p=rep(NA, times=68),
                                                   CT_std=rep(NA, times=68))
a=1
for (i in seq(from = 1, to = 68, by = 2)) {
  data_subsample_CTlat_85_rl_all_GGSEG[i,2]<-parameterEstimates(fit_sem_CTlat_85_all_free)[a+5,7]
  data_subsample_CTlat_85_rl_all_GGSEG[i+1,2]<-parameterEstimates(fit_sem_CTlat_85_all_free)[a+5,7]
  data_subsample_CTlat_85_rl_all_GGSEG[i,3]<-lavInspect(fit_sem_CTlat_85_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_CTlat_85_rl_all_GGSEG[i+1,3]<-lavInspect(fit_sem_CTlat_85_all_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_CTlat_85_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_CTlat_85_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

#Create tables to save the parameters of the regularized model and the model with all the ROIs for cortical thickness
paramest_GMlat_CTreg_85 <- data.frame(lhs=rep(NA, times=22),
                                        op=rep(NA, times=22),
                                        rhs=rep(NA, times=22),
                                        est=rep(NA, times=22),
                                        se=rep(NA, times=22),
                                        z=rep(NA, times=22),
                                        pvalue=rep(NA, times=22),
                                        ci.lower=rep(NA, times=22),
                                        ci.upper=rep(NA, times=22),
                                        std.lv=rep(NA, times=22),
                                        std.all=rep(NA, times=22),
                                        std.nox=rep(NA, times=22))

paramest_GMlat_CTreg_85[,] <-parameterEstimates(fit_sem_CTlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:27,]
write.csv(paramest_GMlat_CTreg_85, "paramest_GMlat_CTreg_85.csv", row.names=FALSE)

paramest_GMlat_CTall_85 <- data.frame(lhs=rep(NA, times=34),
                                        op=rep(NA, times=34),
                                        rhs=rep(NA, times=34),
                                        est=rep(NA, times=34),
                                        se=rep(NA, times=34),
                                        z=rep(NA, times=34),
                                        pvalue=rep(NA, times=34),
                                        ci.lower=rep(NA, times=34),
                                        ci.upper=rep(NA, times=34),
                                        std.lv=rep(NA, times=34),
                                        std.all=rep(NA, times=34),
                                        std.nox=rep(NA, times=34))

paramest_GMlat_CTall_85[,] <-parameterEstimates(fit_sem_CTlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:39,]
write.csv(paramest_GMlat_CTall_85, "paramest_GMlat_CTall_85.csv", row.names=FALSE)


#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Model per metric - surface area
#-------------------------------------------------------------------------------------------------------------------------------
  
#Regularized surface area model
Model_SAlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ caudalanteriorcingulate_sa + cuneus_sa + entorhinal_sa + fusiform_sa + inferiorparietal_sa + isthmuscingulate_sa + lateraloccipital_sa + lingual_sa + parahippocampal_sa + paracentral_sa + parsopercularis_sa + parsorbitalis_sa + pericalcarine_sa + postcentral_sa + precentral_sa + precuneus_sa + rostralmiddlefrontal_sa + transversetemporal_sa'

fit_sem_SAlat_85_reg_free <- sem(Model_SAlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SAlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized surface area model
Model_SAlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*caudalanteriorcingulate_sa + a*cuneus_sa + a*entorhinal_sa + a*fusiform_sa + a*inferiorparietal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*lingual_sa + a*parahippocampal_sa + a*paracentral_sa + a*parsopercularis_sa + a*parsorbitalis_sa + a*pericalcarine_sa + a*postcentral_sa + a*precentral_sa + a*precuneus_sa + a*rostralmiddlefrontal_sa + a*transversetemporal_sa'

fit_sem_SAlat_85_reg_constrained <- sem(Model_SAlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SAlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_SAlat_85_reg_free,fit_sem_SAlat_85_reg_constrained))

#Data Visualization for the regularized surface area model
BIC(fit_sem_SAlat_85_reg_free) 
data_subsample_SAlat_85_rl_reg_GGSEG <- data.frame(label=c("lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lingual","rh_lingual","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_transversetemporal","rh_transversetemporal"),
                                                   SA_p=rep(NA, times=36),
                                                   SA_std=rep(NA, times=36))
a=1
for (i in seq(from = 1, to = 36, by = 2)) {
  data_subsample_SAlat_85_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_SAlat_85_reg_free)[a+5,7]
  data_subsample_SAlat_85_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_SAlat_85_reg_free)[a+5,7]
  data_subsample_SAlat_85_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_SAlat_85_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_SAlat_85_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_SAlat_85_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_SAlat_85_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_SAlat_85_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#Model for surface area with all the regions of interest (34 ROIs)
Model_SAlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ bankssts_sa + caudalanteriorcingulate_sa + caudalmiddlefrontal_sa + cuneus_sa + entorhinal_sa + fusiform_sa + inferiorparietal_sa + inferiortemporal_sa + isthmuscingulate_sa + lateraloccipital_sa + lateralorbitofrontal_sa + lingual_sa + medialorbitofrontal_sa + middletemporal_sa + parahippocampal_sa + paracentral_sa + parsopercularis_sa + parsorbitalis_sa + parstriangularis_sa + pericalcarine_sa + postcentral_sa + posteriorcingulate_sa + precentral_sa + precuneus_sa + rostralanteriorcingulate_sa + rostralmiddlefrontal_sa + superiorfrontal_sa + superiorparietal_sa + superiortemporal_sa + supramarginal_sa + frontalpole_sa + temporalpole_sa + transversetemporal_sa + insula_sa'

fit_sem_SAlat_85_all_free <- sem(Model_SAlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SAlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the surface area model with all the regions
Model_SAlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*bankssts_sa + a*caudalanteriorcingulate_sa + a*caudalmiddlefrontal_sa + a*cuneus_sa + a*entorhinal_sa + a*fusiform_sa + a*inferiorparietal_sa + a*inferiortemporal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*lateralorbitofrontal_sa + a*lingual_sa + a*medialorbitofrontal_sa + a*middletemporal_sa + a*parahippocampal_sa + a*paracentral_sa + a*parsopercularis_sa + a*parsorbitalis_sa + a*parstriangularis_sa + a*pericalcarine_sa + a*postcentral_sa + a*posteriorcingulate_sa + a*precentral_sa + a*precuneus_sa + a*rostralanteriorcingulate_sa + a*rostralmiddlefrontal_sa + a*superiorfrontal_sa + a*superiorparietal_sa + a*superiortemporal_sa + a*supramarginal_sa + a*frontalpole_sa + a*temporalpole_sa + a*transversetemporal_sa + a*insula_sa'

fit_sem_SAlat_85_all_constrained <- sem(Model_SAlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_SAlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_SAlat_85_all_free,fit_sem_SAlat_85_all_constrained))

#Data Visualization for the surface area model with all the regions
BIC(fit_sem_SAlat_85_all_free)
data_subsample_SAlat_85_rl_all_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                   SA_p=rep(NA, times=68),
                                                   SA_std=rep(NA, times=68))
a=1
for (i in seq(from = 1, to = 68, by = 2)) {
  data_subsample_SAlat_85_rl_all_GGSEG[i,2]<-parameterEstimates(fit_sem_SAlat_85_all_free)[a+5,7]
  data_subsample_SAlat_85_rl_all_GGSEG[i+1,2]<-parameterEstimates(fit_sem_SAlat_85_all_free)[a+5,7]
  data_subsample_SAlat_85_rl_all_GGSEG[i,3]<-lavInspect(fit_sem_SAlat_85_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_SAlat_85_rl_all_GGSEG[i+1,3]<-lavInspect(fit_sem_SAlat_85_all_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_SAlat_85_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_SAlat_85_rl_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#Create tables to save the parameters of the regularized model and the model with all the ROIs for surface area
paramest_GMlat_SAall_85 <- data.frame(lhs=rep(NA, times=34),
                                      op=rep(NA, times=34),
                                      rhs=rep(NA, times=34),
                                      est=rep(NA, times=34),
                                      se=rep(NA, times=34),
                                      z=rep(NA, times=34),
                                      pvalue=rep(NA, times=34),
                                      ci.lower=rep(NA, times=34),
                                      ci.upper=rep(NA, times=34),
                                      std.lv=rep(NA, times=34),
                                      std.all=rep(NA, times=34),
                                      std.nox=rep(NA, times=34))

paramest_GMlat_SAall_85[,] <-parameterEstimates(fit_sem_SAlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:39,]
write.csv(paramest_GMlat_SAall_85, "paramest_GMlat_SAall_85.csv", row.names=FALSE)

paramest_GMlat_SAreg_85 <- data.frame(lhs=rep(NA, times=18),
                                      op=rep(NA, times=18),
                                      rhs=rep(NA, times=18),
                                      est=rep(NA, times=18),
                                      se=rep(NA, times=18),
                                      z=rep(NA, times=18),
                                      pvalue=rep(NA, times=18),
                                      ci.lower=rep(NA, times=18),
                                      ci.upper=rep(NA, times=18),
                                      std.lv=rep(NA, times=18),
                                      std.all=rep(NA, times=18),
                                      std.nox=rep(NA, times=18))

paramest_GMlat_SAreg_85[,] <-parameterEstimates(fit_sem_SAlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:23,]
write.csv(paramest_GMlat_SAreg_85, "paramest_GMlat_SAreg_85.csv", row.names=FALSE)


#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Model per metric - grey matter volume
#-------------------------------------------------------------------------------------------------------------------------------
  
#Regularized grey matter volume model
Model_GMVlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ entorhinal_gmv + fusiform_gmv + inferiorparietal_gmv + lateraloccipital_gmv + parahippocampal_gmv + paracentral_gmv + pericalcarine_gmv + precentral_gmv + precuneus_gmv + supramarginal_gmv + insula_gmv'

fit_sem_GMVlat_85_reg_free <- sem(Model_GMVlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMVlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized grey matter volume model
Model_GMVlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                      Cognitive_factor ~ a*entorhinal_gmv + a*fusiform_gmv + a*inferiorparietal_gmv + a*lateraloccipital_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*pericalcarine_gmv + a*precentral_gmv + a*precuneus_gmv + a*supramarginal_gmv + a*insula_gmv'

fit_sem_GMVlat_85_reg_constrained <- sem(Model_GMVlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMVlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_GMVlat_85_reg_free,fit_sem_GMVlat_85_reg_constrained))

#Data Visualization for the regularized grey matter volume model
BIC(fit_sem_GMVlat_85_reg_free)                      
data_subsample_GMVlat_85_rl_reg_GGSEG <- data.frame(label=c("lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_lateraloccipital","rh_lateraloccipital","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_pericalcarine","rh_pericalcarine","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_supramarginal","rh_supramarginal","lh_insula","rh_insula"),
                                                    GMV_p=rep(NA, times=22),
                                                    GMV_std=rep(NA, times=22))
a=1
for (i in seq(from = 1, to = 22, by = 2)) {
  data_subsample_GMVlat_85_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_GMVlat_85_reg_free)[a+5,7]
  data_subsample_GMVlat_85_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMVlat_85_reg_free)[a+5,7]
  data_subsample_GMVlat_85_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_GMVlat_85_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMVlat_85_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_GMVlat_85_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_GMVlat_85_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMVlat_85_rl_reg_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#Model for grey matter volume with all the regions of interest (34 ROIs)
Model_GMVlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ bankssts_gmv + caudalanteriorcingulate_gmv + caudalmiddlefrontal_gmv + cuneus_gmv + entorhinal_gmv + fusiform_gmv + inferiorparietal_gmv + inferiortemporal_gmv + isthmuscingulate_gmv + lateraloccipital_gmv + lateralorbitofrontal_gmv + lingual_gmv + medialorbitofrontal_gmv + middletemporal_gmv + parahippocampal_gmv + paracentral_gmv + parsopercularis_gmv + parsorbitalis_gmv + parstriangularis_gmv + pericalcarine_gmv + postcentral_gmv + posteriorcingulate_gmv + precentral_gmv + precuneus_gmv + rostralanteriorcingulate_gmv + rostralmiddlefrontal_gmv + superiorfrontal_gmv + superiorparietal_gmv + superiortemporal_gmv + supramarginal_gmv + frontalpole_gmv + temporalpole_gmv + transversetemporal_gmv + insula_gmv'

fit_sem_GMVlat_85_all_free <- sem(Model_GMVlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMVlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the grey matter volume model with all the regions
Model_GMVlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ a*bankssts_gmv + a*caudalanteriorcingulate_gmv + a*caudalmiddlefrontal_gmv + a*cuneus_gmv + a*entorhinal_gmv + a*fusiform_gmv + a*inferiorparietal_gmv + a*inferiortemporal_gmv + a*isthmuscingulate_gmv + a*lateraloccipital_gmv + a*lateralorbitofrontal_gmv + a*lingual_gmv + a*medialorbitofrontal_gmv + a*middletemporal_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*parsopercularis_gmv + a*parsorbitalis_gmv + a*parstriangularis_gmv + a*pericalcarine_gmv + a*postcentral_gmv + a*posteriorcingulate_gmv + a*precentral_gmv + a*precuneus_gmv + a*rostralanteriorcingulate_gmv + a*rostralmiddlefrontal_gmv + a*superiorfrontal_gmv + a*superiorparietal_gmv + a*superiortemporal_gmv + a*supramarginal_gmv + a*frontalpole_gmv + a*temporalpole_gmv + a*transversetemporal_gmv + a*insula_gmv'

fit_sem_GMVlat_85_all_constrained <- sem(Model_GMVlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMVlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_GMVlat_85_all_free,fit_sem_GMVlat_85_all_constrained))

#Data Visualization for the grey matter volume model with all the regions
BIC(fit_sem_GMVlat_85_all_free)
data_subsample_GMVlat_85_rl_all_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                    GMV_p=rep(NA, times=68),
                                                    GMV_std=rep(NA, times=68))
a=1
for (i in seq(from = 1, to = 68, by = 2)) {
  data_subsample_GMVlat_85_rl_all_GGSEG[i,2]<-parameterEstimates(fit_sem_GMVlat_85_all_free)[a+5,7]
  data_subsample_GMVlat_85_rl_all_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMVlat_85_all_free)[a+5,7]
  data_subsample_GMVlat_85_rl_all_GGSEG[i,3]<-lavInspect(fit_sem_GMVlat_85_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMVlat_85_rl_all_GGSEG[i+1,3]<-lavInspect(fit_sem_GMVlat_85_all_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_GMVlat_85_rl_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMVlat_85_rl_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#Create tables to save the parameters of the regularized model and the model with all the ROIs for grey matter volume
paramest_GMlat_GMVall_85 <- data.frame(lhs=rep(NA, times=34),
                                       op=rep(NA, times=34),
                                       rhs=rep(NA, times=34),
                                       est=rep(NA, times=34),
                                       se=rep(NA, times=34),
                                       z=rep(NA, times=34),
                                       pvalue=rep(NA, times=34),
                                       ci.lower=rep(NA, times=34),
                                       ci.upper=rep(NA, times=34),
                                       std.lv=rep(NA, times=34),
                                       std.all=rep(NA, times=34),
                                       std.nox=rep(NA, times=34))

paramest_GMlat_GMVall_85[,] <-parameterEstimates(fit_sem_GMVlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:39,]
write.csv(paramest_GMlat_GMVall_85, "paramest_GMlat_GMVall_85.csv", row.names=FALSE)

paramest_GMlat_GMVreg_85 <- data.frame(lhs=rep(NA, times=11),
                                       op=rep(NA, times=11),
                                       rhs=rep(NA, times=11),
                                       est=rep(NA, times=11),
                                       se=rep(NA, times=11),
                                       z=rep(NA, times=11),
                                       pvalue=rep(NA, times=11),
                                       ci.lower=rep(NA, times=11),
                                       ci.upper=rep(NA, times=11),
                                       std.lv=rep(NA, times=11),
                                       std.all=rep(NA, times=11),
                                       std.nox=rep(NA, times=11))

paramest_GMlat_GMVreg_85[,] <-parameterEstimates(fit_sem_GMVlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:16,]
write.csv(paramest_GMlat_GMVreg_85, "paramest_GMlat_GMVreg_85.csv", row.names=FALSE)
save.image("~/Project 1/Data/Workspace_Script_final_Model.RData")

#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Models including the three grey matter metrics
#-------------------------------------------------------------------------------------------------------------------------------

#Regularized grey matter metrics model 
Model_GMlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + entorhinal_ct + fusiform_ct + inferiorparietal_ct + lateraloccipital_ct + lateralorbitofrontal_ct + lingual_ct + middletemporal_ct + parahippocampal_ct + parsopercularis_ct + pericalcarine_ct + superiorfrontal_ct + superiorparietal_ct + frontalpole_ct + temporalpole_ct + insula_ct + caudalmiddlefrontal_sa + entorhinal_sa + inferiorparietal_sa + isthmuscingulate_sa + lateraloccipital_sa + medialorbitofrontal_sa + middletemporal_sa + paracentral_sa + postcentral_sa + precentral_sa + precuneus_sa + superiorfrontal_sa + superiorparietal_sa + superiortemporal_sa + transversetemporal_sa + insula_sa + bankssts_gmv + cuneus_gmv + parahippocampal_gmv + paracentral_gmv + parstriangularis_gmv + precuneus_gmv + rostralmiddlefrontal_gmv + superiorfrontal_gmv + temporalpole_gmv'

fit_sem_GMlat_85_reg_free <- sem(Model_GMlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized grey matter metrics model
Model_GMlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*entorhinal_ct + a*fusiform_ct + a*inferiorparietal_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*lingual_ct + a*middletemporal_ct + a*parahippocampal_ct + a*parsopercularis_ct + a*pericalcarine_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*frontalpole_ct + a*temporalpole_ct + a*insula_ct + a*caudalmiddlefrontal_sa + a*entorhinal_sa + a*inferiorparietal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*medialorbitofrontal_sa + a*middletemporal_sa + a*paracentral_sa + a*postcentral_sa + a*precentral_sa + a*precuneus_sa + a*superiorfrontal_sa + a*superiorparietal_sa + a*superiortemporal_sa + a*transversetemporal_sa + a*insula_sa + a*bankssts_gmv + a*cuneus_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*parstriangularis_gmv + a*precuneus_gmv + a*rostralmiddlefrontal_gmv + a*superiorfrontal_gmv + a*temporalpole_gmv'

fit_sem_GMlat_85_reg_constrained <- sem(Model_GMlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_GMlat_85_reg_free,fit_sem_GMlat_85_reg_constrained))

#Data Visualization for the regularized grey matter metrics model
BIC(fit_sem_GMlat_85_reg_free)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
data_subsample_GMlat_85_regCT_GGSEG <- data.frame(label=c("lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_parsopercularis","rh_parsopercularis","lh_pericalcarine","rh_pericalcarine","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_insula","rh_insula"),
                                                  CT_p=rep(NA, times=34),
                                                  CT_std=rep(NA, times=34))
a=1
for (i in seq(from = 1, to = 34, by = 2)) {
  data_subsample_GMlat_85_regCT_GGSEG[i,2]<-parameterEstimates(fit_sem_GMlat_85_reg_free)[a+5,7]
  data_subsample_GMlat_85_regCT_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMlat_85_reg_free)[a+5,7]
  data_subsample_GMlat_85_regCT_GGSEG[i,3]<-lavInspect(fit_sem_GMlat_85_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMlat_85_regCT_GGSEG[i+1,3]<-lavInspect(fit_sem_GMlat_85_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

data_subsample_GMlat_85_regSA_GGSEG <- data.frame(label=c("lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_entorhinal","rh_entorhinal","lh_inferiorparietal","rh_inferiorparietal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_paracentral","rh_paracentral","lh_postcentral","rh_postcentral","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                  SA_p=rep(NA, times=32),
                                                  SA_std=rep(NA, times=32))
a=1
for (i in seq(from = 1, to = 32, by = 2)) {
  data_subsample_GMlat_85_regSA_GGSEG[i,2]<-parameterEstimates(fit_sem_GMlat_85_reg_free)[a+19,7]
  data_subsample_GMlat_85_regSA_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMlat_85_reg_free)[a+19,7]
  data_subsample_GMlat_85_regSA_GGSEG[i,3]<-lavInspect(fit_sem_GMlat_85_reg_free,what = "std.all")$beta[1,a+15]
  data_subsample_GMlat_85_regSA_GGSEG[i+1,3]<-lavInspect(fit_sem_GMlat_85_reg_free,what = "std.all")$beta[1,a+15]
  a=a+1
}

data_subsample_GMlat_85_regGMV_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_cuneus","rh_cuneus","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parstriangularis","rh_parstriangularis","lh_precuneus","rh_precuneus","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_temporalpole","rh_temporalpole"),
                                                   GMV_p=rep(NA, times=18),
                                                   GMV_std=rep(NA, times=18))
a=1
for (i in seq(from = 1, to = 18, by = 2)) {
  data_subsample_GMlat_85_regGMV_GGSEG[i,2]<-parameterEstimates(fit_sem_GMlat_85_reg_free)[a+30,7]
  data_subsample_GMlat_85_regGMV_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMlat_85_reg_free)[a+30,7]
  data_subsample_GMlat_85_regGMV_GGSEG[i,3]<-lavInspect(fit_sem_GMlat_85_reg_free,what = "std.all")$beta[1,a+26]
  data_subsample_GMlat_85_regGMV_GGSEG[i+1,3]<-lavInspect(fit_sem_GMlat_85_reg_free,what = "std.all")$beta[1,a+26]
  a=a+1
}

data_subsample_GMlat_85_regCT_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_regCT_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
data_subsample_GMlat_85_regSA_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_regSA_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
data_subsample_GMlat_85_regGMV_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_regGMV_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")


#Compare the information provide by the regularized regions in the three grey matter metrics 
data_subsample_GMlat_85_reg_GGSEG_pvalue<-merge(data_subsample_GMlat_85_regCT_GGSEG[,c(1,2)],data_subsample_GMlat_85_regSA_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_GMlat_85_reg_GGSEG_pvalue<-merge(data_subsample_GMlat_85_reg_GGSEG_pvalue,data_subsample_GMlat_85_regGMV_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_GMlat_85_reg_GGSEG_pvalue, "STD_GMlat_85_reg_pvalue.csv", row.names=FALSE)

data_subsample_GMlat_85_reg_GGSEG_long_pvalue <- data_subsample_GMlat_85_reg_GGSEG_pvalue %>%
  rename(c("CT_p" = "CT")) %>%
  rename(c("SA_p" = "SA")) %>%
  rename(c("GMV_p" = "GMV")) %>%
  gather(Metric, pvalue, CT, SA, GMV)
data_subsample_GMlat_85_reg_GGSEG_long_pvalue$Metric<-factor(data_subsample_GMlat_85_reg_GGSEG_long_pvalue$Metric,levels=c("CT","SA","GMV"))
metric_names <- c(CT = "Cortical Thickness",SA = "Surface Area",GMV = "Grey Matter Volume")

data_subsample_GMlat_85_reg_GGSEG_long_pvalue %>%
  group_by(Metric) %>%
  ggseg(hemisphere="left",colour="black",size=.5,
        mapping=aes(fill=pvalue), show.legend = TRUE) +
  theme_classic(base_size = 20)+
  theme(axis.title.x = element_text(size = 15))+
  scale_fill_gradient(low="goldenrod",high="firebrick") +
  facet_wrap(~Metric, ncol=3,labeller = as_labeller(metric_names))+
  labs(fill="P-value")

data_subsample_GMlat_85_reg_GGSEG_STD<-merge(data_subsample_GMlat_85_regCT_GGSEG[,c(1,3)],data_subsample_GMlat_85_regSA_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_GMlat_85_reg_GGSEG_STD<-merge(data_subsample_GMlat_85_reg_GGSEG_STD,data_subsample_GMlat_85_regGMV_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_GMlat_85_reg_GGSEG_STD, "STD_GMlat_85_reg_STD.csv", row.names=FALSE)

data_subsample_GMlat_85_reg_GGSEG_long_STD <- data_subsample_GMlat_85_reg_GGSEG_STD %>%
  rename(c("CT_std" = "CT")) %>%
  rename(c("SA_std" = "SA")) %>%
  rename(c("GMV_std" = "GMV")) %>%
  gather(Metric, STD, CT, SA, GMV) 
data_subsample_GMlat_85_reg_GGSEG_long_STD$Metric<-factor(data_subsample_GMlat_85_reg_GGSEG_long_STD$Metric,levels=c("CT","SA","GMV"))
metric_names <- c(CT = "Cortical Thickness",SA = "Surface Area",GMV = "Grey Matter Volume")

data_subsample_GMlat_85_reg_GGSEG_long_STD %>%
  group_by(Metric) %>%
  ggseg(hemisphere="left",colour="black",size=.5,
        mapping=aes(fill=STD), show.legend = TRUE) +
  theme_classic(base_size = 20)+
  theme(axis.title.x = element_text(size = 15))+
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick") +
  facet_wrap(~Metric, ncol=3,labeller = as_labeller(metric_names))+
  labs(fill="Standardized\nparameter\nestimates")
ggsave("plot_GMlat_85_regall_STD.png")  


#Model for grey matter metrics with all the regions of interest (102 ROIs)
Model_GMlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ bankssts_ct + caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + cuneus_ct + entorhinal_ct + fusiform_ct + inferiorparietal_ct + inferiortemporal_ct + isthmuscingulate_ct + lateraloccipital_ct + lateralorbitofrontal_ct + lingual_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsopercularis_ct + parsorbitalis_ct + parstriangularis_ct + pericalcarine_ct + postcentral_ct + posteriorcingulate_ct + precentral_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiorparietal_ct + superiortemporal_ct + supramarginal_ct + frontalpole_ct + temporalpole_ct + transversetemporal_ct + insula_ct + bankssts_sa + caudalanteriorcingulate_sa + caudalmiddlefrontal_sa + cuneus_sa + entorhinal_sa + fusiform_sa + inferiorparietal_sa + inferiortemporal_sa + isthmuscingulate_sa + lateraloccipital_sa + lateralorbitofrontal_sa + lingual_sa + medialorbitofrontal_sa + middletemporal_sa + parahippocampal_sa + paracentral_sa + parsopercularis_sa + parsorbitalis_sa + parstriangularis_sa + pericalcarine_sa + postcentral_sa + posteriorcingulate_sa + precentral_sa + precuneus_sa + rostralanteriorcingulate_sa + rostralmiddlefrontal_sa + superiorfrontal_sa + superiorparietal_sa + superiortemporal_sa + supramarginal_sa + frontalpole_sa + temporalpole_sa + transversetemporal_sa + insula_sa + bankssts_gmv + caudalanteriorcingulate_gmv + caudalmiddlefrontal_gmv + cuneus_gmv + entorhinal_gmv + fusiform_gmv + inferiorparietal_gmv + inferiortemporal_gmv + isthmuscingulate_gmv + lateraloccipital_gmv + lateralorbitofrontal_gmv + lingual_gmv + medialorbitofrontal_gmv + middletemporal_gmv + parahippocampal_gmv + paracentral_gmv + parsopercularis_gmv + parsorbitalis_gmv + parstriangularis_gmv + pericalcarine_gmv + postcentral_gmv + posteriorcingulate_gmv + precentral_gmv + precuneus_gmv + rostralanteriorcingulate_gmv + rostralmiddlefrontal_gmv + superiorfrontal_gmv + superiorparietal_gmv + superiortemporal_gmv + supramarginal_gmv + frontalpole_gmv + temporalpole_gmv + transversetemporal_gmv + insula_gmv'

fit_sem_GMlat_85_all_free <- sem(Model_GMlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the grey matter metrics model with all the regions
Model_GMlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*bankssts_ct + a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*cuneus_ct + a*entorhinal_ct + a*fusiform_ct + a*inferiorparietal_ct + a*inferiortemporal_ct + a*isthmuscingulate_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*lingual_ct + a*medialorbitofrontal_ct + a*middletemporal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsopercularis_ct + a*parsorbitalis_ct + a*parstriangularis_ct + a*pericalcarine_ct + a*postcentral_ct + a*posteriorcingulate_ct + a*precentral_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*superiortemporal_ct + a*supramarginal_ct + a*frontalpole_ct + a*temporalpole_ct + a*transversetemporal_ct + a*insula_ct + a*bankssts_sa + a*caudalanteriorcingulate_sa + a*caudalmiddlefrontal_sa + a*cuneus_sa + a*entorhinal_sa + a*fusiform_sa + a*inferiorparietal_sa + a*inferiortemporal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*lateralorbitofrontal_sa + a*lingual_sa + a*medialorbitofrontal_sa + a*middletemporal_sa + a*parahippocampal_sa + a*paracentral_sa + a*parsopercularis_sa + a*parsorbitalis_sa + a*parstriangularis_sa + a*pericalcarine_sa + a*postcentral_sa + a*posteriorcingulate_sa + a*precentral_sa + a*precuneus_sa + a*rostralanteriorcingulate_sa + a*rostralmiddlefrontal_sa + a*superiorfrontal_sa + a*superiorparietal_sa + a*superiortemporal_sa + a*supramarginal_sa + a*frontalpole_sa + a*temporalpole_sa + a*transversetemporal_sa + a*insula_sa + a*bankssts_gmv + a*caudalanteriorcingulate_gmv + a*caudalmiddlefrontal_gmv + a*cuneus_gmv + a*entorhinal_gmv + a*fusiform_gmv + a*inferiorparietal_gmv + a*inferiortemporal_gmv + a*isthmuscingulate_gmv + a*lateraloccipital_gmv + a*lateralorbitofrontal_gmv + a*lingual_gmv + a*medialorbitofrontal_gmv + a*middletemporal_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*parsopercularis_gmv + a*parsorbitalis_gmv + a*parstriangularis_gmv + a*pericalcarine_gmv + a*postcentral_gmv + a*posteriorcingulate_gmv + a*precentral_gmv + a*precuneus_gmv + a*rostralanteriorcingulate_gmv + a*rostralmiddlefrontal_gmv + a*superiorfrontal_gmv + a*superiorparietal_gmv + a*superiortemporal_gmv + a*supramarginal_gmv + a*frontalpole_gmv + a*temporalpole_gmv + a*transversetemporal_gmv + a*insula_gmv'

fit_sem_GMlat_85_all_constrained <- sem(Model_GMlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_GMlat_85_all_free,fit_sem_GMlat_85_all_constrained))

#Data Visualization for the grey matter volume model with all the regions
BIC(fit_sem_GMlat_85_all_free)
data_subsample_GMlat_85_all_GGSEG <- data.frame(label=c("lh_bankssts","rh_bankssts","lh_caudalanteriorcingulate","rh_caudalanteriorcingulate","lh_caudalmiddlefrontal","rh_caudalmiddlefrontal","lh_cuneus","rh_cuneus","lh_entorhinal","rh_entorhinal","lh_fusiform","rh_fusiform","lh_inferiorparietal","rh_inferiorparietal","lh_inferiortemporal","rh_inferiortemporal","lh_isthmuscingulate","rh_isthmuscingulate","lh_lateraloccipital","rh_lateraloccipital","lh_lateralorbitofrontal","rh_lateralorbitofrontal","lh_lingual","rh_lingual","lh_medialorbitofrontal","rh_medialorbitofrontal","lh_middletemporal","rh_middletemporal","lh_parahippocampal","rh_parahippocampal","lh_paracentral","rh_paracentral","lh_parsopercularis","rh_parsopercularis","lh_parsorbitalis","rh_parsorbitalis","lh_parstriangularis","rh_parstriangularis","lh_pericalcarine","rh_pericalcarine","lh_postcentral","rh_postcentral","lh_posteriorcingulate","rh_posteriorcingulate","lh_precentral","rh_precentral","lh_precuneus","rh_precuneus","lh_rostralanteriorcingulate","rh_rostralanteriorcingulate","lh_rostralmiddlefrontal","rh_rostralmiddlefrontal","lh_superiorfrontal","rh_superiorfrontal","lh_superiorparietal","rh_superiorparietal","lh_superiortemporal","rh_superiortemporal","lh_supramarginal","rh_supramarginal","lh_frontalpole","rh_frontalpole","lh_temporalpole","rh_temporalpole","lh_transversetemporal","rh_transversetemporal","lh_insula","rh_insula"),
                                                   CT_p=rep(NA, times=68),
                                                   CT_std=rep(NA, times=68),
                                                   SA_p=rep(NA, times=68),
                                                   SA_std=rep(NA, times=68),
                                                   GMV_p=rep(NA, times=68),
                                                   GMV_std=rep(NA, times=68))

a=1
for (i in seq(from = 1, to = 68, by = 2)) {
  data_subsample_GMlat_85_all_GGSEG[i,2]<-parameterEstimates(fit_sem_GMlat_85_all_free)[a+5,7]
  data_subsample_GMlat_85_all_GGSEG[i+1,2]<-parameterEstimates(fit_sem_GMlat_85_all_free)[a+5,7]
  data_subsample_GMlat_85_all_GGSEG[i,4]<-parameterEstimates(fit_sem_GMlat_85_all_free)[a+39,7]
  data_subsample_GMlat_85_all_GGSEG[i+1,4]<-parameterEstimates(fit_sem_GMlat_85_all_free)[a+39,7]
  data_subsample_GMlat_85_all_GGSEG[i,6]<-parameterEstimates(fit_sem_GMlat_85_all_free)[a+73,7]
  data_subsample_GMlat_85_all_GGSEG[i+1,6]<-parameterEstimates(fit_sem_GMlat_85_all_free)[a+73,7]
  data_subsample_GMlat_85_all_GGSEG[i,3]<-lavInspect(fit_sem_GMlat_85_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMlat_85_all_GGSEG[i+1,3]<-lavInspect(fit_sem_GMlat_85_all_free,what = "std.all")$beta[1,a+1]
  data_subsample_GMlat_85_all_GGSEG[i,5]<-lavInspect(fit_sem_GMlat_85_all_free,what = "std.all")$beta[1,a+35]
  data_subsample_GMlat_85_all_GGSEG[i+1,5]<-lavInspect(fit_sem_GMlat_85_all_free,what = "std.all")$beta[1,a+35]
  data_subsample_GMlat_85_all_GGSEG[i,7]<-lavInspect(fit_sem_GMlat_85_all_free,what = "std.all")$beta[1,a+69]
  data_subsample_GMlat_85_all_GGSEG[i+1,7]<-lavInspect(fit_sem_GMlat_85_all_free,what = "std.all")$beta[1,a+69]
  a=a+1
}

data_subsample_GMlat_85_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(CT_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(CT_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
data_subsample_GMlat_85_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(SA_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(SA_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")
data_subsample_GMlat_85_all_GGSEG %>% 
  ggseg(mapping=aes(fill=as.numeric(GMV_p)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient(low="goldenrod",high="firebrick",name="p-value")
data_subsample_GMlat_85_all_GGSEG %>%
  ggseg(mapping=aes(fill=as.numeric(GMV_std)), position="stacked",colour="black",size=.5) +
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick",name="Standardized\nparameter\nestimates")

data_subsample_GMlat_85_all_GGSEG_long_STD <- data_subsample_GMlat_85_all_GGSEG %>%
  rename(c("CT_std" = "CT")) %>%
  rename(c("SA_std" = "SA")) %>%
  rename(c("GMV_std" = "GMV")) %>%
  gather(Metric, STD, CT, SA, GMV)
data_subsample_GMlat_85_all_GGSEG_long_STD$Metric<-factor(data_subsample_GMlat_85_all_GGSEG_long_STD$Metric,levels=c("CT","SA","GMV"))

data_subsample_GMlat_85_all_GGSEG_long_pvalue <- data_subsample_GMlat_85_all_GGSEG %>%
  rename(c("CT_p" = "CT")) %>%
  rename(c("SA_p" = "SA")) %>%
  rename(c("GMV_p" = "GMV")) %>%
  gather(Metric, pvalue, CT, SA, GMV)
data_subsample_GMlat_85_all_GGSEG_long_pvalue$Metric<-factor(data_subsample_GMlat_85_all_GGSEG_long_pvalue$Metric,levels=c("CT","SA","GMV"))
metric_names <- c(CT = "Cortical Thickness",SA = "Surface Area",GMV = "Grey Matter Volume")

data_subsample_GMlat_85_all_GGSEG_long_pvalue %>%
  group_by(Metric) %>%
  ggseg(hemisphere="left",colour="black",size=.5,
        mapping=aes(fill=pvalue), show.legend = TRUE) +
  theme_classic(base_size = 20)+
  theme(axis.title.x = element_text(size = 15))+
  scale_fill_gradient(low="goldenrod",high="firebrick") +
  facet_wrap(~Metric, ncol=3,labeller = as_labeller(metric_names))+
  labs(fill="P-value")

data_subsample_GMlat_85_all_GGSEG_long_STD %>%
  group_by(Metric) %>%
  ggseg(hemisphere="left",colour="black",size=.5,
        mapping=aes(fill=STD), show.legend = TRUE) +
  theme_classic(base_size = 20)+
  theme(axis.title.x = element_text(size = 15))+
  scale_fill_gradient2(midpoint=0,low="goldenrod",mid="white",high="firebrick") +
  facet_wrap(~Metric, ncol=3,labeller = as_labeller(metric_names))+
  labs(fill="Standardized\nparameter\nestimates")


#Create tables to save the parameters of the regularized model and the model with all the ROIs for grey matter metrics
paramest_GMlat_GMall_85 <- data.frame(lhs=rep(NA, times=102),
                                      op=rep(NA, times=102),
                                      rhs=rep(NA, times=102),
                                      est=rep(NA, times=102),
                                      se=rep(NA, times=102),
                                      z=rep(NA, times=102),
                                      pvalue=rep(NA, times=102),
                                      ci.lower=rep(NA, times=102),
                                      ci.upper=rep(NA, times=102),
                                      std.lv=rep(NA, times=102),
                                      std.all=rep(NA, times=102),
                                      std.nox=rep(NA, times=102))

paramest_GMlat_GMall_85[,] <-parameterEstimates(fit_sem_GMlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:107,]
write.csv(paramest_GMlat_GMall_85, "paramest_GMlat_GMall_85.csv", row.names=FALSE)

paramest_GMlat_GMreg_85 <- data.frame(lhs=rep(NA, times=42),
                                      op=rep(NA, times=42),
                                      rhs=rep(NA, times=42),
                                      est=rep(NA, times=42),
                                      se=rep(NA, times=42),
                                      z=rep(NA, times=42),
                                      pvalue=rep(NA, times=42),
                                      ci.lower=rep(NA, times=42),
                                      ci.upper=rep(NA, times=42),
                                      std.lv=rep(NA, times=42),
                                      std.all=rep(NA, times=42),
                                      std.nox=rep(NA, times=42))

paramest_GMlat_GMreg_85[,] <-parameterEstimates(fit_sem_GMlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:47,]
write.csv(paramest_GMlat_GMreg_85, "paramest_GMlat_GMreg_85.csv", row.names=FALSE)


#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Individual models white matter
#-------------------------------------------------------------------------------------------------------------------------------
  
#Create a table to save the regression coefficients or the fit of the ROI model for each ROIs
data_subsample_WMlat_85_GGSEG <- data.frame(label=c("fornix","cingulatecingulum","parahippocampalcingulum","corticospinalpyramidal","anteriorthalamicradiations","uncinate","inferiorlongitudinalfasiculus","inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","corpuscallosum","superiorlongitudinalfasiculus","temporalsuperiorlongitudinalfasiculus","parietalsuperiorlongitudinalfasiculus","superiorcorticostriate","superiorcorticostriatefrontalcortex","superiorcorticostriateparietalcortex","striatalinferiorfrontalcortex","inferiorfrontalsuperiorfrontalcortex","fornix_exfimbria"),
                                            FA_fit=rep(NA, times=20),
                                            FA_p=rep(NA, times=20),
                                            FA_std=rep(NA, times=20),
                                            MD_fit=rep(NA, times=20),
                                            MD_p=rep(NA, times=20),
                                            MD_std=rep(NA, times=20),
                                            WMV_fit=rep(NA, times=20),
                                            WMV_p=rep(NA, times=20),
                                            WMV_std=rep(NA, times=20))

#Create tables to save the parameters of the models for each ROI in each metric
paramest_WMlat_FA_85 <- data.frame(lhs=rep(NA, times=20),
                                   op=rep(NA, times=20),
                                   rhs=rep(NA, times=20),
                                   est=rep(NA, times=20),
                                   se=rep(NA, times=20),
                                   z=rep(NA, times=20),
                                   pvalue=rep(NA, times=20),
                                   ci.lower=rep(NA, times=20),
                                   ci.upper=rep(NA, times=20),
                                   std.lv=rep(NA, times=20),
                                   std.all=rep(NA, times=20),
                                   std.nox=rep(NA, times=20))
paramest_WMlat_MD_85 <- data.frame(lhs=rep(NA, times=20),
                                   op=rep(NA, times=20),
                                   rhs=rep(NA, times=20),
                                   est=rep(NA, times=20),
                                   se=rep(NA, times=20),
                                   z=rep(NA, times=20),
                                   pvalue=rep(NA, times=20),
                                   ci.lower=rep(NA, times=20),
                                   ci.upper=rep(NA, times=20),
                                   std.lv=rep(NA, times=20),
                                   std.all=rep(NA, times=20),
                                   std.nox=rep(NA, times=20))
paramest_WMlat_WMV_85 <- data.frame(lhs=rep(NA, times=20),
                                    op=rep(NA, times=20),
                                    rhs=rep(NA, times=20),
                                    est=rep(NA, times=20),
                                    se=rep(NA, times=20),
                                    z=rep(NA, times=20),
                                    pvalue=rep(NA, times=20),
                                    ci.lower=rep(NA, times=20),
                                    ci.upper=rep(NA, times=20),
                                    std.lv=rep(NA, times=20),
                                    std.all=rep(NA, times=20),
                                    std.nox=rep(NA, times=20))

#Fractional Anisotropy (20 regions)
##fornix
Model_FA_fornix <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fornix_fa'
fit_sem_FA_fornix <- sem(Model_FA_fornix, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_fornix, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[1,2]<-BIC(fit_sem_FA_fornix)
data_subsample_WMlat_85_GGSEG[1,3]<-parameterEstimates(fit_sem_FA_fornix)[6,7]
data_subsample_WMlat_85_GGSEG[1,4]<-lavInspect(fit_sem_FA_fornix,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[1,] <-parameterEstimates(fit_sem_FA_fornix,standardized = TRUE, rsquare=TRUE)[6,]

##cingulatecingulum
Model_FA_cingulatecingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ cingulatecingulum_fa'
fit_sem_FA_cingulatecingulum <- sem(Model_FA_cingulatecingulum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_cingulatecingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[2,2]<-BIC(fit_sem_FA_cingulatecingulum)
data_subsample_WMlat_85_GGSEG[2,3]<-parameterEstimates(fit_sem_FA_cingulatecingulum)[6,7]
data_subsample_WMlat_85_GGSEG[2,4]<-lavInspect(fit_sem_FA_cingulatecingulum,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[2,] <-parameterEstimates(fit_sem_FA_cingulatecingulum,standardized = TRUE, rsquare=TRUE)[6,]

##parahippocampalcingulum
Model_FA_parahippocampalcingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parahippocampalcingulum_fa'
fit_sem_FA_parahippocampalcingulum <- sem(Model_FA_parahippocampalcingulum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_parahippocampalcingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[3,2]<-BIC(fit_sem_FA_parahippocampalcingulum)
data_subsample_WMlat_85_GGSEG[3,3]<-parameterEstimates(fit_sem_FA_parahippocampalcingulum)[6,7]
data_subsample_WMlat_85_GGSEG[3,4]<-lavInspect(fit_sem_FA_parahippocampalcingulum,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[3,] <-parameterEstimates(fit_sem_FA_parahippocampalcingulum,standardized = TRUE, rsquare=TRUE)[6,]

##corticospinalpyramidal
Model_FA_corticospinalpyramidal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ corticospinalpyramidal_fa'
fit_sem_FA_corticospinalpyramidal <- sem(Model_FA_corticospinalpyramidal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_corticospinalpyramidal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[4,2]<-BIC(fit_sem_FA_corticospinalpyramidal)
data_subsample_WMlat_85_GGSEG[4,3]<-parameterEstimates(fit_sem_FA_corticospinalpyramidal)[6,7]
data_subsample_WMlat_85_GGSEG[4,4]<-lavInspect(fit_sem_FA_corticospinalpyramidal,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[4,] <-parameterEstimates(fit_sem_FA_corticospinalpyramidal,standardized = TRUE, rsquare=TRUE)[6,]

##anteriorthalamicradiations
Model_FA_anteriorthalamicradiations <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ anteriorthalamicradiations_fa'
fit_sem_FA_anteriorthalamicradiations <- sem(Model_FA_anteriorthalamicradiations, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_anteriorthalamicradiations, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[5,2]<-BIC(fit_sem_FA_anteriorthalamicradiations)
data_subsample_WMlat_85_GGSEG[5,3]<-parameterEstimates(fit_sem_FA_anteriorthalamicradiations)[6,7]
data_subsample_WMlat_85_GGSEG[5,4]<-lavInspect(fit_sem_FA_anteriorthalamicradiations,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[5,] <-parameterEstimates(fit_sem_FA_anteriorthalamicradiations,standardized = TRUE, rsquare=TRUE)[6,]

##uncinate
Model_FA_uncinate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ uncinate_fa'
fit_sem_FA_uncinate <- sem(Model_FA_uncinate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_uncinate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[6,2]<-BIC(fit_sem_FA_uncinate)
data_subsample_WMlat_85_GGSEG[6,3]<-parameterEstimates(fit_sem_FA_uncinate)[6,7]
data_subsample_WMlat_85_GGSEG[6,4]<-lavInspect(fit_sem_FA_uncinate,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[6,] <-parameterEstimates(fit_sem_FA_uncinate,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorlongitudinalfasiculus
Model_FA_inferiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorlongitudinalfasiculus_fa'
fit_sem_FA_inferiorlongitudinalfasiculus <- sem(Model_FA_inferiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_inferiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[7,2]<-BIC(fit_sem_FA_inferiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[7,3]<-parameterEstimates(fit_sem_FA_inferiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[7,4]<-lavInspect(fit_sem_FA_inferiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[7,] <-parameterEstimates(fit_sem_FA_inferiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorfrontooccipitalfasiculus
Model_FA_inferiorfrontooccipitalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorfrontooccipitalfasiculus_fa'
fit_sem_FA_inferiorfrontooccipitalfasiculus <- sem(Model_FA_inferiorfrontooccipitalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_inferiorfrontooccipitalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[8,2]<-BIC(fit_sem_FA_inferiorfrontooccipitalfasiculus)
data_subsample_WMlat_85_GGSEG[8,3]<-parameterEstimates(fit_sem_FA_inferiorfrontooccipitalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[8,4]<-lavInspect(fit_sem_FA_inferiorfrontooccipitalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[8,] <-parameterEstimates(fit_sem_FA_inferiorfrontooccipitalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##forcepsmajor
Model_FA_forcepsmajor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ forcepsmajor_fa'
fit_sem_FA_forcepsmajor <- sem(Model_FA_forcepsmajor, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_forcepsmajor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[9,2]<-BIC(fit_sem_FA_forcepsmajor)
data_subsample_WMlat_85_GGSEG[9,3]<-parameterEstimates(fit_sem_FA_forcepsmajor)[6,7]
data_subsample_WMlat_85_GGSEG[9,4]<-lavInspect(fit_sem_FA_forcepsmajor,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[9,] <-parameterEstimates(fit_sem_FA_forcepsmajor,standardized = TRUE, rsquare=TRUE)[6,]

##forcepsminor
Model_FA_forcepsminor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ forcepsminor_fa'
fit_sem_FA_forcepsminor <- sem(Model_FA_forcepsminor, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_forcepsminor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[10,2]<-BIC(fit_sem_FA_forcepsminor)
data_subsample_WMlat_85_GGSEG[10,3]<-parameterEstimates(fit_sem_FA_forcepsminor)[6,7]
data_subsample_WMlat_85_GGSEG[10,4]<-lavInspect(fit_sem_FA_forcepsminor,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[10,] <-parameterEstimates(fit_sem_FA_forcepsminor,standardized = TRUE, rsquare=TRUE)[6,]

##corpuscallosum
Model_FA_corpuscallosum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ corpuscallosum_fa'
fit_sem_FA_corpuscallosum <- sem(Model_FA_corpuscallosum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_corpuscallosum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[11,2]<-BIC(fit_sem_FA_corpuscallosum)
data_subsample_WMlat_85_GGSEG[11,3]<-parameterEstimates(fit_sem_FA_corpuscallosum)[6,7]
data_subsample_WMlat_85_GGSEG[11,4]<-lavInspect(fit_sem_FA_corpuscallosum,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[11,] <-parameterEstimates(fit_sem_FA_corpuscallosum,standardized = TRUE, rsquare=TRUE)[6,]

##superiorlongitudinalfasiculus
Model_FA_superiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorlongitudinalfasiculus_fa'
fit_sem_FA_superiorlongitudinalfasiculus <- sem(Model_FA_superiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_superiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[12,2]<-BIC(fit_sem_FA_superiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[12,3]<-parameterEstimates(fit_sem_FA_superiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[12,4]<-lavInspect(fit_sem_FA_superiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[12,] <-parameterEstimates(fit_sem_FA_superiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##temporalsuperiorlongitudinalfasiculus
Model_FA_temporalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalsuperiorlongitudinalfasiculus_fa'
fit_sem_FA_temporalsuperiorlongitudinalfasiculus <- sem(Model_FA_temporalsuperiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_temporalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[13,2]<-BIC(fit_sem_FA_temporalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[13,3]<-parameterEstimates(fit_sem_FA_temporalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[13,4]<-lavInspect(fit_sem_FA_temporalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[13,] <-parameterEstimates(fit_sem_FA_temporalsuperiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##parietalsuperiorlongitudinalfasiculus
Model_FA_parietalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parietalsuperiorlongitudinalfasiculus_fa'
fit_sem_FA_parietalsuperiorlongitudinalfasiculus <- sem(Model_FA_parietalsuperiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_parietalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[14,2]<-BIC(fit_sem_FA_parietalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[14,3]<-parameterEstimates(fit_sem_FA_parietalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[14,4]<-lavInspect(fit_sem_FA_parietalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[14,] <-parameterEstimates(fit_sem_FA_parietalsuperiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriate
Model_FA_superiorcorticostriate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriate_fa'
fit_sem_FA_superiorcorticostriate <- sem(Model_FA_superiorcorticostriate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_superiorcorticostriate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[15,2]<-BIC(fit_sem_FA_superiorcorticostriate)
data_subsample_WMlat_85_GGSEG[15,3]<-parameterEstimates(fit_sem_FA_superiorcorticostriate)[6,7]
data_subsample_WMlat_85_GGSEG[15,4]<-lavInspect(fit_sem_FA_superiorcorticostriate,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[15,] <-parameterEstimates(fit_sem_FA_superiorcorticostriate,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriatefrontalcortex
Model_FA_superiorcorticostriatefrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriatefrontalcortex_fa'
fit_sem_FA_superiorcorticostriatefrontalcortex <- sem(Model_FA_superiorcorticostriatefrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_superiorcorticostriatefrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[16,2]<-BIC(fit_sem_FA_superiorcorticostriatefrontalcortex)
data_subsample_WMlat_85_GGSEG[16,3]<-parameterEstimates(fit_sem_FA_superiorcorticostriatefrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[16,4]<-lavInspect(fit_sem_FA_superiorcorticostriatefrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[16,] <-parameterEstimates(fit_sem_FA_superiorcorticostriatefrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriateparietalcortex
Model_FA_superiorcorticostriateparietalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriateparietalcortex_fa'
fit_sem_FA_superiorcorticostriateparietalcortex <- sem(Model_FA_superiorcorticostriateparietalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_superiorcorticostriateparietalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[17,2]<-BIC(fit_sem_FA_superiorcorticostriateparietalcortex)
data_subsample_WMlat_85_GGSEG[17,3]<-parameterEstimates(fit_sem_FA_superiorcorticostriateparietalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[17,4]<-lavInspect(fit_sem_FA_superiorcorticostriateparietalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[17,] <-parameterEstimates(fit_sem_FA_superiorcorticostriateparietalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##striatalinferiorfrontalcortex
Model_FA_striatalinferiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ striatalinferiorfrontalcortex_fa'
fit_sem_FA_striatalinferiorfrontalcortex <- sem(Model_FA_striatalinferiorfrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_striatalinferiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[18,2]<-BIC(fit_sem_FA_striatalinferiorfrontalcortex)
data_subsample_WMlat_85_GGSEG[18,3]<-parameterEstimates(fit_sem_FA_striatalinferiorfrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[18,4]<-lavInspect(fit_sem_FA_striatalinferiorfrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[18,] <-parameterEstimates(fit_sem_FA_striatalinferiorfrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorfrontalsuperiorfrontalcortex
Model_FA_inferiorfrontalsuperiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorfrontalsuperiorfrontalcortex_fa'
fit_sem_FA_inferiorfrontalsuperiorfrontalcortex <- sem(Model_FA_inferiorfrontalsuperiorfrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[19,2]<-BIC(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex)
data_subsample_WMlat_85_GGSEG[19,3]<-parameterEstimates(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[19,4]<-lavInspect(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[19,] <-parameterEstimates(fit_sem_FA_inferiorfrontalsuperiorfrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##fornix_exfimbria
Model_FA_fornix_exfimbria <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fornix_exfimbria_fa'
fit_sem_FA_fornix_exfimbria <- sem(Model_FA_fornix_exfimbria, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FA_fornix_exfimbria, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[20,2]<-BIC(fit_sem_FA_fornix_exfimbria)
data_subsample_WMlat_85_GGSEG[20,3]<-parameterEstimates(fit_sem_FA_fornix_exfimbria)[6,7]
data_subsample_WMlat_85_GGSEG[20,4]<-lavInspect(fit_sem_FA_fornix_exfimbria,what = "std.all")$beta[1,2]
paramest_WMlat_FA_85[20,] <-parameterEstimates(fit_sem_FA_fornix_exfimbria,standardized = TRUE, rsquare=TRUE)[6,]

write.csv(paramest_WMlat_FA_85, "paramest_WMlat_FA_85.csv", row.names=FALSE) 

#Mean Diffusivity (20 regions)
##fornix
Model_MD_fornix <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fornix_md'
fit_sem_MD_fornix <- sem(Model_MD_fornix, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_fornix, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[1,5]<-BIC(fit_sem_MD_fornix)
data_subsample_WMlat_85_GGSEG[1,6]<-parameterEstimates(fit_sem_MD_fornix)[6,7]
data_subsample_WMlat_85_GGSEG[1,7]<-lavInspect(fit_sem_MD_fornix,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[1,] <-parameterEstimates(fit_sem_MD_fornix,standardized = TRUE, rsquare=TRUE)[6,]

##cingulatecingulum
Model_MD_cingulatecingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ cingulatecingulum_md'
fit_sem_MD_cingulatecingulum <- sem(Model_MD_cingulatecingulum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_cingulatecingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[2,5]<-BIC(fit_sem_MD_cingulatecingulum)
data_subsample_WMlat_85_GGSEG[2,6]<-parameterEstimates(fit_sem_MD_cingulatecingulum)[6,7]
data_subsample_WMlat_85_GGSEG[2,7]<-lavInspect(fit_sem_MD_cingulatecingulum,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[2,] <-parameterEstimates(fit_sem_MD_cingulatecingulum,standardized = TRUE, rsquare=TRUE)[6,]

##parahippocampalcingulum
Model_MD_parahippocampalcingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parahippocampalcingulum_md'
fit_sem_MD_parahippocampalcingulum <- sem(Model_MD_parahippocampalcingulum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_parahippocampalcingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[3,5]<-BIC(fit_sem_MD_parahippocampalcingulum)
data_subsample_WMlat_85_GGSEG[3,6]<-parameterEstimates(fit_sem_MD_parahippocampalcingulum)[6,7]
data_subsample_WMlat_85_GGSEG[3,7]<-lavInspect(fit_sem_MD_parahippocampalcingulum,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[3,] <-parameterEstimates(fit_sem_MD_parahippocampalcingulum,standardized = TRUE, rsquare=TRUE)[6,]

##corticospinalpyramidal
Model_MD_corticospinalpyramidal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ corticospinalpyramidal_md'
fit_sem_MD_corticospinalpyramidal <- sem(Model_MD_corticospinalpyramidal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_corticospinalpyramidal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[4,5]<-BIC(fit_sem_MD_corticospinalpyramidal)
data_subsample_WMlat_85_GGSEG[4,6]<-parameterEstimates(fit_sem_MD_corticospinalpyramidal)[6,7]
data_subsample_WMlat_85_GGSEG[4,7]<-lavInspect(fit_sem_MD_corticospinalpyramidal,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[4,] <-parameterEstimates(fit_sem_MD_corticospinalpyramidal,standardized = TRUE, rsquare=TRUE)[6,]

##anteriorthalamicradiations
Model_MD_anteriorthalamicradiations <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ anteriorthalamicradiations_md'
fit_sem_MD_anteriorthalamicradiations <- sem(Model_MD_anteriorthalamicradiations, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_anteriorthalamicradiations, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[5,5]<-BIC(fit_sem_MD_anteriorthalamicradiations)
data_subsample_WMlat_85_GGSEG[5,6]<-parameterEstimates(fit_sem_MD_anteriorthalamicradiations)[6,7]
data_subsample_WMlat_85_GGSEG[5,7]<-lavInspect(fit_sem_MD_anteriorthalamicradiations,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[5,] <-parameterEstimates(fit_sem_MD_anteriorthalamicradiations,standardized = TRUE, rsquare=TRUE)[6,]

##uncinate
Model_MD_uncinate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ uncinate_md'
fit_sem_MD_uncinate <- sem(Model_MD_uncinate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_uncinate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[6,5]<-BIC(fit_sem_MD_uncinate)
data_subsample_WMlat_85_GGSEG[6,6]<-parameterEstimates(fit_sem_MD_uncinate)[6,7]
data_subsample_WMlat_85_GGSEG[6,7]<-lavInspect(fit_sem_MD_uncinate,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[6,] <-parameterEstimates(fit_sem_MD_uncinate,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorlongitudinalfasiculus
Model_MD_inferiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorlongitudinalfasiculus_md'
fit_sem_MD_inferiorlongitudinalfasiculus <- sem(Model_MD_inferiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_inferiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[7,5]<-BIC(fit_sem_MD_inferiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[7,6]<-parameterEstimates(fit_sem_MD_inferiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[7,7]<-lavInspect(fit_sem_MD_inferiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[7,] <-parameterEstimates(fit_sem_MD_inferiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorfrontooccipitalfasiculus
Model_MD_inferiorfrontooccipitalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorfrontooccipitalfasiculus_md'
fit_sem_MD_inferiorfrontooccipitalfasiculus <- sem(Model_MD_inferiorfrontooccipitalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_inferiorfrontooccipitalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[8,5]<-BIC(fit_sem_MD_inferiorfrontooccipitalfasiculus)
data_subsample_WMlat_85_GGSEG[8,6]<-parameterEstimates(fit_sem_MD_inferiorfrontooccipitalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[8,7]<-lavInspect(fit_sem_MD_inferiorfrontooccipitalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[8,] <-parameterEstimates(fit_sem_MD_inferiorfrontooccipitalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##forcepsmajor
Model_MD_forcepsmajor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ forcepsmajor_md'
fit_sem_MD_forcepsmajor <- sem(Model_MD_forcepsmajor, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_forcepsmajor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[9,5]<-BIC(fit_sem_MD_forcepsmajor)
data_subsample_WMlat_85_GGSEG[9,6]<-parameterEstimates(fit_sem_MD_forcepsmajor)[6,7]
data_subsample_WMlat_85_GGSEG[9,7]<-lavInspect(fit_sem_MD_forcepsmajor,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[9,] <-parameterEstimates(fit_sem_MD_forcepsmajor,standardized = TRUE, rsquare=TRUE)[6,]

##forcepsminor
Model_MD_forcepsminor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ forcepsminor_md'
fit_sem_MD_forcepsminor <- sem(Model_MD_forcepsminor, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_forcepsminor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[10,5]<-BIC(fit_sem_MD_forcepsminor)
data_subsample_WMlat_85_GGSEG[10,6]<-parameterEstimates(fit_sem_MD_forcepsminor)[6,7]
data_subsample_WMlat_85_GGSEG[10,7]<-lavInspect(fit_sem_MD_forcepsminor,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[10,] <-parameterEstimates(fit_sem_MD_forcepsminor,standardized = TRUE, rsquare=TRUE)[6,]

##corpuscallosum
Model_MD_corpuscallosum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ corpuscallosum_md'
fit_sem_MD_corpuscallosum <- sem(Model_MD_corpuscallosum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_corpuscallosum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[11,5]<-BIC(fit_sem_MD_corpuscallosum)
data_subsample_WMlat_85_GGSEG[11,6]<-parameterEstimates(fit_sem_MD_corpuscallosum)[6,7]
data_subsample_WMlat_85_GGSEG[11,7]<-lavInspect(fit_sem_MD_corpuscallosum,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[11,] <-parameterEstimates(fit_sem_MD_corpuscallosum,standardized = TRUE, rsquare=TRUE)[6,]

##superiorlongitudinalfasiculus
Model_MD_superiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorlongitudinalfasiculus_md'
fit_sem_MD_superiorlongitudinalfasiculus <- sem(Model_MD_superiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_superiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[12,5]<-BIC(fit_sem_MD_superiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[12,6]<-parameterEstimates(fit_sem_MD_superiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[12,7]<-lavInspect(fit_sem_MD_superiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[12,] <-parameterEstimates(fit_sem_MD_superiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##temporalsuperiorlongitudinalfasiculus
Model_MD_temporalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalsuperiorlongitudinalfasiculus_md'
fit_sem_MD_temporalsuperiorlongitudinalfasiculus <- sem(Model_MD_temporalsuperiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_temporalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[13,5]<-BIC(fit_sem_MD_temporalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[13,6]<-parameterEstimates(fit_sem_MD_temporalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[13,7]<-lavInspect(fit_sem_MD_temporalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[13,] <-parameterEstimates(fit_sem_MD_temporalsuperiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##parietalsuperiorlongitudinalfasiculus
Model_MD_parietalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parietalsuperiorlongitudinalfasiculus_md'
fit_sem_MD_parietalsuperiorlongitudinalfasiculus <- sem(Model_MD_parietalsuperiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_parietalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[14,5]<-BIC(fit_sem_MD_parietalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[14,6]<-parameterEstimates(fit_sem_MD_parietalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[14,7]<-lavInspect(fit_sem_MD_parietalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[14,] <-parameterEstimates(fit_sem_MD_parietalsuperiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriate
Model_MD_superiorcorticostriate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriate_md'
fit_sem_MD_superiorcorticostriate <- sem(Model_MD_superiorcorticostriate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_superiorcorticostriate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[15,5]<-BIC(fit_sem_MD_superiorcorticostriate)
data_subsample_WMlat_85_GGSEG[15,6]<-parameterEstimates(fit_sem_MD_superiorcorticostriate)[6,7]
data_subsample_WMlat_85_GGSEG[15,7]<-lavInspect(fit_sem_MD_superiorcorticostriate,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[15,] <-parameterEstimates(fit_sem_MD_superiorcorticostriate,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriatefrontalcortex
Model_MD_superiorcorticostriatefrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriatefrontalcortex_md'
fit_sem_MD_superiorcorticostriatefrontalcortex <- sem(Model_MD_superiorcorticostriatefrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_superiorcorticostriatefrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[16,5]<-BIC(fit_sem_MD_superiorcorticostriatefrontalcortex)
data_subsample_WMlat_85_GGSEG[16,6]<-parameterEstimates(fit_sem_MD_superiorcorticostriatefrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[16,7]<-lavInspect(fit_sem_MD_superiorcorticostriatefrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[16,] <-parameterEstimates(fit_sem_MD_superiorcorticostriatefrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriateparietalcortex
Model_MD_superiorcorticostriateparietalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriateparietalcortex_md'
fit_sem_MD_superiorcorticostriateparietalcortex <- sem(Model_MD_superiorcorticostriateparietalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_superiorcorticostriateparietalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[17,5]<-BIC(fit_sem_MD_superiorcorticostriateparietalcortex)
data_subsample_WMlat_85_GGSEG[17,6]<-parameterEstimates(fit_sem_MD_superiorcorticostriateparietalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[17,7]<-lavInspect(fit_sem_MD_superiorcorticostriateparietalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[17,] <-parameterEstimates(fit_sem_MD_superiorcorticostriateparietalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##striatalinferiorfrontalcortex
Model_MD_striatalinferiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ striatalinferiorfrontalcortex_md'
fit_sem_MD_striatalinferiorfrontalcortex <- sem(Model_MD_striatalinferiorfrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_striatalinferiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[18,5]<-BIC(fit_sem_MD_striatalinferiorfrontalcortex)
data_subsample_WMlat_85_GGSEG[18,6]<-parameterEstimates(fit_sem_MD_striatalinferiorfrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[18,7]<-lavInspect(fit_sem_MD_striatalinferiorfrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[18,] <-parameterEstimates(fit_sem_MD_striatalinferiorfrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorfrontalsuperiorfrontalcortex
Model_MD_inferiorfrontalsuperiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorfrontalsuperiorfrontalcortex_md'
fit_sem_MD_inferiorfrontalsuperiorfrontalcortex <- sem(Model_MD_inferiorfrontalsuperiorfrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[19,5]<-BIC(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex)
data_subsample_WMlat_85_GGSEG[19,6]<-parameterEstimates(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[19,7]<-lavInspect(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[19,] <-parameterEstimates(fit_sem_MD_inferiorfrontalsuperiorfrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##fornix_exfimbria
Model_MD_fornix_exfimbria <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fornix_exfimbria_md'
fit_sem_MD_fornix_exfimbria <- sem(Model_MD_fornix_exfimbria, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MD_fornix_exfimbria, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[20,5]<-BIC(fit_sem_MD_fornix_exfimbria)
data_subsample_WMlat_85_GGSEG[20,6]<-parameterEstimates(fit_sem_MD_fornix_exfimbria)[6,7]
data_subsample_WMlat_85_GGSEG[20,7]<-lavInspect(fit_sem_MD_fornix_exfimbria,what = "std.all")$beta[1,2]
paramest_WMlat_MD_85[20,] <-parameterEstimates(fit_sem_MD_fornix_exfimbria,standardized = TRUE, rsquare=TRUE)[6,]

write.csv(paramest_WMlat_MD_85, "paramest_WMlat_MD_85.csv", row.names=FALSE) 

#White Matter Volume (20 regions)
##fornix
Model_WMV_fornix <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fornix_wmv'
fit_sem_WMV_fornix <- sem(Model_WMV_fornix, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_fornix, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[1,8]<-BIC(fit_sem_WMV_fornix)
data_subsample_WMlat_85_GGSEG[1,9]<-parameterEstimates(fit_sem_WMV_fornix)[6,7]
data_subsample_WMlat_85_GGSEG[1,10]<-lavInspect(fit_sem_WMV_fornix,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[1,] <-parameterEstimates(fit_sem_WMV_fornix,standardized = TRUE, rsquare=TRUE)[6,]

##cingulatecingulum
Model_WMV_cingulatecingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ cingulatecingulum_wmv'
fit_sem_WMV_cingulatecingulum <- sem(Model_WMV_cingulatecingulum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_cingulatecingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[2,8]<-BIC(fit_sem_WMV_cingulatecingulum)
data_subsample_WMlat_85_GGSEG[2,9]<-parameterEstimates(fit_sem_WMV_cingulatecingulum)[6,7]
data_subsample_WMlat_85_GGSEG[2,10]<-lavInspect(fit_sem_WMV_cingulatecingulum,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[2,] <-parameterEstimates(fit_sem_WMV_cingulatecingulum,standardized = TRUE, rsquare=TRUE)[6,]

##parahippocampalcingulum
Model_WMV_parahippocampalcingulum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parahippocampalcingulum_wmv'
fit_sem_WMV_parahippocampalcingulum <- sem(Model_WMV_parahippocampalcingulum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_parahippocampalcingulum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[3,8]<-BIC(fit_sem_WMV_parahippocampalcingulum)
data_subsample_WMlat_85_GGSEG[3,9]<-parameterEstimates(fit_sem_WMV_parahippocampalcingulum)[6,7]
data_subsample_WMlat_85_GGSEG[3,10]<-lavInspect(fit_sem_WMV_parahippocampalcingulum,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[3,] <-parameterEstimates(fit_sem_WMV_parahippocampalcingulum,standardized = TRUE, rsquare=TRUE)[6,]

##corticospinalpyramidal
Model_WMV_corticospinalpyramidal <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ corticospinalpyramidal_wmv'
fit_sem_WMV_corticospinalpyramidal <- sem(Model_WMV_corticospinalpyramidal, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_corticospinalpyramidal, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[4,8]<-BIC(fit_sem_WMV_corticospinalpyramidal)
data_subsample_WMlat_85_GGSEG[4,9]<-parameterEstimates(fit_sem_WMV_corticospinalpyramidal)[6,7]
data_subsample_WMlat_85_GGSEG[4,10]<-lavInspect(fit_sem_WMV_corticospinalpyramidal,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[4,] <-parameterEstimates(fit_sem_WMV_corticospinalpyramidal,standardized = TRUE, rsquare=TRUE)[6,]

##anteriorthalamicradiations
Model_WMV_anteriorthalamicradiations <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ anteriorthalamicradiations_wmv'
fit_sem_WMV_anteriorthalamicradiations <- sem(Model_WMV_anteriorthalamicradiations, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_anteriorthalamicradiations, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[5,8]<-BIC(fit_sem_WMV_anteriorthalamicradiations)
data_subsample_WMlat_85_GGSEG[5,9]<-parameterEstimates(fit_sem_WMV_anteriorthalamicradiations)[6,7]
data_subsample_WMlat_85_GGSEG[5,10]<-lavInspect(fit_sem_WMV_anteriorthalamicradiations,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[5,] <-parameterEstimates(fit_sem_WMV_anteriorthalamicradiations,standardized = TRUE, rsquare=TRUE)[6,]

##uncinate
Model_WMV_uncinate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ uncinate_wmv'
fit_sem_WMV_uncinate <- sem(Model_WMV_uncinate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_uncinate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[6,8]<-BIC(fit_sem_WMV_uncinate)
data_subsample_WMlat_85_GGSEG[6,9]<-parameterEstimates(fit_sem_WMV_uncinate)[6,7]
data_subsample_WMlat_85_GGSEG[6,10]<-lavInspect(fit_sem_WMV_uncinate,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[6,] <-parameterEstimates(fit_sem_WMV_uncinate,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorlongitudinalfasiculus
Model_WMV_inferiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorlongitudinalfasiculus_wmv'
fit_sem_WMV_inferiorlongitudinalfasiculus <- sem(Model_WMV_inferiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_inferiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[7,8]<-BIC(fit_sem_WMV_inferiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[7,9]<-parameterEstimates(fit_sem_WMV_inferiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[7,10]<-lavInspect(fit_sem_WMV_inferiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[7,] <-parameterEstimates(fit_sem_WMV_inferiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorfrontooccipitalfasiculus
Model_WMV_inferiorfrontooccipitalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorfrontooccipitalfasiculus_wmv'
fit_sem_WMV_inferiorfrontooccipitalfasiculus <- sem(Model_WMV_inferiorfrontooccipitalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_inferiorfrontooccipitalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[8,8]<-BIC(fit_sem_WMV_inferiorfrontooccipitalfasiculus)
data_subsample_WMlat_85_GGSEG[8,9]<-parameterEstimates(fit_sem_WMV_inferiorfrontooccipitalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[8,10]<-lavInspect(fit_sem_WMV_inferiorfrontooccipitalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[8,] <-parameterEstimates(fit_sem_WMV_inferiorfrontooccipitalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##forcepsmajor
Model_WMV_forcepsmajor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ forcepsmajor_wmv'
fit_sem_WMV_forcepsmajor <- sem(Model_WMV_forcepsmajor, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_forcepsmajor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[9,8]<-BIC(fit_sem_WMV_forcepsmajor)
data_subsample_WMlat_85_GGSEG[9,9]<-parameterEstimates(fit_sem_WMV_forcepsmajor)[6,7]
data_subsample_WMlat_85_GGSEG[9,10]<-lavInspect(fit_sem_WMV_forcepsmajor,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[9,] <-parameterEstimates(fit_sem_WMV_forcepsmajor,standardized = TRUE, rsquare=TRUE)[6,]

##forcepsminor
Model_WMV_forcepsminor <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ forcepsminor_wmv'
fit_sem_WMV_forcepsminor <- sem(Model_WMV_forcepsminor, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_forcepsminor, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[10,8]<-BIC(fit_sem_WMV_forcepsminor)
data_subsample_WMlat_85_GGSEG[10,9]<-parameterEstimates(fit_sem_WMV_forcepsminor)[6,7]
data_subsample_WMlat_85_GGSEG[10,10]<-lavInspect(fit_sem_WMV_forcepsminor,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[10,] <-parameterEstimates(fit_sem_WMV_forcepsminor,standardized = TRUE, rsquare=TRUE)[6,]

##corpuscallosum
Model_WMV_corpuscallosum <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ corpuscallosum_wmv'
fit_sem_WMV_corpuscallosum <- sem(Model_WMV_corpuscallosum, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_corpuscallosum, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[11,8]<-BIC(fit_sem_WMV_corpuscallosum)
data_subsample_WMlat_85_GGSEG[11,9]<-parameterEstimates(fit_sem_WMV_corpuscallosum)[6,7]
data_subsample_WMlat_85_GGSEG[11,10]<-lavInspect(fit_sem_WMV_corpuscallosum,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[11,] <-parameterEstimates(fit_sem_WMV_corpuscallosum,standardized = TRUE, rsquare=TRUE)[6,]

##superiorlongitudinalfasiculus
Model_WMV_superiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorlongitudinalfasiculus_wmv'
fit_sem_WMV_superiorlongitudinalfasiculus <- sem(Model_WMV_superiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_superiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[12,8]<-BIC(fit_sem_WMV_superiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[12,9]<-parameterEstimates(fit_sem_WMV_superiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[12,10]<-lavInspect(fit_sem_WMV_superiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[12,] <-parameterEstimates(fit_sem_WMV_superiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##temporalsuperiorlongitudinalfasiculus
Model_WMV_temporalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ temporalsuperiorlongitudinalfasiculus_wmv'
fit_sem_WMV_temporalsuperiorlongitudinalfasiculus <- sem(Model_WMV_temporalsuperiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[13,8]<-BIC(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[13,9]<-parameterEstimates(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[13,10]<-lavInspect(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[13,] <-parameterEstimates(fit_sem_WMV_temporalsuperiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##parietalsuperiorlongitudinalfasiculus
Model_WMV_parietalsuperiorlongitudinalfasiculus <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ parietalsuperiorlongitudinalfasiculus_wmv'
fit_sem_WMV_parietalsuperiorlongitudinalfasiculus <- sem(Model_WMV_parietalsuperiorlongitudinalfasiculus, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[14,8]<-BIC(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus)
data_subsample_WMlat_85_GGSEG[14,9]<-parameterEstimates(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus)[6,7]
data_subsample_WMlat_85_GGSEG[14,10]<-lavInspect(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[14,] <-parameterEstimates(fit_sem_WMV_parietalsuperiorlongitudinalfasiculus,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriate
Model_WMV_superiorcorticostriate <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriate_wmv'
fit_sem_WMV_superiorcorticostriate <- sem(Model_WMV_superiorcorticostriate, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_superiorcorticostriate, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[15,8]<-BIC(fit_sem_WMV_superiorcorticostriate)
data_subsample_WMlat_85_GGSEG[15,9]<-parameterEstimates(fit_sem_WMV_superiorcorticostriate)[6,7]
data_subsample_WMlat_85_GGSEG[15,10]<-lavInspect(fit_sem_WMV_superiorcorticostriate,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[15,] <-parameterEstimates(fit_sem_WMV_superiorcorticostriate,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriatefrontalcortex
Model_WMV_superiorcorticostriatefrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriatefrontalcortex_wmv'
fit_sem_WMV_superiorcorticostriatefrontalcortex <- sem(Model_WMV_superiorcorticostriatefrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_superiorcorticostriatefrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[16,8]<-BIC(fit_sem_WMV_superiorcorticostriatefrontalcortex)
data_subsample_WMlat_85_GGSEG[16,9]<-parameterEstimates(fit_sem_WMV_superiorcorticostriatefrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[16,10]<-lavInspect(fit_sem_WMV_superiorcorticostriatefrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[16,] <-parameterEstimates(fit_sem_WMV_superiorcorticostriatefrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##superiorcorticostriateparietalcortex
Model_WMV_superiorcorticostriateparietalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ superiorcorticostriateparietalcortex_wmv'
fit_sem_WMV_superiorcorticostriateparietalcortex <- sem(Model_WMV_superiorcorticostriateparietalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_superiorcorticostriateparietalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[17,8]<-BIC(fit_sem_WMV_superiorcorticostriateparietalcortex)
data_subsample_WMlat_85_GGSEG[17,9]<-parameterEstimates(fit_sem_WMV_superiorcorticostriateparietalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[17,10]<-lavInspect(fit_sem_WMV_superiorcorticostriateparietalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[17,] <-parameterEstimates(fit_sem_WMV_superiorcorticostriateparietalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##striatalinferiorfrontalcortex
Model_WMV_striatalinferiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ striatalinferiorfrontalcortex_wmv'
fit_sem_WMV_striatalinferiorfrontalcortex <- sem(Model_WMV_striatalinferiorfrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_striatalinferiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[18,8]<-BIC(fit_sem_WMV_striatalinferiorfrontalcortex)
data_subsample_WMlat_85_GGSEG[18,9]<-parameterEstimates(fit_sem_WMV_striatalinferiorfrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[18,10]<-lavInspect(fit_sem_WMV_striatalinferiorfrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[18,] <-parameterEstimates(fit_sem_WMV_striatalinferiorfrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##inferiorfrontalsuperiorfrontalcortex
Model_WMV_inferiorfrontalsuperiorfrontalcortex <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ inferiorfrontalsuperiorfrontalcortex_wmv'
fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex <- sem(Model_WMV_inferiorfrontalsuperiorfrontalcortex, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[19,8]<-BIC(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex)
data_subsample_WMlat_85_GGSEG[19,9]<-parameterEstimates(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex)[6,7]
data_subsample_WMlat_85_GGSEG[19,10]<-lavInspect(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[19,] <-parameterEstimates(fit_sem_WMV_inferiorfrontalsuperiorfrontalcortex,standardized = TRUE, rsquare=TRUE)[6,]

##fornix_exfimbria
Model_WMV_fornix_exfimbria <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                          Cognitive_factor ~ fornix_exfimbria_wmv'
fit_sem_WMV_fornix_exfimbria <- sem(Model_WMV_fornix_exfimbria, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMV_fornix_exfimbria, fit.measures=TRUE,rsquare=T,standardized=T)
data_subsample_WMlat_85_GGSEG[20,8]<-BIC(fit_sem_WMV_fornix_exfimbria)
data_subsample_WMlat_85_GGSEG[20,9]<-parameterEstimates(fit_sem_WMV_fornix_exfimbria)[6,7]
data_subsample_WMlat_85_GGSEG[20,10]<-lavInspect(fit_sem_WMV_fornix_exfimbria,what = "std.all")$beta[1,2]
paramest_WMlat_WMV_85[20,] <-parameterEstimates(fit_sem_WMV_fornix_exfimbria,standardized = TRUE, rsquare=TRUE)[6,]

write.csv(paramest_WMlat_WMV_85, "paramest_WMlat_WMV_85.csv", row.names=FALSE) 

#Data Visualization for the individual models in the three white matter metrics
#Select all the std for the white matter metrics
data_subsample_WMlat_85_GGSEG_FA <- data_subsample_WMlat_85_GGSEG %>%
  dplyr::select(-c("FA_fit","FA_p","MD_fit","MD_p","MD_std","WMV_fit","WMV_p","WMV_std")) %>%
  rename(c("FA_std" = "Std")) %>%
  add_column(Metric = "FA")

data_subsample_WMlat_85_GGSEG_MD <- data_subsample_WMlat_85_GGSEG %>%
  dplyr::select(-c("FA_fit","FA_p","FA_std","MD_fit","MD_p","WMV_fit","WMV_p","WMV_std")) %>%
  rename(c("MD_std" = "Std")) %>%
  add_column(Metric = "MD")

data_subsample_WMlat_85_GGSEG_WMV <- data_subsample_WMlat_85_GGSEG %>%
  dplyr::select(-c("FA_fit","FA_p","FA_std","MD_fit","MD_p","MD_std","WMV_fit","WMV_p")) %>%
  rename(c("WMV_std" = "Std")) %>%
  add_column(Metric = "WMV")

#Compute the maximum and the minimum
min(data_subsample_WMlat_85_GGSEG_FA[,2])
max(data_subsample_WMlat_85_GGSEG_FA[,2])
min(data_subsample_WMlat_85_GGSEG_MD[,2])
max(data_subsample_WMlat_85_GGSEG_MD[,2])
min(data_subsample_WMlat_85_GGSEG_WMV[,2])
max(data_subsample_WMlat_85_GGSEG_WMV[,2])

data_subsample_STD_WMlat_85 <- rbind(data_subsample_WMlat_85_GGSEG_FA,data_subsample_WMlat_85_GGSEG_MD,by=c("label"))
data_subsample_STD_WMlat_85 <- rbind(data_subsample_STD_WMlat_85,data_subsample_WMlat_85_GGSEG_WMV,by=c("label"))
data_subsample_STD_WMlat_85<- data_subsample_STD_WMlat_85[-c(41,62),] 

data_subsample_STD_WMlat_85$mean <- ave(as.numeric(data_subsample_STD_WMlat_85$Std), data_subsample_STD_WMlat_85$Metric)

ggplot(data_subsample_STD_WMlat_85,aes(as.numeric(Std),fill=Metric)) + 
  geom_histogram(binwidth=0.005)+
  geom_density(adjust=2,alpha=0.2) +
  facet_grid(Metric~.) +
  geom_vline(aes(xintercept = mean, group = Metric), linetype="dashed", size=1,alpha=0.4) +
  theme_classic(base_size = 30)+
  xlab('Standardized estimated model parameters') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_fill_viridis(discrete=TRUE, option="mako",begin=0.3,end=0.9,name="Metrics",labels=c("Fractional Anisotropy","Mean Diffusivity","White Matter Volume"))
ggsave("plot_WMlat_85_std.png")


#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Model per metric - fractional anisotropy
#-------------------------------------------------------------------------------------------------------------------------------

#Regularized fractional anisotropy model
Model_FAlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ forcepsminor_fa + corpuscallosum_fa + cingulatecingulum_fa + corticospinalpyramidal_fa + anteriorthalamicradiations_fa + inferiorfrontooccipitalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa'

fit_sem_FAlat_85_reg_free <- sem(Model_FAlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FAlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized fractional anisotropy model
Model_FAlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*forcepsminor_fa + a*corpuscallosum_fa + a*cingulatecingulum_fa + a*corticospinalpyramidal_fa + a*anteriorthalamicradiations_fa + a*inferiorfrontooccipitalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa'

fit_sem_FAlat_85_reg_constrained <- sem(Model_FAlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FAlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_FAlat_85_reg_free,fit_sem_FAlat_85_reg_constrained))

#Table for the regularized fractional anisotropy model
BIC(fit_sem_FAlat_85_reg_free)
data_subsample_FAlat_85_rl_reg_GGSEG <- data.frame(label=c("forcepsminor","corpuscallosum","lh_cingulatecingulum","rh_cingulatecingulum","lh_corticospinalpyramidal","rh_corticospinalpyramidal","lh_anteriorthalamicradiations","rh_anteriorthalamicradiations","lh_inferiorfrontooccipitalfasiculus","rh_inferiorfrontooccipitalfasiculus","lh_temporalsuperiorlongitudinalfasiculus","rh_temporalsuperiorlongitudinalfasiculus"),
                                                   FA_p=rep(NA, times=12),
                                                   FA_std=rep(NA, times=12))
data_subsample_FAlat_85_rl_reg_GGSEG[1,2]<-parameterEstimates(fit_sem_FAlat_85_reg_free)[6,7]
data_subsample_FAlat_85_rl_reg_GGSEG[1,3]<-lavInspect(fit_sem_FAlat_85_reg_free,what = "std.all")$beta[1,2]
data_subsample_FAlat_85_rl_reg_GGSEG[2,2]<-parameterEstimates(fit_sem_FAlat_85_reg_free)[7,7]
data_subsample_FAlat_85_rl_reg_GGSEG[2,3]<-lavInspect(fit_sem_FAlat_85_reg_free,what = "std.all")$beta[1,3]

a=3
for (i in seq(from = 3, to = 12, by = 2)) {
  data_subsample_FAlat_85_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_FAlat_85_reg_free)[a+5,7]
  data_subsample_FAlat_85_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_FAlat_85_reg_free)[a+5,7]
  data_subsample_FAlat_85_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_FAlat_85_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_FAlat_85_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_FAlat_85_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}

#Model for fractional anisotropy with all the regions of interest (20 ROIs)
Model_FAlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ fornix_fa + cingulatecingulum_fa + parahippocampalcingulum_fa + corticospinalpyramidal_fa + anteriorthalamicradiations_fa + uncinate_fa + inferiorlongitudinalfasiculus_fa + inferiorfrontooccipitalfasiculus_fa + forcepsmajor_fa + forcepsminor_fa + corpuscallosum_fa + superiorlongitudinalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + parietalsuperiorlongitudinalfasiculus_fa + superiorcorticostriate_fa + superiorcorticostriatefrontalcortex_fa + superiorcorticostriateparietalcortex_fa + striatalinferiorfrontalcortex_fa + inferiorfrontalsuperiorfrontalcortex_fa + fornix_exfimbria_fa'

fit_sem_FAlat_85_all_free <- sem(Model_FAlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FAlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the fractional anisotropy model with all the regions
Model_FAlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ a*fornix_fa + a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*corticospinalpyramidal_fa + a*anteriorthalamicradiations_fa + a*uncinate_fa + a*inferiorlongitudinalfasiculus_fa + a*inferiorfrontooccipitalfasiculus_fa + a*forcepsmajor_fa + a*forcepsminor_fa + a*corpuscallosum_fa + a*superiorlongitudinalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*parietalsuperiorlongitudinalfasiculus_fa + a*superiorcorticostriate_fa + a*superiorcorticostriatefrontalcortex_fa + a*superiorcorticostriateparietalcortex_fa + a*striatalinferiorfrontalcortex_fa + a*inferiorfrontalsuperiorfrontalcortex_fa + a*fornix_exfimbria_fa'

fit_sem_FAlat_85_all_constrained <- sem(Model_FAlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_FAlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_FAlat_85_all_free,fit_sem_FAlat_85_all_constrained))

#Create tables to save the parameters of the regularized model and the model with all the ROIs for fractional anisotropy
paramest_WMlat_FAall_85 <- data.frame(lhs=rep(NA, times=20),
                                      op=rep(NA, times=20),
                                      rhs=rep(NA, times=20),
                                      est=rep(NA, times=20),
                                      se=rep(NA, times=20),
                                      z=rep(NA, times=20),
                                      pvalue=rep(NA, times=20),
                                      ci.lower=rep(NA, times=20),
                                      ci.upper=rep(NA, times=20),
                                      std.lv=rep(NA, times=20),
                                      std.all=rep(NA, times=20),
                                      std.nox=rep(NA, times=20))

paramest_WMlat_FAall_85[,] <-parameterEstimates(fit_sem_FAlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:25,]
write.csv(paramest_WMlat_FAall_85, "paramest_WMlat_FAall_85.csv", row.names=FALSE)

paramest_WMlat_FAreg_85 <- data.frame(lhs=rep(NA, times=7),
                                      op=rep(NA, times=7),
                                      rhs=rep(NA, times=7),
                                      est=rep(NA, times=7),
                                      se=rep(NA, times=7),
                                      z=rep(NA, times=7),
                                      pvalue=rep(NA, times=7),
                                      ci.lower=rep(NA, times=7),
                                      ci.upper=rep(NA, times=7),
                                      std.lv=rep(NA, times=7),
                                      std.all=rep(NA, times=7),
                                      std.nox=rep(NA, times=7))

paramest_WMlat_FAreg_85[,] <-parameterEstimates(fit_sem_FAlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:12,]
write.csv(paramest_WMlat_FAreg_85, "paramest_WMlat_FAreg_85.csv", row.names=FALSE)

#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Model per metric - mean diffusivity
#-------------------------------------------------------------------------------------------------------------------------------
  
#Regularized mean diffusivity model
Model_MDlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ forcepsminor_md + corpuscallosum_md + parahippocampalcingulum_md + anteriorthalamicradiations_md + superiorcorticostriate_md'

fit_sem_MDlat_85_reg_free <- sem(Model_MDlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MDlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized mean diffusivity model
Model_MDlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*forcepsminor_md + a*corpuscallosum_md + a*parahippocampalcingulum_md + a*anteriorthalamicradiations_md + a*superiorcorticostriate_md'

fit_sem_MDlat_85_reg_constrained <- sem(Model_MDlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MDlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_MDlat_85_reg_free,fit_sem_MDlat_85_reg_constrained))

#Table for the regularized mean diffusivity model
BIC(fit_sem_MDlat_85_reg_free)
data_subsample_MDlat_85_rl_reg_GGSEG <- data.frame(label=c("forcepsminor","corpuscallosum","lh_parahippocampalcingulum","rh_parahippocampalcingulum","lh_anteriorthalamicradiations","rh_anteriorthalamicradiations","lh_superiorcorticostriate","rh_superiorcorticostriate"),
                                                   MD_p=rep(NA, times=8),
                                                   MD_std=rep(NA, times=8))
data_subsample_MDlat_85_rl_reg_GGSEG[1,2]<-parameterEstimates(fit_sem_MDlat_85_reg_free)[6,7]
data_subsample_MDlat_85_rl_reg_GGSEG[1,3]<-lavInspect(fit_sem_MDlat_85_reg_free,what = "std.all")$beta[1,2]
data_subsample_MDlat_85_rl_reg_GGSEG[2,2]<-parameterEstimates(fit_sem_MDlat_85_reg_free)[7,7]
data_subsample_MDlat_85_rl_reg_GGSEG[2,3]<-lavInspect(fit_sem_MDlat_85_reg_free,what = "std.all")$beta[1,3]

a=3
for (i in seq(from = 3, to = 8, by = 2)) {
  data_subsample_MDlat_85_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_MDlat_85_reg_free)[a+5,7]
  data_subsample_MDlat_85_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_MDlat_85_reg_free)[a+5,7]
  data_subsample_MDlat_85_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_MDlat_85_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_MDlat_85_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_MDlat_85_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}


#Model for mean diffusivity with all the regions of interest (20 ROIs)
Model_MDlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ fornix_md + cingulatecingulum_md + parahippocampalcingulum_md + corticospinalpyramidal_md + anteriorthalamicradiations_md + uncinate_md + inferiorlongitudinalfasiculus_md + inferiorfrontooccipitalfasiculus_md + forcepsmajor_md + forcepsminor_md + corpuscallosum_md + superiorlongitudinalfasiculus_md + temporalsuperiorlongitudinalfasiculus_md + parietalsuperiorlongitudinalfasiculus_md + superiorcorticostriate_md + superiorcorticostriatefrontalcortex_md + superiorcorticostriateparietalcortex_md + striatalinferiorfrontalcortex_md + inferiorfrontalsuperiorfrontalcortex_md + fornix_exfimbria_md'

fit_sem_MDlat_85_all_free <- sem(Model_MDlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MDlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the mean diffusivity model with all the regions
Model_MDlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                Cognitive_factor ~ a*fornix_md + a*cingulatecingulum_md + a*parahippocampalcingulum_md + a*corticospinalpyramidal_md + a*anteriorthalamicradiations_md + a*uncinate_md + a*inferiorlongitudinalfasiculus_md + a*inferiorfrontooccipitalfasiculus_md + a*forcepsmajor_md + a*forcepsminor_md + a*corpuscallosum_md + a*superiorlongitudinalfasiculus_md + a*temporalsuperiorlongitudinalfasiculus_md + a*parietalsuperiorlongitudinalfasiculus_md + a*superiorcorticostriate_md + a*superiorcorticostriatefrontalcortex_md + a*superiorcorticostriateparietalcortex_md + a*striatalinferiorfrontalcortex_md + a*inferiorfrontalsuperiorfrontalcortex_md + a*fornix_exfimbria_md'

fit_sem_MDlat_85_all_constrained <- sem(Model_MDlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_MDlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_MDlat_85_all_free,fit_sem_MDlat_85_all_constrained))

#Create tables to save the parameters of the regularized model and the model with all the ROIs for mean diffusivity  
paramest_WMlat_MDall_85 <- data.frame(lhs=rep(NA, times=20),
                                      op=rep(NA, times=20),
                                      rhs=rep(NA, times=20),
                                      est=rep(NA, times=20),
                                      se=rep(NA, times=20),
                                      z=rep(NA, times=20),
                                      pvalue=rep(NA, times=20),
                                      ci.lower=rep(NA, times=20),
                                      ci.upper=rep(NA, times=20),
                                      std.lv=rep(NA, times=20),
                                      std.all=rep(NA, times=20),
                                      std.nox=rep(NA, times=20))

paramest_WMlat_MDall_85[,] <-parameterEstimates(fit_sem_MDlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:25,]
write.csv(paramest_WMlat_MDall_85, "paramest_WMlat_MDall_85.csv", row.names=FALSE)

paramest_WMlat_MDreg_85 <- data.frame(lhs=rep(NA, times=5),
                                      op=rep(NA, times=5),
                                      rhs=rep(NA, times=5),
                                      est=rep(NA, times=5),
                                      se=rep(NA, times=5),
                                      z=rep(NA, times=5),
                                      pvalue=rep(NA, times=5),
                                      ci.lower=rep(NA, times=5),
                                      ci.upper=rep(NA, times=5),
                                      std.lv=rep(NA, times=5),
                                      std.all=rep(NA, times=5),
                                      std.nox=rep(NA, times=5))

paramest_WMlat_MDreg_85[,] <-parameterEstimates(fit_sem_MDlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:10,]
write.csv(paramest_WMlat_MDreg_85, "paramest_WMlat_MDreg_85.csv", row.names=FALSE)

#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Model per metric - white matter volume
#-------------------------------------------------------------------------------------------------------------------------------
  
#Regularized white matter volume model
Model_WMVlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ forcepsminor_wmv + cingulatecingulum_wmv + parahippocampalcingulum_wmv + corticospinalpyramidal_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + superiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv'

fit_sem_WMVlat_85_reg_free <- sem(Model_WMVlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMVlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized white matter volume model
Model_WMVlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                      Cognitive_factor ~ a*forcepsminor_wmv + a*cingulatecingulum_wmv + a*parahippocampalcingulum_wmv + a*corticospinalpyramidal_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*superiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv'

fit_sem_WMVlat_85_reg_constrained <- sem(Model_WMVlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMVlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_WMVlat_85_reg_free,fit_sem_WMVlat_85_reg_constrained))

#Table for the regularized white matter volume model
BIC(fit_sem_WMVlat_85_reg_free)
data_subsample_WMVlat_85_rl_reg_GGSEG <- data.frame(label=c("forcepsminor","corpuscallosum","lh_cingulatecingulum","rh_cingulatecingulum","lh_parahippocampalcingulum","rh_parahippocampalcingulum","lh_corticospinalpyramidal","rh_corticospinalpyramidal","lh_anteriorthalamicradiations","rh_anteriorthalamicradiations","lh_uncinate","rh_uncinate","lh_inferiorlongitudinalfasiculus","rh_inferiorlongitudinalfasiculus","lh_inferiorfrontooccipitalfasiculus","rh_inferiorfrontooccipitalfasiculus","lh_superiorlongitudinalfasiculus","rh_superiorlongitudinalfasiculus","lh_parietalsuperiorlongitudinalfasiculus","rh_parietalsuperiorlongitudinalfasiculus","lh_superiorcorticostriateparietalcortex","rh_superiorcorticostriateparietalcortex","lh_striatalinferiorfrontalcortex","rh_striatalinferiorfrontalcortex","lh_inferiorfrontalsuperiorfrontalcortex","rh_inferiorfrontalsuperiorfrontalcortex","lh_fornix_exfimbria","rh_fornix_exfimbria"),
                                                    MD_p=rep(NA, times=28),
                                                    MD_std=rep(NA, times=28))
data_subsample_WMVlat_85_rl_reg_GGSEG[1,2]<-parameterEstimates(fit_sem_WMVlat_85_reg_free)[6,7]
data_subsample_WMVlat_85_rl_reg_GGSEG[1,3]<-lavInspect(fit_sem_WMVlat_85_reg_free,what = "std.all")$beta[1,2]
data_subsample_WMVlat_85_rl_reg_GGSEG[2,2]<-parameterEstimates(fit_sem_WMVlat_85_reg_free)[7,7]
data_subsample_WMVlat_85_rl_reg_GGSEG[2,3]<-lavInspect(fit_sem_WMVlat_85_reg_free,what = "std.all")$beta[1,3]

a=3
for (i in seq(from = 3, to = 28, by = 2)) {
  data_subsample_WMVlat_85_rl_reg_GGSEG[i,2]<-parameterEstimates(fit_sem_WMVlat_85_reg_free)[a+5,7]
  data_subsample_WMVlat_85_rl_reg_GGSEG[i+1,2]<-parameterEstimates(fit_sem_WMVlat_85_reg_free)[a+5,7]
  data_subsample_WMVlat_85_rl_reg_GGSEG[i,3]<-lavInspect(fit_sem_WMVlat_85_reg_free,what = "std.all")$beta[1,a+1]
  data_subsample_WMVlat_85_rl_reg_GGSEG[i+1,3]<-lavInspect(fit_sem_WMVlat_85_reg_free,what = "std.all")$beta[1,a+1]
  a=a+1
}


#Model for white matter volume with all the regions of interest (20 ROIs)
Model_WMVlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                             Cognitive_factor ~ fornix_wmv + cingulatecingulum_wmv + parahippocampalcingulum_wmv + corticospinalpyramidal_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + forcepsmajor_wmv + forcepsminor_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + temporalsuperiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + superiorcorticostriate_wmv + superiorcorticostriatefrontalcortex_wmv + superiorcorticostriateparietalcortex_wmv + striatalinferiorfrontalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv + fornix_exfimbria_wmv'

fit_sem_WMVlat_85_all_free <- sem(Model_WMVlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMVlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the white matter volume model with all the regions
Model_WMVlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                    Cognitive_factor ~ a*fornix_wmv + a*cingulatecingulum_wmv + a*parahippocampalcingulum_wmv + a*corticospinalpyramidal_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*forcepsmajor_wmv + a*forcepsminor_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*temporalsuperiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*superiorcorticostriate_wmv + a*superiorcorticostriatefrontalcortex_wmv + a*superiorcorticostriateparietalcortex_wmv + a*striatalinferiorfrontalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv + a*fornix_exfimbria_wmv'

fit_sem_WMVlat_85_all_constrained <- sem(Model_WMVlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMVlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_WMVlat_85_all_free,fit_sem_WMVlat_85_all_constrained))

#Create tables to save the parameters of the regularized model and the model with all the ROIs for white matter volume  
paramest_WMlat_WMVall_85 <- data.frame(lhs=rep(NA, times=20),
                                       op=rep(NA, times=20),
                                       rhs=rep(NA, times=20),
                                       est=rep(NA, times=20),
                                       se=rep(NA, times=20),
                                       z=rep(NA, times=20),
                                       pvalue=rep(NA, times=20),
                                       ci.lower=rep(NA, times=20),
                                       ci.upper=rep(NA, times=20),
                                       std.lv=rep(NA, times=20),
                                       std.all=rep(NA, times=20),
                                       std.nox=rep(NA, times=20))

paramest_WMlat_WMVall_85[,] <-parameterEstimates(fit_sem_WMVlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:25,]
write.csv(paramest_WMlat_WMVall_85, "paramest_WMlat_WMVall_85.csv", row.names=FALSE)

paramest_WMlat_WMVreg_85 <- data.frame(lhs=rep(NA, times=10),
                                       op=rep(NA, times=10),
                                       rhs=rep(NA, times=10),
                                       est=rep(NA, times=10),
                                       se=rep(NA, times=10),
                                       z=rep(NA, times=10),
                                       pvalue=rep(NA, times=10),
                                       ci.lower=rep(NA, times=10),
                                       ci.upper=rep(NA, times=10),
                                       std.lv=rep(NA, times=10),
                                       std.all=rep(NA, times=10),
                                       std.nox=rep(NA, times=10))

paramest_WMlat_WMVreg_85[,] <-parameterEstimates(fit_sem_WMVlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:15,]
write.csv(paramest_WMlat_WMVreg_85, "paramest_WMlat_WMVreg_85.csv", row.names=FALSE)


#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Models including the three white matter metrics
#-------------------------------------------------------------------------------------------------------------------------------

#Regularized white matter metrics model 
Model_WMlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ cingulatecingulum_fa + parahippocampalcingulum_fa + uncinate_fa + inferiorfrontooccipitalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + inferiorfrontalsuperiorfrontalcortex_fa + fornix_exfimbria_fa + uncinate_md + inferiorfrontooccipitalfasiculus_md + forcepsmajor_md + forcepsminor_md + striatalinferiorfrontalcortex_md + fornix_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + forcepsmajor_wmv + forcepsminor_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_WMlat_85_reg_free <- sem(Model_WMlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized white matter metrics model
Model_WMlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*uncinate_fa + a*inferiorfrontooccipitalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*inferiorfrontalsuperiorfrontalcortex_fa + a*fornix_exfimbria_fa + a*uncinate_md + a*inferiorfrontooccipitalfasiculus_md + a*forcepsmajor_md + a*forcepsminor_md + a*striatalinferiorfrontalcortex_md + a*fornix_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*forcepsmajor_wmv + a*forcepsminor_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_WMlat_85_reg_constrained <- sem(Model_WMlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_WMlat_85_reg_free,fit_sem_WMlat_85_reg_constrained))

#Table for the regularized white matter metrics model
BIC(fit_sem_WMlat_85_reg_free)
data_subsample_WMlat_85_regFA_GGSEG <- data.frame(label=c("cingulatecingulum","parahippocampalcingulum","uncinate","inferiorfrontooccipitalfasiculus","temporalsuperiorlongitudinalfasiculus","inferiorfrontalsuperiorfrontalcortex","fornix_exfimbria"),
                                                  FA_p=rep(NA, times=7),
                                                  FA_std=rep(NA, times=7))
a=0
for (i in seq(from = 1, to =7)) {
  a=a+1
  data_subsample_WMlat_85_regFA_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_85_reg_free)[a+5,7]
  data_subsample_WMlat_85_regFA_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_85_reg_free,what = "std.all")$beta[1,a+1]
}

data_subsample_WMlat_85_regMD_GGSEG <- data.frame(label=c("uncinate","inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","striatalinferiorfrontalcortex"),
                                                  MD_p=rep(NA, times=5),
                                                  MD_std=rep(NA, times=5))
for (i in seq(from = 1, to = 5)) {
  a=a+1
  data_subsample_WMlat_85_regMD_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_85_reg_free)[a+5,7]
  data_subsample_WMlat_85_regMD_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_85_reg_free,what = "std.all")$beta[1,a+1]
  
}

data_subsample_WMlat_85_regWMV_GGSEG <- data.frame(label=c("fornix","anteriorthalamicradiations","uncinate","inferiorlongitudinalfasiculus","inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","corpuscallosum","superiorlongitudinalfasiculus","parietalsuperiorlongitudinalfasiculus","inferiorfrontalsuperiorfrontalcortex"),
                                                   WMV_p=rep(NA, times=11),
                                                   WMV_std=rep(NA, times=11))

for (i in seq(from = 1, to = 11)) {
  a=a+1
  data_subsample_WMlat_85_regWMV_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_85_reg_free)[a+5,7]
  data_subsample_WMlat_85_regWMV_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_85_reg_free,what = "std.all")$beta[1,a+1]
}

#Compare the information provide by the regularized regions in the three white matter metrics 
data_subsample_WMlat_85_reg_GGSEG_pvalue<-merge(data_subsample_WMlat_85_regFA_GGSEG[,c(1,2)],data_subsample_WMlat_85_regMD_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_WMlat_85_reg_GGSEG_pvalue<-merge(data_subsample_WMlat_85_reg_GGSEG_pvalue,data_subsample_WMlat_85_regWMV_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_WMlat_85_reg_GGSEG_pvalue, "STD_WMlat_85_reg_pvalue.csv", row.names=FALSE)

data_subsample_WMlat_85_reg_GGSEG_STD<-merge(data_subsample_WMlat_85_regFA_GGSEG[,c(1,3)],data_subsample_WMlat_85_regMD_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_WMlat_85_reg_GGSEG_STD<-merge(data_subsample_WMlat_85_reg_GGSEG_STD,data_subsample_WMlat_85_regWMV_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_WMlat_85_reg_GGSEG_STD, "STD_WMlat_85_reg_STD.csv", row.names=FALSE)

#Model for white matter metrics with all the regions of interest (60 ROIs)
Model_WMlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                            Cognitive_factor ~ fornix_fa + cingulatecingulum_fa + parahippocampalcingulum_fa + corticospinalpyramidal_fa + anteriorthalamicradiations_fa + uncinate_fa + inferiorlongitudinalfasiculus_fa + inferiorfrontooccipitalfasiculus_fa + forcepsmajor_fa + forcepsminor_fa + corpuscallosum_fa + superiorlongitudinalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + parietalsuperiorlongitudinalfasiculus_fa + superiorcorticostriate_fa + superiorcorticostriatefrontalcortex_fa + superiorcorticostriateparietalcortex_fa + striatalinferiorfrontalcortex_fa + inferiorfrontalsuperiorfrontalcortex_fa + fornix_exfimbria_fa + fornix_md + cingulatecingulum_md + parahippocampalcingulum_md + corticospinalpyramidal_md + anteriorthalamicradiations_md + uncinate_md + inferiorlongitudinalfasiculus_md + inferiorfrontooccipitalfasiculus_md + forcepsmajor_md + forcepsminor_md + corpuscallosum_md + superiorlongitudinalfasiculus_md + temporalsuperiorlongitudinalfasiculus_md + parietalsuperiorlongitudinalfasiculus_md + superiorcorticostriate_md + superiorcorticostriatefrontalcortex_md + superiorcorticostriateparietalcortex_md + striatalinferiorfrontalcortex_md + inferiorfrontalsuperiorfrontalcortex_md + fornix_exfimbria_md + fornix_wmv + cingulatecingulum_wmv + parahippocampalcingulum_wmv + corticospinalpyramidal_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + forcepsmajor_wmv + forcepsminor_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + temporalsuperiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + superiorcorticostriate_wmv + superiorcorticostriatefrontalcortex_wmv + superiorcorticostriateparietalcortex_wmv + striatalinferiorfrontalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv + fornix_exfimbria_wmv'

fit_sem_WMlat_85_all_free <- sem(Model_WMlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the white matter metrics model with all the regions
Model_WMlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ a*fornix_fa + a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*corticospinalpyramidal_fa + a*anteriorthalamicradiations_fa + a*uncinate_fa + a*inferiorlongitudinalfasiculus_fa + a*inferiorfrontooccipitalfasiculus_fa + a*forcepsmajor_fa + a*forcepsminor_fa + a*corpuscallosum_fa + a*superiorlongitudinalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*parietalsuperiorlongitudinalfasiculus_fa + a*superiorcorticostriate_fa + a*superiorcorticostriatefrontalcortex_fa + a*superiorcorticostriateparietalcortex_fa + a*striatalinferiorfrontalcortex_fa + a*inferiorfrontalsuperiorfrontalcortex_fa + a*fornix_exfimbria_fa + a*fornix_md + a*cingulatecingulum_md + a*parahippocampalcingulum_md + a*corticospinalpyramidal_md + a*anteriorthalamicradiations_md + a*uncinate_md + a*inferiorlongitudinalfasiculus_md + a*inferiorfrontooccipitalfasiculus_md + a*forcepsmajor_md + a*forcepsminor_md + a*corpuscallosum_md + a*superiorlongitudinalfasiculus_md + a*temporalsuperiorlongitudinalfasiculus_md + a*parietalsuperiorlongitudinalfasiculus_md + a*superiorcorticostriate_md + a*superiorcorticostriatefrontalcortex_md + a*superiorcorticostriateparietalcortex_md + a*striatalinferiorfrontalcortex_md + a*inferiorfrontalsuperiorfrontalcortex_md + a*fornix_exfimbria_md + a*fornix_wmv + a*cingulatecingulum_wmv + a*parahippocampalcingulum_wmv + a*corticospinalpyramidal_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*forcepsmajor_wmv + a*forcepsminor_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*temporalsuperiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*superiorcorticostriate_wmv + a*superiorcorticostriatefrontalcortex_wmv + a*superiorcorticostriateparietalcortex_wmv + a*striatalinferiorfrontalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv + a*fornix_exfimbria_wmv'

fit_sem_WMlat_85_all_constrained <- sem(Model_WMlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_WMlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_WMlat_85_all_free,fit_sem_WMlat_85_all_constrained))

#Table for the white matter metrics model
BIC(fit_sem_WMlat_85_all_free)
data_subsample_WMlat_85_allFA_GGSEG <- data.frame(label=c("fornix","cingulatecingulum","parahippocampalcingulum","corticospinalpyramidal","anteriorthalamicradiations","uncinate","inferiorlongitudinalfasiculus","inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","corpuscallosum","superiorlongitudinalfasiculus","temporalsuperiorlongitudinalfasiculus","parietalsuperiorlongitudinalfasiculus","superiorcorticostriate","superiorcorticostriatefrontalcortex","superiorcorticostriateparietalcortex","striatalinferiorfrontalcortex","inferiorfrontalsuperiorfrontalcortex","fornix_exfimbria"),
                                                  FA_p=rep(NA, times=20),
                                                  FA_std=rep(NA, times=20))
a=0
for (i in seq(from = 1, to = 20)) {
  a=a+1
  data_subsample_WMlat_85_allFA_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_85_all_free)[a+5,7]
  data_subsample_WMlat_85_allFA_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_85_all_free,what = "std.all")$beta[1,a+1]
}

data_subsample_WMlat_85_allMD_GGSEG <- data.frame(label=c("fornix","cingulatecingulum","parahippocampalcingulum","corticospinalpyramidal","anteriorthalamicradiations","uncinate","inferiorlongitudinalfasiculus","inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","corpuscallosum","superiorlongitudinalfasiculus","temporalsuperiorlongitudinalfasiculus","parietalsuperiorlongitudinalfasiculus","superiorcorticostriate","superiorcorticostriatefrontalcortex","superiorcorticostriateparietalcortex","striatalinferiorfrontalcortex","inferiorfrontalsuperiorfrontalcortex","fornix_exfimbria"),
                                                  MD_p=rep(NA, times=20),
                                                  MD_std=rep(NA, times=20))

for (i in seq(from = 1, to = 20)) {
  a=a+1
  data_subsample_WMlat_85_allMD_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_85_all_free)[a+5,7]
  data_subsample_WMlat_85_allMD_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_85_all_free,what = "std.all")$beta[1,a+1]
}

data_subsample_WMlat_85_allWMV_GGSEG <- data.frame(label=c("fornix","cingulatecingulum","parahippocampalcingulum","corticospinalpyramidal","anteriorthalamicradiations","uncinate","inferiorlongitudinalfasiculus","inferiorfrontooccipitalfasiculus","forcepsmajor","forcepsminor","corpuscallosum","superiorlongitudinalfasiculus","temporalsuperiorlongitudinalfasiculus","parietalsuperiorlongitudinalfasiculus","superiorcorticostriate","superiorcorticostriatefrontalcortex","superiorcorticostriateparietalcortex","striatalinferiorfrontalcortex","inferiorfrontalsuperiorfrontalcortex","fornix_exfimbria"),
                                                   WMV_p=rep(NA, times=20),
                                                   WMV_std=rep(NA, times=20))

for (i in seq(from = 1, to = 20)) {
  a=a+1
  data_subsample_WMlat_85_allWMV_GGSEG[i,2]<-parameterEstimates(fit_sem_WMlat_85_all_free)[a+5,7]
  data_subsample_WMlat_85_allWMV_GGSEG[i,3]<-lavInspect(fit_sem_WMlat_85_all_free,what = "std.all")$beta[1,a+1]
}

#Compare the information provide by the regions in the three white matter metrics 
data_subsample_WMlat_85_all_GGSEG_pvalue<-merge(data_subsample_WMlat_85_allFA_GGSEG[,c(1,2)],data_subsample_WMlat_85_allMD_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_WMlat_85_all_GGSEG_pvalue<-merge(data_subsample_WMlat_85_all_GGSEG_pvalue,data_subsample_WMlat_85_allWMV_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_WMlat_85_all_GGSEG_STD<-merge(data_subsample_WMlat_85_allFA_GGSEG[,c(1,3)],data_subsample_WMlat_85_allMD_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_WMlat_85_all_GGSEG_STD<-merge(data_subsample_WMlat_85_all_GGSEG_STD,data_subsample_WMlat_85_allWMV_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)

#Create tables to save the parameters of the regularized model and the model with the regularized ROIs for white matter metrics
paramest_WMlat_WMall_85 <- data.frame(lhs=rep(NA, times=60),
                                      op=rep(NA, times=60),
                                      rhs=rep(NA, times=60),
                                      est=rep(NA, times=60),
                                      se=rep(NA, times=60),
                                      z=rep(NA, times=60),
                                      pvalue=rep(NA, times=60),
                                      ci.lower=rep(NA, times=60),
                                      ci.upper=rep(NA, times=60),
                                      std.lv=rep(NA, times=60),
                                      std.all=rep(NA, times=60),
                                      std.nox=rep(NA, times=60))

paramest_WMlat_WMall_85[,] <-parameterEstimates(fit_sem_WMlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:65,]
write.csv(paramest_WMlat_WMall_85, "paramest_WMlat_WMall_85.csv", row.names=FALSE)

paramest_WMlat_WMreg_85 <- data.frame(lhs=rep(NA, times=23),
                                      op=rep(NA, times=23),
                                      rhs=rep(NA, times=23),
                                      est=rep(NA, times=23),
                                      se=rep(NA, times=23),
                                      z=rep(NA, times=23),
                                      pvalue=rep(NA, times=23),
                                      ci.lower=rep(NA, times=23),
                                      ci.upper=rep(NA, times=23),
                                      std.lv=rep(NA, times=23),
                                      std.all=rep(NA, times=23),
                                      std.nox=rep(NA, times=23))

paramest_WMlat_WMreg_85[,] <-parameterEstimates(fit_sem_WMlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:28,]
write.csv(paramest_WMlat_WMreg_85, "paramest_WMlat_WMreg_T085.csv", row.names=FALSE)


#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Models including grey & white matter metrics
#-------------------------------------------------------------------------------------------------------------------------------

#Regularized grey and white matter metrics model 
Model_GMWMlat_85_reg_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ caudalmiddlefrontal_ct + fusiform_ct + lateralorbitofrontal_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsorbitalis_ct + pericalcarine_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiortemporal_ct + insula_ct + bankssts_sa + cuneus_sa + entorhinal_sa + lateraloccipital_sa + parsorbitalis_sa + postcentral_sa + precentral_sa + superiortemporal_sa + fusiform_gmv + inferiorparietal_gmv + parahippocampal_gmv + paracentral_gmv + pericalcarine_gmv + precentral_gmv + precuneus_gmv + rostralmiddlefrontal_gmv + supramarginal_gmv + cingulatecingulum_fa + parahippocampalcingulum_fa + uncinate_fa + inferiorfrontooccipitalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + superiorcorticostriateparietalcortex_fa + fornix_exfimbria_fa + parahippocampalcingulum_md + inferiorlongitudinalfasiculus_md + forcepsmajor_md + temporalsuperiorlongitudinalfasiculus_md + striatalinferiorfrontalcortex_md + anteriorthalamicradiations_wmv + inferiorlongitudinalfasiculus_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + superiorcorticostriateparietalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_GMWMlat_85_reg_free <- sem(Model_GMWMlat_85_reg_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMWMlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the regularized grey & white matter metrics model
Model_GMWMlat_85_reg_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*caudalmiddlefrontal_ct + a*fusiform_ct + a*lateralorbitofrontal_ct + a*medialorbitofrontal_ct + a*middletemporal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsorbitalis_ct + a*pericalcarine_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiortemporal_ct + a*insula_ct + a*bankssts_sa + a*cuneus_sa + a*entorhinal_sa + a*lateraloccipital_sa + a*parsorbitalis_sa + a*postcentral_sa + a*precentral_sa + a*superiortemporal_sa + a*fusiform_gmv + a*inferiorparietal_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*pericalcarine_gmv + a*precentral_gmv + a*precuneus_gmv + a*rostralmiddlefrontal_gmv + a*supramarginal_gmv + a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*uncinate_fa + a*inferiorfrontooccipitalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*superiorcorticostriateparietalcortex_fa + a*fornix_exfimbria_fa + a*parahippocampalcingulum_md + a*inferiorlongitudinalfasiculus_md + a*forcepsmajor_md + a*temporalsuperiorlongitudinalfasiculus_md + a*striatalinferiorfrontalcortex_md + a*anteriorthalamicradiations_wmv + a*inferiorlongitudinalfasiculus_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*superiorcorticostriateparietalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_GMWMlat_85_reg_constrained <- sem(Model_GMWMlat_85_reg_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMWMlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

kable(anova(fit_sem_GMWMlat_85_reg_free,fit_sem_GMWMlat_85_reg_constrained))

#Table for the regularized white matter metrics model
BIC(fit_sem_GMWMlat_85_reg_free)
data_subsample_GMWMlat_85_regCT_GGSEG <- data.frame(label=c("caudalmiddlefrontal","fusiform","lateralorbitofrontal","medialorbitofrontal","middletemporal","parahippocampal","paracentral","parsorbitalis","pericalcarine","precuneus","rostralanteriorcingulate","rostralmiddlefrontal","superiorfrontal","superiortemporal","insula"),
                                                  CT_p=rep(NA, times=15),
                                                  CT_std=rep(NA, times=15))
a=0
for (i in seq(from = 1, to = 15)) {
  a=a+1
  data_subsample_GMWMlat_85_regCT_GGSEG[i,2]<-parameterEstimates(fit_sem_GMWMlat_85_reg_free)[a+5,7]
  data_subsample_GMWMlat_85_regCT_GGSEG[i,3]<-lavInspect(fit_sem_GMWMlat_85_reg_free,what = "std.all")$beta[1,a+1]
}

data_subsample_GMWMlat_85_regSA_GGSEG <- data.frame(label=c("bankssts","cuneus","entorhinal","lateraloccipital","parsorbitalis","postcentral","precentral","superiortemporal"),
                                                  SA_p=rep(NA, times=8),
                                                  SA_std=rep(NA, times=8))
for (i in seq(from = 1, to = 8)) {
  a=a+1
  data_subsample_GMWMlat_85_regSA_GGSEG[i,2]<-parameterEstimates(fit_sem_GMWMlat_85_reg_free)[a+5,7]
  data_subsample_GMWMlat_85_regSA_GGSEG[i,3]<-lavInspect(fit_sem_GMWMlat_85_reg_free,what = "std.all")$beta[1,a+1]
}

data_subsample_GMWMlat_85_regGMV_GGSEG <- data.frame(label=c("fusiform","inferiorparietal","parahippocampal","paracentral","pericalcarine","precentral","precuneus","rostralmiddlefrontal","supramarginal"),
                                                   GMV_p=rep(NA, times=9),
                                                   GMV_std=rep(NA, times=9))
for (i in seq(from = 1, to = 9)) {
  a=a+1  
  data_subsample_GMWMlat_85_regGMV_GGSEG[i,2]<-parameterEstimates(fit_sem_GMWMlat_85_reg_free)[a+5,7]
  data_subsample_GMWMlat_85_regGMV_GGSEG[i,3]<-lavInspect(fit_sem_GMWMlat_85_reg_free,what = "std.all")$beta[1,a+1]
}

data_subsample_GMWMlat_85_regFA_GGSEG <- data.frame(label=c("cingulatecingulum","parahippocampalcingulum","uncinate","inferiorfrontooccipitalfasiculus","temporalsuperiorlongitudinalfasiculus","superiorcorticostriateparietalcortex","fornix_exfimbria" ),
                                                  FA_p=rep(NA, times=7),
                                                  FA_std=rep(NA, times=7))
for (i in seq(from = 1, to =7)) {
  a=a+1
  data_subsample_GMWMlat_85_regFA_GGSEG[i,2]<-parameterEstimates(fit_sem_GMWMlat_85_reg_free)[a+5,7]
  data_subsample_GMWMlat_85_regFA_GGSEG[i,3]<-lavInspect(fit_sem_GMWMlat_85_reg_free,what = "std.all")$beta[1,a+1]
}

data_subsample_GMWMlat_85_regMD_GGSEG <- data.frame(label=c("parahippocampalcingulum","inferiorlongitudinalfasiculus","forcepsmajor","temporalsuperiorlongitudinalfasiculus","striatalinferiorfrontalcortex"),
                                                  MD_p=rep(NA, times=5),
                                                  MD_std=rep(NA, times=5))
for (i in seq(from = 1, to = 5)) {
  a=a+1
  data_subsample_GMWMlat_85_regMD_GGSEG[i,2]<-parameterEstimates(fit_sem_GMWMlat_85_reg_free)[a+5,7]
  data_subsample_GMWMlat_85_regMD_GGSEG[i,3]<-lavInspect(fit_sem_GMWMlat_85_reg_free,what = "std.all")$beta[1,a+1]
  
}

data_subsample_GMWMlat_85_regWMV_GGSEG <- data.frame(label=c("anteriorthalamicradiations","inferiorlongitudinalfasiculus","corpuscallosum","superiorlongitudinalfasiculus","superiorcorticostriateparietalcortex","inferiorfrontalsuperiorfrontalcortex"),
                                                   WMV_p=rep(NA, times=6),
                                                   WMV_std=rep(NA, times=6))

for (i in seq(from = 1, to = 6)) {
  a=a+1
  data_subsample_GMWMlat_85_regWMV_GGSEG[i,2]<-parameterEstimates(fit_sem_GMWMlat_85_reg_free)[a+5,7]
  data_subsample_GMWMlat_85_regWMV_GGSEG[i,3]<-lavInspect(fit_sem_GMWMlat_85_reg_free,what = "std.all")$beta[1,a+1]
}

#Compare the information provide by the regularized regions in the three white matter metrics 
data_subsample_GMWMlat_85_reg_GGSEG_GMpvalue<-merge(data_subsample_GMWMlat_85_regCT_GGSEG[,c(1,2)],data_subsample_GMWMlat_85_regSA_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_GMWMlat_85_reg_GGSEG_GMpvalue<-merge(data_subsample_GMWMlat_85_reg_GGSEG_GMpvalue,data_subsample_GMWMlat_85_regGMV_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_GMWMlat_85_reg_GGSEG_GMpvalue, "STD_GMWMlat_85_reg_GMpvalue.csv", row.names=FALSE)

data_subsample_GMWMlat_85_reg_GGSEG_WMpvalue<-merge(data_subsample_GMWMlat_85_regFA_GGSEG[,c(1,2)],data_subsample_GMWMlat_85_regMD_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_GMWMlat_85_reg_GGSEG_WMpvalue<-merge(data_subsample_GMWMlat_85_reg_GGSEG_WMpvalue,data_subsample_GMWMlat_85_regWMV_GGSEG[,c(1,2)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_GMWMlat_85_reg_GGSEG_WMpvalue, "STD_GMWMlat_85_reg_WMpvalue.csv", row.names=FALSE)

data_subsample_GMWMlat_85_reg_GGSEG_GMSTD<-merge(data_subsample_GMWMlat_85_regCT_GGSEG[,c(1,3)],data_subsample_GMWMlat_85_regSA_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_GMWMlat_85_reg_GGSEG_GMSTD<-merge(data_subsample_GMWMlat_85_reg_GGSEG_GMSTD,data_subsample_GMWMlat_85_regGMV_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_GMWMlat_85_reg_GGSEG_GMSTD, "STD_GMWMlat_85_reg_GMSTD.csv", row.names=FALSE)

data_subsample_GMWMlat_85_reg_GGSEG_WMSTD<-merge(data_subsample_GMWMlat_85_regFA_GGSEG[,c(1,3)],data_subsample_GMWMlat_85_regMD_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
data_subsample_GMWMlat_85_reg_GGSEG_WMSTD<-merge(data_subsample_GMWMlat_85_reg_GGSEG_WMSTD,data_subsample_GMWMlat_85_regWMV_GGSEG[,c(1,3)],by=c("label"),all.x=TRUE,all.y=TRUE)
write.csv(data_subsample_GMWMlat_85_reg_GGSEG_WMSTD, "STD_GMWMlat_85_reg_WMSTD.csv", row.names=FALSE)

#Model for grey & white matter metrics with all the regions of interest (162 ROIs) 
Model_GMWMlat_85_all_free <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                              Cognitive_factor ~ bankssts_ct + caudalanteriorcingulate_ct + caudalmiddlefrontal_ct + cuneus_ct + entorhinal_ct + fusiform_ct + inferiorparietal_ct + inferiortemporal_ct + isthmuscingulate_ct + lateraloccipital_ct + lateralorbitofrontal_ct + lingual_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsopercularis_ct + parsorbitalis_ct + parstriangularis_ct + pericalcarine_ct + postcentral_ct + posteriorcingulate_ct + precentral_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiorparietal_ct + superiortemporal_ct + supramarginal_ct + frontalpole_ct + temporalpole_ct + transversetemporal_ct + insula_ct + bankssts_sa + caudalanteriorcingulate_sa + caudalmiddlefrontal_sa + cuneus_sa + entorhinal_sa + fusiform_sa + inferiorparietal_sa + inferiortemporal_sa + isthmuscingulate_sa + lateraloccipital_sa + lateralorbitofrontal_sa + lingual_sa + medialorbitofrontal_sa + middletemporal_sa + parahippocampal_sa + paracentral_sa + parsopercularis_sa + parsorbitalis_sa + parstriangularis_sa + pericalcarine_sa + postcentral_sa + posteriorcingulate_sa + precentral_sa + precuneus_sa + rostralanteriorcingulate_sa + rostralmiddlefrontal_sa + superiorfrontal_sa + superiorparietal_sa + superiortemporal_sa + supramarginal_sa + frontalpole_sa + temporalpole_sa + transversetemporal_sa + insula_sa + bankssts_gmv + caudalanteriorcingulate_gmv + caudalmiddlefrontal_gmv + cuneus_gmv + entorhinal_gmv + fusiform_gmv + inferiorparietal_gmv + inferiortemporal_gmv + isthmuscingulate_gmv + lateraloccipital_gmv + lateralorbitofrontal_gmv + lingual_gmv + medialorbitofrontal_gmv + middletemporal_gmv + parahippocampal_gmv + paracentral_gmv + parsopercularis_gmv + parsorbitalis_gmv + parstriangularis_gmv + pericalcarine_gmv + postcentral_gmv + posteriorcingulate_gmv + precentral_gmv + precuneus_gmv + rostralanteriorcingulate_gmv + rostralmiddlefrontal_gmv + superiorfrontal_gmv + superiorparietal_gmv + superiortemporal_gmv + supramarginal_gmv + frontalpole_gmv + temporalpole_gmv + transversetemporal_gmv + insula_gmv + fornix_fa + cingulatecingulum_fa + parahippocampalcingulum_fa + corticospinalpyramidal_fa + anteriorthalamicradiations_fa + uncinate_fa + inferiorlongitudinalfasiculus_fa + inferiorfrontooccipitalfasiculus_fa + forcepsmajor_fa + forcepsminor_fa + corpuscallosum_fa + superiorlongitudinalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + parietalsuperiorlongitudinalfasiculus_fa + superiorcorticostriate_fa + superiorcorticostriatefrontalcortex_fa + superiorcorticostriateparietalcortex_fa + striatalinferiorfrontalcortex_fa + inferiorfrontalsuperiorfrontalcortex_fa + fornix_exfimbria_fa + fornix_md + cingulatecingulum_md + parahippocampalcingulum_md + corticospinalpyramidal_md + anteriorthalamicradiations_md + uncinate_md + inferiorlongitudinalfasiculus_md + inferiorfrontooccipitalfasiculus_md + forcepsmajor_md + forcepsminor_md + corpuscallosum_md + superiorlongitudinalfasiculus_md + temporalsuperiorlongitudinalfasiculus_md + parietalsuperiorlongitudinalfasiculus_md + superiorcorticostriate_md + superiorcorticostriatefrontalcortex_md + superiorcorticostriateparietalcortex_md + striatalinferiorfrontalcortex_md + inferiorfrontalsuperiorfrontalcortex_md + fornix_exfimbria_md + fornix_wmv + cingulatecingulum_wmv + parahippocampalcingulum_wmv + corticospinalpyramidal_wmv + anteriorthalamicradiations_wmv + uncinate_wmv + inferiorlongitudinalfasiculus_wmv + inferiorfrontooccipitalfasiculus_wmv + forcepsmajor_wmv + forcepsminor_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + temporalsuperiorlongitudinalfasiculus_wmv + parietalsuperiorlongitudinalfasiculus_wmv + superiorcorticostriate_wmv + superiorcorticostriatefrontalcortex_wmv + superiorcorticostriateparietalcortex_wmv + striatalinferiorfrontalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv + fornix_exfimbria_wmv'

fit_sem_GMWMlat_85_all_free <- sem(Model_GMWMlat_85_all_free, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMWMlat_85_all_free, fit.measures=TRUE,rsquare=T,standardized=T)

#Comparison between the model with freely estimated parameters and a model with constrained parameters for the grey & white matter metrics model with all the regions
Model_GMWMlat_85_all_constrained <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                     Cognitive_factor ~ a*bankssts_ct + a*caudalanteriorcingulate_ct + a*caudalmiddlefrontal_ct + a*cuneus_ct + a*entorhinal_ct + a*fusiform_ct + a*inferiorparietal_ct + a*inferiortemporal_ct + a*isthmuscingulate_ct + a*lateraloccipital_ct + a*lateralorbitofrontal_ct + a*lingual_ct + a*medialorbitofrontal_ct + a*middletemporal_ct + a*parahippocampal_ct + a*paracentral_ct + a*parsopercularis_ct + a*parsorbitalis_ct + a*parstriangularis_ct + a*pericalcarine_ct + a*postcentral_ct + a*posteriorcingulate_ct + a*precentral_ct + a*precuneus_ct + a*rostralanteriorcingulate_ct + a*rostralmiddlefrontal_ct + a*superiorfrontal_ct + a*superiorparietal_ct + a*superiortemporal_ct + a*supramarginal_ct + a*frontalpole_ct + a*temporalpole_ct + a*transversetemporal_ct + a*insula_ct + a*bankssts_sa + a*caudalanteriorcingulate_sa + a*caudalmiddlefrontal_sa + a*cuneus_sa + a*entorhinal_sa + a*fusiform_sa + a*inferiorparietal_sa + a*inferiortemporal_sa + a*isthmuscingulate_sa + a*lateraloccipital_sa + a*lateralorbitofrontal_sa + a*lingual_sa + a*medialorbitofrontal_sa + a*middletemporal_sa + a*parahippocampal_sa + a*paracentral_sa + a*parsopercularis_sa + a*parsorbitalis_sa + a*parstriangularis_sa + a*pericalcarine_sa + a*postcentral_sa + a*posteriorcingulate_sa + a*precentral_sa + a*precuneus_sa + a*rostralanteriorcingulate_sa + a*rostralmiddlefrontal_sa + a*superiorfrontal_sa + a*superiorparietal_sa + a*superiortemporal_sa + a*supramarginal_sa + a*frontalpole_sa + a*temporalpole_sa + a*transversetemporal_sa + a*insula_sa + a*bankssts_gmv + a*caudalanteriorcingulate_gmv + a*caudalmiddlefrontal_gmv + a*cuneus_gmv + a*entorhinal_gmv + a*fusiform_gmv + a*inferiorparietal_gmv + a*inferiortemporal_gmv + a*isthmuscingulate_gmv + a*lateraloccipital_gmv + a*lateralorbitofrontal_gmv + a*lingual_gmv + a*medialorbitofrontal_gmv + a*middletemporal_gmv + a*parahippocampal_gmv + a*paracentral_gmv + a*parsopercularis_gmv + a*parsorbitalis_gmv + a*parstriangularis_gmv + a*pericalcarine_gmv + a*postcentral_gmv + a*posteriorcingulate_gmv + a*precentral_gmv + a*precuneus_gmv + a*rostralanteriorcingulate_gmv + a*rostralmiddlefrontal_gmv + a*superiorfrontal_gmv + a*superiorparietal_gmv + a*superiortemporal_gmv + a*supramarginal_gmv + a*frontalpole_gmv + a*temporalpole_gmv + a*transversetemporal_gmv + a*insula_gmv + a*fornix_fa + a*cingulatecingulum_fa + a*parahippocampalcingulum_fa + a*corticospinalpyramidal_fa + a*anteriorthalamicradiations_fa + a*uncinate_fa + a*inferiorlongitudinalfasiculus_fa + a*inferiorfrontooccipitalfasiculus_fa + a*forcepsmajor_fa + a*forcepsminor_fa + a*corpuscallosum_fa + a*superiorlongitudinalfasiculus_fa + a*temporalsuperiorlongitudinalfasiculus_fa + a*parietalsuperiorlongitudinalfasiculus_fa + a*superiorcorticostriate_fa + a*superiorcorticostriatefrontalcortex_fa + a*superiorcorticostriateparietalcortex_fa + a*striatalinferiorfrontalcortex_fa + a*inferiorfrontalsuperiorfrontalcortex_fa + a*fornix_exfimbria_fa + a*fornix_md + a*cingulatecingulum_md + a*parahippocampalcingulum_md + a*corticospinalpyramidal_md + a*anteriorthalamicradiations_md + a*uncinate_md + a*inferiorlongitudinalfasiculus_md + a*inferiorfrontooccipitalfasiculus_md + a*forcepsmajor_md + a*forcepsminor_md + a*corpuscallosum_md + a*superiorlongitudinalfasiculus_md + a*temporalsuperiorlongitudinalfasiculus_md + a*parietalsuperiorlongitudinalfasiculus_md + a*superiorcorticostriate_md + a*superiorcorticostriatefrontalcortex_md + a*superiorcorticostriateparietalcortex_md + a*striatalinferiorfrontalcortex_md + a*inferiorfrontalsuperiorfrontalcortex_md + a*fornix_exfimbria_md + a*fornix_wmv + a*cingulatecingulum_wmv + a*parahippocampalcingulum_wmv + a*corticospinalpyramidal_wmv + a*anteriorthalamicradiations_wmv + a*uncinate_wmv + a*inferiorlongitudinalfasiculus_wmv + a*inferiorfrontooccipitalfasiculus_wmv + a*forcepsmajor_wmv + a*forcepsminor_wmv + a*corpuscallosum_wmv + a*superiorlongitudinalfasiculus_wmv + a*temporalsuperiorlongitudinalfasiculus_wmv + a*parietalsuperiorlongitudinalfasiculus_wmv + a*superiorcorticostriate_wmv + a*superiorcorticostriatefrontalcortex_wmv + a*superiorcorticostriateparietalcortex_wmv + a*striatalinferiorfrontalcortex_wmv + a*inferiorfrontalsuperiorfrontalcortex_wmv + a*fornix_exfimbria_wmv'

fit_sem_GMWMlat_85_all_constrained <- sem(Model_GMWMlat_85_all_constrained, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMWMlat_85_all_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_GMWMlat_85_all_free,fit_sem_GMWMlat_85_all_constrained)

#Create tables to save the parameters of the regularized model and the model with the regularized ROIs for white matter metrics
paramest_GMWMlat_all_85 <- data.frame(lhs=rep(NA, times=162),
                                        op=rep(NA, times=162),
                                        rhs=rep(NA, times=162),
                                        est=rep(NA, times=162),
                                        se=rep(NA, times=162),
                                        z=rep(NA, times=162),
                                        pvalue=rep(NA, times=162),
                                        ci.lower=rep(NA, times=162),
                                        ci.upper=rep(NA, times=162),
                                        std.lv=rep(NA, times=162),
                                        std.all=rep(NA, times=162),
                                        std.nox=rep(NA, times=162))

paramest_GMWMlat_all_85[,] <-parameterEstimates(fit_sem_GMWMlat_85_all_free,standardized = TRUE, rsquare=TRUE)[6:167,]
write.csv(paramest_GMWMlat_all_85, "paramest_GMWMlat_GMWMall_85.csv", row.names=FALSE)

paramest_GMWMlat_reg_85 <- data.frame(lhs=rep(NA, times=50),
                                      op=rep(NA, times=50),
                                      rhs=rep(NA, times=50),
                                      est=rep(NA, times=50),
                                      se=rep(NA, times=50),
                                      z=rep(NA, times=50),
                                      pvalue=rep(NA, times=50),
                                      ci.lower=rep(NA, times=50),
                                      ci.upper=rep(NA, times=50),
                                      std.lv=rep(NA, times=50),
                                      std.all=rep(NA, times=50),
                                      std.nox=rep(NA, times=50))

paramest_GMWMlat_reg_85[,] <-parameterEstimates(fit_sem_GMWMlat_85_reg_free,standardized = TRUE, rsquare=TRUE)[6:55,]
write.csv(paramest_GMWMlat_reg_85, "paramest_GMWMlat_reg_T085.csv", row.names=FALSE)


#-------------------------------------------------------------------------------------------------------------------------------
#Sample 85%: Compare the regularized models with grey and white matter
#-------------------------------------------------------------------------------------------------------------------------------

fit_vector_85<-data.frame(label=c("Model Grey & White Matter","Model Grey Matter","Model White Matter"),
                            AIC=rep(NA,times=3),
                            BIC=rep(NA,times=3))
#Grey & White Matter
summary(fit_sem_GMWMlat_85_reg_free, fit.measures=TRUE,rsquare=T,standardized=T)
fit_vector_85[1,2]<-AIC(fit_sem_GMWMlat_85_reg_free)
fit_vector_85[1,3]<-BIC(fit_sem_GMWMlat_85_reg_free)

#'I don't need the white matter if I already have the grey matter'
Model_GMWM_WMconstrainedlat_85 <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ caudalmiddlefrontal_ct + fusiform_ct + lateralorbitofrontal_ct + medialorbitofrontal_ct + middletemporal_ct + parahippocampal_ct + paracentral_ct + parsorbitalis_ct + pericalcarine_ct + precuneus_ct + rostralanteriorcingulate_ct + rostralmiddlefrontal_ct + superiorfrontal_ct + superiortemporal_ct + insula_ct + bankssts_sa + cuneus_sa + entorhinal_sa + lateraloccipital_sa + parsorbitalis_sa + postcentral_sa + precentral_sa + superiortemporal_sa + fusiform_gmv + inferiorparietal_gmv + parahippocampal_gmv + paracentral_gmv + pericalcarine_gmv + precentral_gmv + precuneus_gmv + rostralmiddlefrontal_gmv + supramarginal_gmv + 
                                                      0*cingulatecingulum_fa + 0*parahippocampalcingulum_fa + 0*uncinate_fa + 0*inferiorfrontooccipitalfasiculus_fa + 0*temporalsuperiorlongitudinalfasiculus_fa + 0*superiorcorticostriateparietalcortex_fa + 0*fornix_exfimbria_fa + 0*parahippocampalcingulum_md + 0*inferiorlongitudinalfasiculus_md + 0*forcepsmajor_md + 0*temporalsuperiorlongitudinalfasiculus_md + 0*striatalinferiorfrontalcortex_md + 0*anteriorthalamicradiations_wmv + 0*inferiorlongitudinalfasiculus_wmv + 0*corpuscallosum_wmv + 0*superiorlongitudinalfasiculus_wmv + 0*superiorcorticostriateparietalcortex_wmv + 0*inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_GMWM_WMconstrainedlat_85 <- sem(Model_GMWM_WMconstrainedlat_85, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMWM_WMconstrainedlat_85, fit.measures=TRUE,rsquare=T,standardized=T)
fit_vector_85[2,2]<-AIC(fit_sem_GMWM_WMconstrainedlat_85)
fit_vector_85[2,3]<-BIC(fit_sem_GMWM_WMconstrainedlat_85)

anova(fit_sem_GMWMlat_85_reg_free,fit_sem_GMWM_WMconstrainedlat_85)

#'I don't need the grey matter if I already have the white matter'
Model_GMWM_GMconstrainedlat_85 <- 'Cognitive_factor =~ Picture_Vocabulary + Flanker + Oral_Reading_Recognition + Rey_Auditory_Verbal_Learning + Little_Man
                                   Cognitive_factor ~ 0*caudalmiddlefrontal_ct + 0*fusiform_ct + 0*lateralorbitofrontal_ct + 0*medialorbitofrontal_ct + 0*middletemporal_ct + 0*parahippocampal_ct + 0*paracentral_ct + 0*parsorbitalis_ct + 0*pericalcarine_ct + 0*precuneus_ct + 0*rostralanteriorcingulate_ct + 0*rostralmiddlefrontal_ct + 0*superiorfrontal_ct + 0*superiortemporal_ct + 0*insula_ct + 0*bankssts_sa + 0*cuneus_sa + 0*entorhinal_sa + 0*lateraloccipital_sa + 0*parsorbitalis_sa + 0*postcentral_sa + 0*precentral_sa + 0*superiortemporal_sa + 0*fusiform_gmv + 0*inferiorparietal_gmv + 0*parahippocampal_gmv + 0*paracentral_gmv + 0*pericalcarine_gmv + 0*precentral_gmv + 0*precuneus_gmv + 0*rostralmiddlefrontal_gmv + 0*supramarginal_gmv +
                                                      cingulatecingulum_fa + parahippocampalcingulum_fa + uncinate_fa + inferiorfrontooccipitalfasiculus_fa + temporalsuperiorlongitudinalfasiculus_fa + superiorcorticostriateparietalcortex_fa + fornix_exfimbria_fa + parahippocampalcingulum_md + inferiorlongitudinalfasiculus_md + forcepsmajor_md + temporalsuperiorlongitudinalfasiculus_md + striatalinferiorfrontalcortex_md + anteriorthalamicradiations_wmv + inferiorlongitudinalfasiculus_wmv + corpuscallosum_wmv + superiorlongitudinalfasiculus_wmv + superiorcorticostriateparietalcortex_wmv + inferiorfrontalsuperiorfrontalcortex_wmv'

fit_sem_GMWM_GMconstrainedlat_85 <- sem(Model_GMWM_GMconstrainedlat_85, data=data_total_subsample85,estimator="mlr",missing="fiml")
summary(fit_sem_GMWM_GMconstrainedlat_85, fit.measures=TRUE,rsquare=T,standardized=T)
fit_vector_85[3,2]<-AIC(fit_sem_GMWM_GMconstrainedlat_85)
fit_vector_85[3,3]<-BIC(fit_sem_GMWM_GMconstrainedlat_85)

anova(fit_sem_GMWMlat_85_reg_free,fit_sem_GMWM_GMconstrainedlat_85)

#'Grey and white matter have equally sized, complementary effects'
summary(fit_sem_GMWMlat_85_reg_constrained, fit.measures=TRUE,rsquare=T,standardized=T)

anova(fit_sem_GMWMlat_85_reg_free,fit_sem_GMWMlat_85_reg_constrained)

#Data Visualization for the comparison between the grey matter model, the white matter model and the model with grey & white matter
#Create a data.frame with model names
AIC_object_85 <- data.frame(Model = fit_vector_85[,1],
                            AIC = factor(rep('AIC', each = 3)),
                            values = fit_vector_85[,2])

BIC_object_85 <- data.frame(Model = fit_vector_85[,1],
                            BIC = factor(rep('BIC', each = 3)),
                            values = fit_vector_85[,3])

#Plot the 'raw' AIC's
plot_fitraw_aic_85<-ggplot(AIC_object_85, aes(values, AIC, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1),colour="black") +
  coord_cartesian(xlim = c(122800, 123300)) +
  xlab('Information criterion (AIC)') +
  ylab('') +
  theme_grey(base_size = 30) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
        axis.line = element_line(colour = "black"))+  # adding a black line for x and y axis
  geom_text(aes(label=c("Grey & White\nMatter Model","Grey Matter Model","White Matter Model")), size = 6, color="white",position = position_dodge(width = 1),hjust = 1.1) +
  theme(legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_fill_manual("legend", values = c("Model Grey & White Matter" = "black", "Model Grey Matter" = "coral1", "Model White Matter" = "darkcyan"))+
  annotate("text", x = 123273, y = 1.33, label = "123253.3",size=6)+
  annotate("text", x = 123184, y = 1, label = "123164.3",size=6)+
  annotate("text", x = 122975, y = 0.66, label = "122955.0",size=6)
ggsave("plot_fitraw_aic85.png")

#Plot the 'raw' BIC's
plot_fitraw_bic_85<-ggplot(BIC_object_85, aes(values, BIC, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1),colour="black") +
  coord_cartesian(xlim = c(123380, 123510)) +
  xlab('Information criterion (BIC)') +
  ylab('') +
  theme_grey(base_size = 30) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
        axis.line = element_line(colour = "black"))+  # adding a black line for x and y axis
  geom_text(aes(label=c("Grey & White\nMatter Model","Grey Matter Model","White Matter Model")), size = 6, color="white",position = position_dodge(width = 1),hjust = 1.1) +
  theme(legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_fill_manual("legend", values = c("Model Grey & White Matter" = "black", "Model Grey Matter" = "coral1", "Model White Matter" = "darkcyan"))+
  annotate("text", x = 123495, y = 1.33, label = "123489.4",size=6)+
  annotate("text", x = 123506, y = 1, label = "123500.6",size=6)+
  annotate("text", x = 123426, y = 0.66, label = "123420.1",size=6)
ggsave("plot_fitraw_bic85.png")  

#Plot the adjusted Rsquare
Rsquare_all_85 <- data.frame(Model = c("Model Grey Matter","Model White Matter","Model Grey & White Matter"),
                             Rsquare = factor(rep('Rsquare', each = 3)),
                             Values_Rsquare = c(0.127,0.112,0.159),
                             Values_Rsquare.adjusted=c(0.123,0.109,0.154))
#adj(R^2) = 1 - (1-R^2) * (n - 1)/(n - p - 1)
Rsquare_all_85$Model <- factor(Rsquare_all_85$Model,levels = c("Model Grey Matter","Model White Matter","Model Grey & White Matter"))


plot_rsquare_adj_85<-ggplot(Rsquare_all_85,aes(Model,Values_Rsquare.adjusted,fill=Model))+
  geom_col(color="black") +
  coord_cartesian(ylim = c(0.1, 0.16)) +
  xlab(' ') +
  ylab('Adjusted Rsquare') +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
        axis.line = element_line(colour = "black")) +  # adding a black line for x and y axis
  geom_text(aes(label=c("12.3%","10.9%","15.4%")), size = 8, color="black",position = position_dodge(width = 0.5),hjust = 0.5,vjust=-0.5) +
  theme(legend.position="none",axis.ticks.x=element_blank(),axis.title=element_text(size=20),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 15)) +
  scale_x_discrete(breaks=c("Model Grey Matter", "Model White Matter", "Model Grey & White Matter"),
                   labels=c("Model Grey \n Matter", "Model White \n Matter", "Model Grey & \n White Matter")) +
  scale_fill_manual(breaks=c("Model Grey Matter", "Model White Matter", "Model Grey & White Matter"),values = c("coral1","darkcyan","black"))
ggsave("plot_rsquare_adj_85.png")  

#-------------------------------------------------------------------------------------------------------------------------------



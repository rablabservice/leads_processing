library(xlsx)
library(psych)
library(outliers)
library(reshape2)
library(OptimalCutpoints)
library(pROC)
library(statquotes)
library(ggplot2)
library(cutpointr)
library(ggpubr)
library(leonaRdo)
library(lubridate)
library(flexdashboard)
library(plotly)
library(ggseg)
library(knitr)
library(rmarkdown)
library(ggthemes)
library(rcompanion)
library(effsize)
library(sjstats)

setwd("/Volumes/petcore/Projects/LEADS/data_f7p1/LONI_uploads/")

#### Loading FBB and FTP metadata and extraction values ####

fbbextfiles=file.info(list.files("service/",pattern="FBB_ROI_Extraction*", full.names = T))
fbb_ext=read.csv(rownames(fbbextfiles)[order(fbbextfiles$mtime)][nrow(fbbextfiles)])
ftpextfiles=file.info(list.files("service/",pattern="FTP_ROI_Extraction*", full.names = T))
ftp_ext=read.csv(rownames(ftpextfiles)[order(ftpextfiles$mtime)][nrow(ftpextfiles)])
fbbmetafiles=file.info(list.files("service/",pattern="LEADS_ALLFBBs*", full.names = T))
fbb_meta=read.csv(rownames(fbbmetafiles)[order(fbbmetafiles$mtime)][nrow(fbbmetafiles)])
ftpmetafiles=file.info(list.files("service/",pattern="LEADS_ALLFTPs*", full.names = T))
ftp_meta=read.csv(rownames(ftpmetafiles)[order(ftpmetafiles$mtime)][nrow(ftpmetafiles)])

#### Adding FDG data ####

fdgextfiles=file.info(list.files("service/",pattern="FDG_ROI_Extraction*", full.names = T))
fdg_ext=read.csv(rownames(fdgextfiles)[order(fdgextfiles$mtime)][nrow(fdgextfiles)])
fdgmetafiles=file.info(list.files("service/",pattern="LEADS_ALLFDGs*", full.names = T))
fdg_meta=read.csv(rownames(fdgmetafiles)[order(fdgmetafiles$mtime)][nrow(fdgmetafiles)])

#### Loading File with cohort assignments from EDC (only CN/PT) ####

# Updated June 15th 2020 to use the Participants Registration file from EDC (CRF data)
idsedcfiles=file.info(list.files("service/",pattern="Participa*", full.names = T))
idsedc=read.csv(rownames(idsedcfiles)[order(idsedcfiles$mtime)][nrow(idsedcfiles)])
idsedc=idsedc[,c("subject.label","dd_revision_field.value")]
colnames(idsedc)=c("ID","Cohort")
idsedc$Cohort=ifelse(idsedc$Cohort==1,"CN","Symptomatic") ## keep it for later

# Edited 9/26/23 by DRS - I removed columns from this spreadsheet that
# aren't kept up-to-date
#### Loading and working on Mayo MRI QC data  ####
mrimetafiles=file.info(list.files("service/",pattern="Mayo_*", full.names = T))
mri_meta=read.csv(rownames(mrimetafiles)[order(mrimetafiles$mtime)][nrow(mrimetafiles)])
mri_meta=subset(mri_meta, (grepl("MPRAGE", mri_meta$series_description) | grepl("IR-FSPGR", mri_meta$series_description)))
mri_included=subset(mri_meta, !mri_meta$release_from_quarantine==0)
mri_included=mri_included[1:2]
colnames(mri_included)=c("ID","MRI_Date")
mri_included$MRI_Date=as.Date(as.character(mri_included$MRI_Date, format = "%m/%d/%Y"))
mri_included$MRI_ScanIdentifier=paste(mri_included$ID,"_",mri_included$MRI_Date,sep="")
mri_included=mri_included[!duplicated(mri_included$MRI_ScanIdentifier),]
mri_included=subset(mri_included, !grepl("LDP",mri_included$ID))

#### Loading updated PETCore amyloid-PET screening file  ####

fbb_pt_reads=xlsx::read.xlsx("service/LEADS_Internal_PET-Screening.xlsx",1,header=T)
fbb_missing_reads=subset(fbb_pt_reads, is.na(fbb_pt_reads$Final.Read))
fbb_pt_reads=subset(fbb_pt_reads, !is.na(fbb_pt_reads$Final.Read))
fbb_pt_reads$Cohort=trimws(fbb_pt_reads$Cohort, which=c("both")) # remove leading and trailing spaces to make functions work in the script
fbb_pt_reads$disagreement=trimws(fbb_pt_reads$disagreement, which=c("both")) # remove weird spaces that maybe make this cell not an NA


#### Harmonizing ID+Date identifiers and dates across files  ####

fbb_ext$ScanIdentifier=paste(fbb_ext$ID,"_",fbb_ext$FBBPET_Date,sep="")
ftp_ext$ScanIdentifier=paste(ftp_ext$ID,"_",ftp_ext$FTPPET_Date,sep="")
fbb_pt_reads$ScanIdentifier=paste(fbb_pt_reads$ID,"_",fbb_pt_reads$FBB.Scan.Date,sep="")

fbb_meta$Acq.Date=as.Date(fbb_meta$Acq.Date, format = "%m/%d/%Y")
ftp_meta$Acq.Date=as.Date(ftp_meta$Acq.Date, format = "%m/%d/%Y")

fbb_meta$ScanIdentifier=paste(fbb_meta$Subject,"_",fbb_meta$Acq.Date,sep="")
ftp_meta$ScanIdentifier=paste(ftp_meta$Subject,"_",ftp_meta$Acq.Date,sep="")

#### Doing the same for FDG-PET now ####

fdg_ext$ScanIdentifier=paste(fdg_ext$ID,"_",fdg_ext$FDGPET_Date,sep="")
fdg_meta$Acq.Date=as.Date(fdg_meta$Acq.Date, format = "%m/%d/%Y")
fdg_meta$ScanIdentifier=paste(fdg_meta$Subject,"_",fdg_meta$Acq.Date,sep="")

#### Recode and clean data in PETCore screening file  ####

fbb_pt_reads=fbb_pt_reads[c("ID","FBB.Scan.Date","ScanIdentifier","SUVR", "Final.Read","disagreement","Cohort","BAPL")]
colnames(fbb_pt_reads)[c(2,4,5,6)]=c("FBBPET_Date","Screening_PETONLY_Composite_SUVR","Screening_PETONLY_Final_Read","Screening_PETONLY_Disagreement")
fbb_pt_reads$Screening_PETONLY_Final_Read=ifelse(grepl("pos",fbb_pt_reads$Screening_PETONLY_Final_Read)=="TRUE","1","0")
fbb_pt_reads$Screening_PETONLY_Disagreement=ifelse(is.na(fbb_pt_reads$Screening_PETONLY_Disagreement)==FALSE,"1",0)
fbb_pt_reads=fbb_pt_reads[c(1,3:ncol(fbb_pt_reads))]
fbb_pt_reads$Screening_PETONLY_AmyPos_Quantification_1p18=ifelse(fbb_pt_reads$Screening_PETONLY_Composite_SUVR>1.18,"1","0")
fbb_pt_reads$Screening_PETONLY_VisualRead=ifelse(fbb_pt_reads$BAPL>1,"1","0") # BAPL scored is used to call a scan visually pos or neg
fbb_pt_reads=fbb_pt_reads[c(1:6,8:9)]

fbb_pt_reads=fbb_pt_reads[c(2,1,3,7,8,5,4,6)]

#### Merge final FBB-PET database ####

fbb_meta=fbb_meta[c(1,ncol(fbb_meta))]
colnames(fbb_meta)[1]=c("ImageID")
fbb_ext=fbb_ext[c(1:3,6:ncol(fbb_ext))]
fbb_meta_ext=merge(fbb_meta, fbb_ext, by="ScanIdentifier")
fbbdf=merge(fbb_pt_reads, fbb_meta_ext, by="ID",all.y = TRUE)

check_missing_cohort=subset(fbbdf, fbbdf$ID %in% fbb_missing_reads$ID)

if (nrow(check_missing_cohort)>0) {
  
  message("These subjects were included in the processing spreadsheet and in the internal spreadsheet, but are missing the read")
  print(check_missing_cohort[c(1,9:11)])
  message("I will exclude them from all the files. Press [enter] to acknowledge and continue")
  line <- readline()
  fbbdf=subset(fbbdf, !(fbbdf$ID %in% fbb_missing_reads$ID)) # New line to not crash when a subject was processed (and so is in the extraction files) but has not been read yet
  
}

fbbdf$Cohort=ifelse(is.na(fbbdf$Cohort),"CN",as.character(fbbdf$Cohort))
fbbdf$Cohort=factor(fbbdf$Cohort, levels=c("EOAD","EOnonAD","CN"))

#### Subset FBB-PET database based on MRI passing QC from Mayo ####

fbbdf$MRI_ScanIdentifier=paste(fbbdf$ID,"_",fbbdf$MRI_Date,sep="")
excluded_fbb_mriqc=subset(fbbdf, !(MRI_ScanIdentifier %in% mri_included$MRI_ScanIdentifier))
colnames(excluded_fbb_mriqc)[8]="CohortAssgn"
fname1=paste("logs/FBB_PETCore_Excluded_MRIQC_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
write.csv(excluded_fbb_mriqc,fname1, row.names = F, na="")
fbbdf=subset(fbbdf, MRI_ScanIdentifier %in% mri_included$MRI_ScanIdentifier)

#### Additional QC on FBB-PET database: remove rows with NA for Global score ####

failed_globalsuvrs=length(which(is.na(fbbdf$GlobalSUVR)))

if (failed_globalsuvrs>0) {
  failed_suvrs=subset(fbbdf, is.na(fbbdf$GlobalSUVR))
  fname2=paste("logs/FBB_PETCore_FailedGlobalSUVR_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(failed_suvrs,fname2, row.names = F, na="")
  message("**********************")
  message(paste("Global SUVR failed to be calculated for",failed_globalsuvrs,"subject(s).","Check log file for details!"))
  message("**********************")
  fbbdf=subset(fbbdf, !is.na(fbbdf$GlobalSUVR))
}

#### Additional QC on FBB-PET database: remove rows with NA for Global score Type 2 ####

failed_globalsuvrs_t2=length(which(is.na(fbbdf$GlobalSUVR_Type2)))

if (failed_globalsuvrs_t2>0) {
  failed_suvrs=subset(fbbdf, is.na(fbbdf$GlobalSUVR_Type2))
  fname2=paste("logs/FBB_PETCore_FailedGlobalSUVR_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(failed_suvrs,fname2, row.names = F, na="")
  message("**********************")
  message(paste("Global SUVR Type 2 failed to be calculated for",failed_globalsuvrs_t2,"subject(s).","Check log file for details!"))
  message("**********************")
  fbbdf=subset(fbbdf, !is.na(fbbdf$GlobalSUVR))
}

#### Additional QC on FBB-PET database: Dynamic check on outliers in reference region size ####

outlier_refregsize=fbbdf[c("ID","ScanIdentifier.y","ScalingFactor_WholeCereb_ClustSize")]
outlier_refregsize$ScalingFactor_WholeCereb_ClustSize=scale(outlier_refregsize$ScalingFactor_WholeCereb_ClustSize)
outlier_refregsize=subset(outlier_refregsize, abs(outlier_refregsize$ScalingFactor_WholeCereb_ClustSize)>3)

# if (nrow(outlier_refregsize)>0) {
#   message("**********************")
#   message(paste("The cluster size for the Whole Cerebellum was >3SD from the mean for",nrow(outlier_refregsize),"subject(s).","Check log file for details!"))
#   message("**********************")
#   for (i in 1:nrow(outlier_refregsize)) {
#     tempid=outlier_refregsize[i,2]
#     dec<-readline(prompt=paste("Should I keep",tempid,"in the database? Please answer y/n: "))
#     outlier_refregsize[i,4]=dec
# 
#     if (dec=="n") {
# 
#       fbbdf=subset(fbbdf, !fbbdf$ScanIdentifier.y==tempid)
# 
#     }
# 
#   }
# 
#   colnames(outlier_refregsize)[4]="Present_in_LONI_database"
#   fname3=paste("logs/FBB_PETCore_OutlierWholeCBL_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
#   write.csv(outlier_refregsize,fname3, row.names = F, na="")
# 
# }

#### Sort final FBB-PET database, add centiloids ####

fbbdf=fbbdf[c(1,11,10,3:8,13:(ncol(fbbdf)-1))]
colnames(fbbdf)[c(11,12)]=c("MRIBASED_Composite_SUVR_Type2","MRIBASED_Composite_SUVR")

fbbdf$ScalingFactor_WholeCereb=round(fbbdf$ScalingFactor_WholeCereb, 4)

colnames(fbbdf)[13:130]=paste(names(fbbdf[13:130]),"_MRIBASED_SUVR",sep="")

fbbdf$MRIBASED_Composite_Centiloids=(157.15*fbbdf$MRIBASED_Composite_SUVR_Type2)-151.87 # New formula for Freesurfer 7
fbbdf=fbbdf[c(1:12,250,13:249)]

#### Additional QC: Check that broad cohort assignment CN and Symptomatic fits with EDC data ####

idsedc=subset(idsedc, ID %in% fbbdf$ID)
idsedc=merge(idsedc, fbbdf[c("ID","Cohort")], by="ID")
idsedc$Cohort.QC=ifelse(idsedc$Cohort.y=="CN","CN","Symptomatic")

cohortqc=all(idsedc$Cohort.x == idsedc$Cohort.QC)
idsedc_sbs=subset(idsedc, !idsedc$Cohort.x==idsedc$Cohort.QC)
idsedc_sbs=idsedc_sbs[c(1,2,4)]
colnames(idsedc_sbs)=c("ID","Cohort_EDC","Cohort_ReadsFile")

#### Success for QC: save database ####

if (cohortqc==TRUE) {

  colnames(fbbdf)[9]="CohortAssgn"
  
  ## Small module to merge in the longit reference extraction data ##
  
  fbbextfiles_compwm=file.info(list.files("service/",pattern="FBB_CompWM_Extraction*", full.names = T))
  fbb_ext_compwm=read.csv(rownames(fbbextfiles_compwm)[order(fbbextfiles_compwm$mtime)][nrow(fbbextfiles_compwm)])
  fbb_ext_compwm$ScanIdentifier=paste(fbb_ext_compwm$ID,"_",fbb_ext_compwm$FBBPET_Date,sep="")
  
  #### Additional QC on FBB-PET Longitudinal extraction database: Dynamic check on outliers in reference region size ####
  
  outlier_refregsize_compwm=fbb_ext_compwm[c("ID","ScanIdentifier","ScalingFactor_CompWM_ClustSize")]
  outlier_refregsize_compwm$ScalingFactor_CompWM_ClustSize=scale(outlier_refregsize_compwm$ScalingFactor_CompWM_ClustSize)
  outlier_refregsize_compwm=subset(outlier_refregsize_compwm, abs(outlier_refregsize_compwm$ScalingFactor_CompWM_ClustSize)>3)
  
  # if (nrow(outlier_refregsize_compwm)>0) {
  #   message("**********************")
  #   message(paste("The cluster size for the Composite WM was >3SD from the mean for",nrow(outlier_refregsize_compwm),"subject(s).","Check log file for details!"))
  #   message("**********************")
  #   
  #   for (i in 1:nrow(outlier_refregsize_compwm)) {
  #     tempid=outlier_refregsize_compwm[i,2]
  #     dec<-readline(prompt=paste("Should I keep",tempid,"in the database? Please answer y/n: "))
  #     outlier_refregsize_compwm[i,4]=dec
  #     
  #     if (dec=="n") {
  #       
  #       fbb_ext_compwm=subset(fbb_ext_compwm, !fbb_ext_compwm$ScanIdentifier==tempid)
  #       
  #     }
  #     
  #   }
  #   
  #   colnames(outlier_refregsize_compwm)[4]="Present_in_LONI_database"
  #   fname_compwmFBB=paste("logs/FBB_PETCore_CompositeWM_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  #   write.csv(outlier_refregsize_compwm,fname_compwmFBB, row.names = F, na="")
  #   
  # }
  
  
  ### we are theoretically ready to merge in the final fbbdf
  
  fbb_ext_compwm=fbb_ext_compwm[c(1:2,6:8)]
  colnames(fbb_ext_compwm)[c(3:5)]=c("ScalingFactor_CompositeWM","CompositeWM_MRIBASED_SUVR","ScalingFactor_CompositeWM_ClustSize")
  fbb_ext_compwm$ScalingFactor_CompositeWM=round(fbb_ext_compwm$ScalingFactor_CompositeWM, 4)
  fbbdf=merge(fbbdf, fbb_ext_compwm, by=c("ID","FBBPET_Date"))
  fbbdf=fbbdf[c(1:10,251,11:18,252,19:132,253,133:137,253,138:(ncol(fbbdf)-3))]
  colnames(fbbdf)[141]="CompositeWM_ClustSize"
  
  ### Done, ready to export!
  ##

  fname=paste("ready_datasets/FBB_PETCore_Analysis_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(fbbdf,fname, row.names = F, na="")

  #### Module for FBB-PET Quarterly report: To be done with report from MRICore as input as well ####
  #### To be modified according to Freesurfer 7 new data, e.g. one more column in FBB dataset ####
  # load report for MRI sent by Jessica Collins - Now Alexandra Touratouglu
  
  mrireportfiles=file.info(list.files("service/",pattern="MGH_LEADS*", full.names = T))
  mrireport=read.csv(rownames(mrireportfiles)[order(mrireportfiles$mtime)][nrow(mrireportfiles)])
  mrireport=mrireport[c("subject_code","Acq.Date")]
  colnames(mrireport)=c("ID","MRI_Date")
  mrireport$MRI_Date=as.Date(mrireport$MRI_Date)

  fbbdf_qrtrep=merge(fbb_ext[1:3],fbbdf,by=c("ID","FBBPET_Date"))
  fbbdf_qrtrep$MRIident=paste(fbbdf_qrtrep$ID,"_",fbbdf_qrtrep$MRI_Date,sep="")
  mrireport$MRIident=paste(mrireport$ID,"_",mrireport$MRI_Date,sep="")

  fbbdf_qrtrep=subset(fbbdf_qrtrep, MRIident %in% mrireport$MRIident)
  fbbdf_qrtrep=fbbdf_qrtrep[c(1,2,4:13,15:255)] # modified to only include "new" amyloid computation

  colnames(fbbdf_qrtrep)[12]="MRIBASED_Composite_SUVR"

  fname2=paste("ready_datasets/QuarterlyReport",format(Sys.time(), "%m%Y"),"_FBB_PETCore_Analysis_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(fbbdf_qrtrep,fname2, row.names = F, na="")

  #### Preparing FTP-PET database ####

  ftp_meta=ftp_meta[c(1,ncol(ftp_meta))]
  colnames(ftp_meta)[1]=c("ImageID")
  ftp_ext=ftp_ext[c(1:3,6:ncol(ftp_ext))]
  ftp_meta_ext=merge(ftp_meta, ftp_ext, by="ScanIdentifier")
  ptcohort=fbbdf[c("ID","CohortAssgn")]
  ptcohort=ptcohort[!duplicated(ptcohort$ID),]
  ftpdf=merge(ptcohort, ftp_meta_ext, by="ID",all.y = TRUE)
  ftpdf$CohortAssgn=ifelse(is.na(ftpdf$CohortAssgn),"NotFound_CheckEDC",as.character(ftpdf$CohortAssgn))
  ftpdf$CohortAssgn=factor(ftpdf$CohortAssgn, levels=c("EOAD","EOnonAD","CN","NotFound_CheckEDC"))

  #### Subset FTP-PET database based on MRI passing QC from Mayo ####

  ftpdf$MRI_ScanIdentifier=paste(ftpdf$ID,"_",ftpdf$MRI_Date,sep="")
  excluded_ftp_mriqc=subset(ftpdf, !(MRI_ScanIdentifier %in% mri_included$MRI_ScanIdentifier))
  fname1=paste("logs/FTP_PETCore_Excluded_MRIQC_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(excluded_ftp_mriqc,fname1, row.names = F, na="")
  ftpdf=subset(ftpdf, MRI_ScanIdentifier %in% mri_included$MRI_ScanIdentifier)
  
  #### Additional QC on FTP-PET database: remove rows with NA for metaROI and Braak regions ####

  failed_metaroi=length(which(is.na(ftpdf$MetaROI)))

  if (failed_metaroi>0) {
    failed_mroi=subset(ftpdf, is.na(ftpdf$MetaROI))
    fname2=paste("logs/FTP_PETCore_FailedMetaROI_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
    write.csv(failed_mroi,fname2, row.names = F, na="")
    message("**********************")
    message(paste("MetaROI SUVR failed to be calculated for",failed_metaroi,"subject(s).","Check log file for details!"))
    message("**********************")
    ftpdf=subset(ftpdf, !is.na(ftpdf$MetaROI))
  }

  failed_braak12=length(which(is.na(ftpdf$Braak_12)))

  if (failed_braak12>0) {
    failed_b12roi=subset(ftpdf, is.na(ftpdf$Braak_12))
    fname2=paste("logs/FTP_PETCore_FailedBraak12_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
    write.csv(failed_b12roi,fname2, row.names = F, na="")
    message("**********************")
    message(paste("Braak_12 SUVR failed to be calculated for",failed_braak12,"subject(s).","Check log file for details!"))
    message("**********************")
    ftpdf=subset(ftpdf, !is.na(ftpdf$Braak_12))
  }

  failed_braak34=length(which(is.na(ftpdf$Braak_34)))

  if (failed_braak34>0) {
    failed_b34roi=subset(ftpdf, is.na(ftpdf$Braak_34))
    fname2=paste("logs/FTP_PETCore_FailedBraak34_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
    write.csv(failed_b34roi,fname2, row.names = F, na="")
    message("**********************")
    message(paste("Braak_34 SUVR failed to be calculated for",failed_braak34,"subject(s).","Check log file for details!"))
    message("**********************")
    ftpdf=subset(ftpdf, !is.na(ftpdf$Braak_34))
  }

  failed_braak56=length(which(is.na(ftpdf$Braak_56)))

  if (failed_braak56>0) {
    failed_b56roi=subset(ftpdf, is.na(ftpdf$Braak_56))
    fname2=paste("logs/FTP_PETCore_FailedBraak56_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
    write.csv(failed_b56roi,fname2, row.names = F, na="")
    message("**********************")
    message(paste("Braak_56 SUVR failed to be calculated for",failed_braak56,"subject(s).","Check log file for details!"))
    message("**********************")
    ftpdf=subset(ftpdf, !is.na(ftpdf$Braak_56))
  }


  #### Additional QC on FTP-PET database: Dynamic check on outliers in reference region size ####

  outlier_refregsize=ftpdf[c("ID","ScanIdentifier","ScalingFactor_InfCerebGray_ClustSize")]
  outlier_refregsize$ScalingFactor_InfCerebGray_ClustSize=scale(outlier_refregsize$ScalingFactor_InfCerebGray_ClustSize)
  outlier_refregsize=subset(outlier_refregsize, abs(outlier_refregsize$ScalingFactor_InfCerebGray_ClustSize)>3)

  # if (nrow(outlier_refregsize)>0) {
  #   message("**********************")
  #   message(paste("The cluster size for the Inferior Cerebellar GM  was >3SD from the mean for",nrow(outlier_refregsize),"subject(s).","Check log file for details!"))
  #   message("**********************")
  #   for (i in 1:nrow(outlier_refregsize)) {
  #     tempid=outlier_refregsize[i,2]
  #     dec<-readline(prompt=paste("Should I keep",tempid,"in the database? Please answer y/n: "))
  #     outlier_refregsize[i,4]=dec
  # 
  #     if (dec=="n") {
  # 
  #       ftpdf=subset(ftpdf, !ftpdf$ScanIdentifier==tempid)
  # 
  #     }
  # 
  #   }
  # 
  #   colnames(outlier_refregsize)[4]="Present_in_LONI_database"
  #   fname3=paste("logs/FTP_PETCore_OutlierinfCBLGM_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  #   write.csv(outlier_refregsize,fname3, row.names = F, na="")
  # 
  # }

  #### Sort final FTP-PET database and save it ####

  ftpdf=ftpdf[c(1,5,4,2,7:(ncol(ftpdf)-1))]
  ftpdf$ScalingFactor_InfCerebGray=round(ftpdf$ScalingFactor_InfCerebGray, 4)
  colnames(ftpdf)[c(9,128)]=c("Braak_1","Braak_1_ClustSize")
  ftpdf=ftpdf[-8]
  colnames(ftpdf)[7:124]=paste(names(ftpdf)[7:124],"_MRIBASED_SUVR",sep="")
  colnames(ftpdf)[6]="Assigned_MRIBASED_MetaROI_ADNIcutoff_1p2"
  
  ## Small module to add the longitudinal extractions
  
  ftpextfiles_erodedwm=file.info(list.files("service/",pattern="FTP_ErodedWM_Extraction*", full.names = T))
  ftp_ext_erodedwm=read.csv(rownames(ftpextfiles_erodedwm)[order(ftpextfiles_erodedwm$mtime)][nrow(ftpextfiles_erodedwm)])
  ftp_ext_erodedwm$ScanIdentifier=paste(ftp_ext_erodedwm$ID,"_",ftp_ext_erodedwm$FTPPET_Date,sep="")
  
  #### Additional QC on FTP-PET Longitudinal extractionn database: Dynamic check on outliers in reference region size ####
  
  outlier_refregsize_erodedwm=ftp_ext_erodedwm[c("ID","ScanIdentifier","ScalingFactor_ErodedWM_ClustSize")]
  outlier_refregsize_erodedwm$ScalingFactor_ErodedWM_ClustSize=scale(outlier_refregsize_erodedwm$ScalingFactor_ErodedWM_ClustSize)
  outlier_refregsize_erodedwm=subset(outlier_refregsize_erodedwm, abs(outlier_refregsize_erodedwm$ScalingFactor_ErodedWM_ClustSize)>3)
  
  # if (nrow(outlier_refregsize_erodedwm)>0) {
  #   message("**********************")
  #   message(paste("The cluster size for the Eroded WM was >3SD from the mean for",nrow(outlier_refregsize_erodedwm),"subject(s).","Check log file for details!"))
  #   message("**********************")
  #   
  #   for (i in 1:nrow(outlier_refregsize_erodedwm)) {
  #     tempid=outlier_refregsize_erodedwm[i,2]
  #     dec<-readline(prompt=paste("Should I keep",tempid,"in the database? Please answer y/n: "))
  #     outlier_refregsize_erodedwm[i,4]=dec
  #     
  #     if (dec=="n") {
  #       
  #       ftp_ext_erodedwm=subset(ftp_ext_erodedwm, !ftp_ext_erodedwm$ScanIdentifier==tempid)
  #       
  #     }
  #     
  #   }
  #   
  #   colnames(outlier_refregsize_erodedwm)[4]="Present_in_LONI_database"
  #   fname_erodedwmFTP=paste("logs/FTP_PETCore_ErodedWM_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  #   write.csv(outlier_refregsize_erodedwm,fname_erodedwmFTP, row.names = F, na="")
  #   
  # }
  
  
  ### we are theoretically ready to merge in the final ftpdf
  ftp_ext_erodedwm=ftp_ext_erodedwm[c(1:2,6:8)]
  colnames(ftp_ext_erodedwm)[c(3:5)]=c("ScalingFactor_ErodedWM","ErodedWM_MRIBASED_SUVR","ScalingFactor_ErodedWM_ClustSize")
  ftp_ext_erodedwm$ScalingFactor_ErodedWM=round(ftp_ext_erodedwm$ScalingFactor_ErodedWM, 4)
  ftpdf=merge(ftpdf, ftp_ext_erodedwm, by=c("ID","FTPPET_Date"))
  ftpdf=ftpdf[c(1:5,244,6:11,245,12:125,246,126:130,246,131:(ncol(ftpdf)-3))]
  colnames(ftpdf)[134]="ErodedWM_ClustSize"
  
  ### Done, ready to export!
  
  ##
  
  fname2=paste("ready_datasets/FTP_PETCore_Analysis_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(ftpdf,fname2, row.names = F, na="")

  #### Module for FTP-PET Quarterly report: To be done with report from MRICore as input as well ####

  ftpdf_qrtrep=merge(ftp_ext[1:3],ftpdf,by=c("ID","FTPPET_Date"))
  ftpdf_qrtrep$MRIident=paste(ftpdf_qrtrep$ID,"_",ftpdf_qrtrep$MRI_Date,sep="")

  ftpdf_qrtrep=subset(ftpdf_qrtrep, MRIident %in% mrireport$MRIident)
  ftpdf_qrtrep=ftpdf_qrtrep[c(1,2,4:248)]

  fname2=paste("ready_datasets/QuarterlyReport",format(Sys.time(), "%m%Y"),"_FTP_PETCore_Analysis_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(ftpdf_qrtrep,fname2, row.names = F, na="")
  
  #### Preparing FDG-PET database ####
  
  fdg_meta=fdg_meta[c(1,ncol(fdg_meta))]
  colnames(fdg_meta)[1]=c("ImageID")
  fdg_ext=fdg_ext[c(1:3,6:ncol(fdg_ext))]
  fdg_meta_ext=merge(fdg_meta, fdg_ext, by="ScanIdentifier")
  ptcohort=fbbdf[c("ID","CohortAssgn")]
  ptcohort=ptcohort[!duplicated(ptcohort$ID),]
  fdgdf=merge(ptcohort, fdg_meta_ext, by="ID",all.y = TRUE)
  fdgdf$CohortAssgn=ifelse(is.na(fdgdf$CohortAssgn),"NotFound_CheckEDC",as.character(fdgdf$CohortAssgn))
  fdgdf$CohortAssgn=factor(fdgdf$CohortAssgn, levels=c("EOAD","EOnonAD","CN","NotFound_CheckEDC"))
  
  #### Subset FDG-PET database based on MRI passing QC from Mayo ####
  
  fdgdf$MRI_ScanIdentifier=paste(fdgdf$ID,"_",fdgdf$MRI_Date,sep="")
  excluded_fdg_mriqc=subset(fdgdf, !(MRI_ScanIdentifier %in% mri_included$MRI_ScanIdentifier))
  fname1=paste("logs/FDG_PETCore_Excluded_MRIQC_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(excluded_fdg_mriqc,fname1, row.names = F, na="")
  fdgdf=subset(fdgdf, MRI_ScanIdentifier %in% mri_included$MRI_ScanIdentifier)
  
  #### Additional QC on FDG-PET database: Dynamic check on outliers in reference region size ####
  
  outlier_refregsize=fdgdf[c("ID","ScanIdentifier","ScalingFactor_Pons_ClustSize")]
  outlier_refregsize$ScalingFactor_Pons_ClustSize=scale(outlier_refregsize$ScalingFactor_Pons_ClustSize)
  outlier_refregsize=subset(outlier_refregsize, abs(outlier_refregsize$ScalingFactor_Pons_ClustSize)>3)
  
  # if (nrow(outlier_refregsize)>0) {
  #   message("**********************")
  #   message(paste("The cluster size for the Pons  was >3SD from the mean for",nrow(outlier_refregsize),"subject(s).","Check log file for details!"))
  #   message("**********************")
  #   for (i in 1:nrow(outlier_refregsize)) {
  #     tempid=outlier_refregsize[i,2]
  #     dec<-readline(prompt=paste("Should I keep",tempid,"in the database? Please answer y/n: "))
  #     outlier_refregsize[i,4]=dec
  #     
  #     if (dec=="n") {
  #       
  #       fdgdf=subset(fdgdf, !fdgdf$ScanIdentifier==tempid)
  #       
  #     }
  #     
  #   }
  #   
  #   colnames(outlier_refregsize)[4]="Present_in_LONI_database"
  #   fname3=paste("logs/FDG_PETCore_OutlierPons_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  #   write.csv(outlier_refregsize,fname3, row.names = F, na="")
  #   
  # }
  
  #### Sort final FDG-PET database and save it ####
  
  fdgdf=fdgdf[c(1,5,4,2,7:(ncol(fdgdf)-1))]
  fdgdf$ScalingFactor_Pons=round(fdgdf$ScalingFactor_Pons, 4)
  colnames(fdgdf)[6:118]=paste(names(fdgdf)[6:118],"_MRIBASED_SUVR",sep="")
  fname2=paste("ready_datasets/FDG_PETCore_Analysis_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(fdgdf,fname2, row.names = F, na="")
  
  #### Module for FDG-PET Quarterly report: To be done with report from MRICore as input as well ####
  
  fdgdf_qrtrep=merge(fdg_ext[1:3],fdgdf,by=c("ID","FDGPET_Date"))
  fdgdf_qrtrep$MRIident=paste(fdgdf_qrtrep$ID,"_",fdgdf_qrtrep$MRI_Date,sep="")

  fdgdf_qrtrep=subset(fdgdf_qrtrep, MRIident %in% mrireport$MRIident)
  fdgdf_qrtrep=fdgdf_qrtrep[c(1,2,4:233)]

  fname2=paste("ready_datasets/QuarterlyReport",format(Sys.time(), "%m%Y"),"_FDG_PETCore_Analysis_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(fdgdf_qrtrep,fname2, row.names = F, na="")

  #### Save workspace for later QC and Greetings ####

  rdatafname=paste("logs/Datasets_Generation_LONI_Step1_",format(Sys.time(), "%Y-%m-%d"),".Rdata",sep = "")
  save.image(file=rdatafname)
  
  message("**********************")
  message(paste("Finished! I saved",nrow(fbbdf),"rows for FBB,",nrow(ftpdf),"rows for FTP and",nrow(fdgdf),"for FDG. Enjoy!"))
  message("**********************")
  
  #### Report: Saving databases for MRI-based A+ CN or A+ EOnonAD ####

  amyposCNEOnonAD=subset(fbbdf, fbbdf$MRIBASED_Composite_SUVR_Type2>1.20 & (fbbdf$CohortAssgn=="CN" | fbbdf$CohortAssgn=="EOnonAD"))

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Database_CNEOnonADamypos1p20.xlsx",sep = "")

  if (nrow(amyposCNEOnonAD)>0) {
    write.xlsx(amyposCNEOnonAD, dbname, row.names=FALSE, col.names = TRUE)
  }

  amyposCNEOnonAD=subset(fbbdf, fbbdf$MRIBASED_Composite_SUVR_Type2>1.08 & (fbbdf$CohortAssgn=="CN" | fbbdf$CohortAssgn=="EOnonAD"))

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Database_CNEOnonADamypos1p08.xlsx",sep = "")

  if (nrow(amyposCNEOnonAD)>0) {
    write.xlsx(amyposCNEOnonAD, dbname, row.names=FALSE, col.names = TRUE)
  }

  #### Report: Saving database for MRI-based A- EOAD ####

  amynegEOAD=subset(fbbdf, (fbbdf$MRIBASED_Composite_SUVR_Type2<=1.20 & fbbdf$CohortAssgn=="EOAD"))

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Database_EOAD_ADamyneg.xlsx",sep = "")

  if (nrow(amynegEOAD)>0) {
    write.xlsx(amynegEOAD, dbname, row.names=FALSE, col.names = TRUE)
  }

  #### Report: Saving databases for discordant cases (visual vs. quantification) ####

  PTQneg1p20Vpos=subset(fbbdf, (fbbdf$MRIBASED_Composite_SUVR_Type2<=1.20 & fbbdf$Screening_PETONLY_VisualRead=="1"))

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Database_PT_Qneg1p20Vpos.xlsx",sep = "")

  if (nrow(PTQneg1p20Vpos)>0) {
    write.xlsx(PTQneg1p20Vpos, dbname, row.names=FALSE, col.names = TRUE)
  }

  PTQneg1p08Vpos=subset(fbbdf, (fbbdf$MRIBASED_Composite_SUVR_Type2<=1.08 & fbbdf$Screening_PETONLY_VisualRead=="1"))

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Database_PT_Qneg1p08Vpos.xlsx",sep = "")

  if (nrow(PTQneg1p08Vpos)>0) {
    write.xlsx(PTQneg1p08Vpos, dbname, row.names=FALSE, col.names = TRUE)
  }

  PTQpos1p20Vneg=subset(fbbdf, (fbbdf$MRIBASED_Composite_SUVR_Type2>1.20 & fbbdf$Screening_PETONLY_VisualRead=="0"))

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Database_PT_Qpos1p20Vneg.xlsx",sep = "")

  if (nrow(PTQpos1p20Vneg)>0) {
    write.xlsx(PTQpos1p20Vneg, dbname, row.names=FALSE, col.names = TRUE)
  }

  PTQpos1p08Vneg=subset(fbbdf, (fbbdf$MRIBASED_Composite_SUVR_Type2>1.08 & fbbdf$Screening_PETONLY_VisualRead=="0"))

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Database_PT_Qpos1p08Vneg.xlsx",sep = "")

  if (nrow(PTQpos1p08Vneg)>0) {
    write.xlsx(PTQpos1p08Vneg, dbname, row.names=FALSE, col.names = TRUE)
  }

  #### Report: Preparing baseline datasets ####

  fbbdf_bl=data.frame()
  allids=unique(fbbdf$ID)

  for (i in 1:length(unique(fbbdf$ID))) {

    test=allids[i]
    tempfbbdf=subset(fbbdf, ID %in% test)
    tempfbbdf$FBBPET_Date=as.Date.factor(tempfbbdf$FBBPET_Date)

    if (nrow(tempfbbdf)==1) {

      fbbdf_bl=rbind(fbbdf_bl, tempfbbdf)

    } else if (nrow(tempfbbdf)>1) {

      tempfbbdf=subset(tempfbbdf, tempfbbdf$FBBPET_Date==min(as.Date(tempfbbdf$FBBPET_Date)))
      fbbdf_bl=rbind(fbbdf_bl, tempfbbdf)
    }

  }

  ftpdf_bl=data.frame()
  allids=unique(ftpdf$ID)

  for (i in 1:length(unique(ftpdf$ID))) {

    test=allids[i]
    tempftpdf=subset(ftpdf, ID %in% test)
    tempftpdf$FTPPET_Date=as.Date.factor(tempftpdf$FTPPET_Date)

    if (nrow(tempftpdf)==1) {

      ftpdf_bl=rbind(ftpdf_bl, tempftpdf)

    } else if (nrow(tempftpdf)>1) {

      tempftpdf=subset(tempftpdf, tempftpdf$FTPPET_Date==min(as.Date(tempftpdf$FTPPET_Date)))
      ftpdf_bl=rbind(ftpdf_bl, tempftpdf)
    }

  }

  fbbftpdf_bl=merge(fbbdf_bl[c(1:14,20)], ftpdf_bl[1:13], by="ID")

  #### Report: ROC analysis for amyloid-positivity based on visual read in our patients ####

  #https://cran.r-project.org/web/packages/cutpointr/vignettes/cutpointr.html

  fbbdf_bl_pt=subset(fbbdf_bl, !fbbdf_bl$CohortAssgn=="CN")

  opt_cut <- cutpointr(fbbdf_bl_pt, MRIBASED_Composite_Centiloids, Screening_PETONLY_VisualRead, direction = ">=", pos_class = "1",
                       neg_class = "0", method = maximize_metric, metric = youden, use_midpoints = FALSE, na.rm=T)
  opt_cut_b <- cutpointr(fbbdf_bl_pt, MRIBASED_Composite_Centiloids, Screening_PETONLY_VisualRead, direction = ">=", pos_class = "1",
                         neg_class = "0", method = maximize_metric, metric = youden, boot_runs=1000, na.rm=T)

  opt_cut_p1=plot_metric(opt_cut_b)
  pname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Centiloids_by_VisualRead_ThresYouden.jpg",sep = "")
  ggsave(pname,plot=opt_cut_p1, height = 10, width=14, units = "cm")

  df_optcuts=t(as.data.frame(as.vector(summary(opt_cut_b$boot[[1]]$optimal_cutpoint))))
  df_optcuts=as.data.frame.character(df_optcuts)
  colnames(df_optcuts)="Youden Cutpoint"
  rownames(df_optcuts)=c("Min","1st_Quartile","Median","Mean","3rd_Quartile","Max")
  df_optcuts$`Youden Cutpoint`=as.numeric(as.character(df_optcuts$`Youden Cutpoint`))
  df_optcuts=rbind(df_optcuts,1000)
  rownames(df_optcuts)[nrow(df_optcuts)]=c("Nboot")

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Centiloids_by_VisualRead_BootsThresYouden.xlsx",sep = "")
  xlsx::write.xlsx(df_optcuts, dbname, row.names=TRUE, col.names = TRUE)

  #### Report: Plotting centiloids by visual read and cohort in patients ####

  ggdf_fbb_pt=fbbdf_bl[c(1,9,6,14)]
  ggdf_fbb_pt=subset(ggdf_fbb_pt, !ggdf_fbb_pt$CohortAssgn=="CN")
  ggdf_fbb_pt$Screening_PETONLY_VisualRead=factor(ggdf_fbb_pt$Screening_PETONLY_VisualRead, levels=c("1","0"))
  ggdfmelt_fbb_pt=melt(ggdf_fbb_pt, id.var=c("CohortAssgn","ID","Screening_PETONLY_VisualRead"))
  p=ggplot(data=ggdfmelt_fbb_pt, aes(x=variable, y=value, color=Screening_PETONLY_VisualRead)) +
    geom_point(aes(x=variable, y=value, fill=Screening_PETONLY_VisualRead), position=position_jitterdodge(), size=2.5, alpha=0.7) +
    ylab("18F-Florbetaben Centiloids") +
    xlab("18F-Florbetaben Visual Read") +
    theme_minimal() + labs(col="18F-Florbetaben \nVisual Read", fill="18F-Florbetaben \nVisual Read") +
    theme(axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold")) +
    # add mean lines + datapoints
    geom_point(aes(group=Screening_PETONLY_VisualRead),stat="summary", fun="median", size=4, position=position_dodge(0.75), color="white", shape=17) +
    # add mean lines + datapoints
    geom_point(aes(group=Screening_PETONLY_VisualRead),stat="summary", fun="median", size=3, position=position_dodge(0.75), color="black", shape=17) +
    scale_color_manual(labels=c("Amyloid Positive","Amyloid Negative"), values=c("#716FB2","#EB093C")) +
    scale_fill_manual(labels=c("Amyloid Positive","Amyloid Negative"), values=c("#716FB2","#EB093C")) +
    scale_x_discrete(labels=c("MRIBASED_Composite_Centiloids" = ""), expand=c(-0.5, 0)) + geom_hline(yintercept=opt_cut$optimal_cutpoint, lty="dashed", col="firebrick") +
    annotate(geom="text", x=1.375, y=(opt_cut$optimal_cutpoint-8), label=paste("Cutoff:",round(opt_cut$optimal_cutpoint,2),"\nAUC:",round(opt_cut$AUC,3), sep=""), color="firebrick", size=2)
  pname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Centiloids_by_PTCohort_WithThresDiscr.jpg",sep = "")
  ggsave(pname,plot=p, height = 14, width=12, units = "cm")

  # 2a. Same thing but with positivity threshold instead

  p=ggplot(data=ggdfmelt_fbb_pt, aes(x=variable, y=value, color=Screening_PETONLY_VisualRead)) +
    geom_point(aes(x=variable, y=value, fill=Screening_PETONLY_VisualRead), position=position_jitterdodge(), size=2.5, alpha=0.7) +
    ylab("18F-Florbetaben Centiloids") +
    xlab("18F-Florbetaben Visual Read") +
    theme_minimal() + labs(col="18F-Florbetaben \nVisual Read", fill="18F-Florbetaben \nVisual Read") +
    theme(axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold")) +
    # add mean lines + datapoints
    geom_point(aes(group=Screening_PETONLY_VisualRead),stat="summary", fun="median", size=4, position=position_dodge(0.75), color="white", shape=17) +
    # add mean lines + datapoints
    geom_point(aes(group=Screening_PETONLY_VisualRead),stat="summary", fun="median", size=3, position=position_dodge(0.75), color="black", shape=17) +
    scale_color_manual(labels=c("Amyloid Positive","Amyloid Negative"), values=c("#716FB2","#EB093C")) +
    scale_fill_manual(labels=c("Amyloid Positive","Amyloid Negative"), values=c("#716FB2","#EB093C")) +
    scale_x_discrete(labels=c("MRIBASED_Composite_Centiloids" = ""), expand=c(-0.5, 0)) + geom_hline(yintercept=36.7, lty="dashed", col="black") +
    annotate(geom="text", x=1.2, y=(36.7+8), label=paste("Pos. Threshold\nCL=36.7 | SUVR=1.20"),color="black", size=3)
  p_ctlsvsvisread_poscp=p
  pname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Centiloids_by_PTCohort_WithPosThres.jpg",sep = "")
  ggsave(pname,plot=p, height = 14, width=12, units = "cm")

  # 3. Plotting centiloids by cohort

  ggdf_fbb=fbbdf_bl[c(1,9,12,14)]
  ggdfmelt_fbb=melt(ggdf_fbb, id.var=c("CohortAssgn","ID"))
  p=ggplot(data=subset(ggdfmelt_fbb, ggdfmelt_fbb$variable=="MRIBASED_Composite_Centiloids"), aes(x=variable, y=value, color=CohortAssgn)) +
    geom_point(aes(x=variable, y=value, fill=CohortAssgn), position=position_jitterdodge(), size=2.5, alpha=0.7) +
    ylab("18F-Florbetaben Centiloids") +
    xlab("Cohort") +
    theme_minimal() + labs(col="Subgroup", fill="Subgroup") +
    theme(axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold")) +
    # add mean lines + datapoints
    geom_point(aes(group=CohortAssgn),stat="summary", fun="median", size=4, position=position_dodge(0.75), color="white", shape=17) +
    # add mean lines + datapoints
    geom_point(aes(group=CohortAssgn),stat="summary", fun="median", size=3, position=position_dodge(0.75), color="#F26D04", shape=17) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF")) + geom_hline(yintercept=36.7, lty="dashed", col="black") +
    scale_x_discrete(labels=c("MRIBASED_Composite_Centiloids" = ""), expand=c(-0.5, 0)) +
    annotate(geom="text", x=1.2, y=(36.7+8), label=paste("Pos. Threshold\nCL=36.7 | SUVR=1.20"),color="black", size=3)
  p_ctlsbycohort_poscp=p
  pname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Centiloids_by_Cohort.jpg",sep = "")
  ggsave(pname,plot=p, height = 14, width=12, units = "cm")

  # 4. Plotting centiloids by cohort + Pos threshold

  p=ggplot(data=subset(ggdfmelt_fbb, ggdfmelt_fbb$variable=="MRIBASED_Composite_Centiloids"), aes(x=variable, y=value, color=CohortAssgn)) +
    geom_point(aes(x=variable, y=value, fill=CohortAssgn), position=position_jitterdodge(), size=2.5, alpha=0.7) +
    ylab("18F-Florbetaben Centiloids") +
    xlab("Cohort") +
    theme_minimal() + labs(col="Subgroup", fill="Subgroup") +
    theme(axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold")) +
    # add mean lines + datapoints
    geom_point(aes(group=CohortAssgn),stat="summary", fun="median", size=4, position=position_dodge(0.75), color="white", shape=17) +
    # add mean lines + datapoints
    geom_point(aes(group=CohortAssgn),stat="summary", fun="median", size=3, position=position_dodge(0.75), color="#F26D04", shape=17) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF")) +
    scale_x_discrete(labels=c("MRIBASED_Composite_Centiloids" = ""), expand=c(-0.5, 0)) + geom_hline(yintercept=36.7, lty="dashed", col="firebrick") +
    annotate(geom="text", x=0.68, y=(36.7+4), label=paste("Quant Thres:36.7"), color="firebrick", size=3)

  pname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Centiloids_by_Cohort_WithThres.jpg",sep = "")
  ggsave(pname,plot=p, height = 14, width=12, units = "cm")

  #### Report: Reading available clinical data ####

  ### Module to add clinical data and create summaries
  ### Added July 1st 2020

  ### read MMSE

  mmsefiles=file.info(list.files("service/",pattern="LEADS_MMSE*", full.names = T))
  mmsedb=read.csv(rownames(mmsefiles)[order(mmsefiles$mtime)][nrow(mmsefiles)])
  mmsedb=mmsedb[c("subject_code","event_code","mmscore")] #isolate what's needed
  mmsedb=subset(mmsedb, mmsedb$event_code=="sc") #exclude follow-up scores
  mmsedb$subject_code=gsub("lds","LDS",mmsedb$subject_code) #IDs now conform to other databases
  mmsedb$mmscore=ifelse(mmsedb$mmscore>30,NA,
                        ifelse(mmsedb$mmscore<0,NA,mmsedb$mmscore)) # clean possible weird values

  ### read NACC UDS

  udsfiles=file.info(list.files("service/",pattern="NACC_UDS*", full.names = T))
  udsdb=read.csv(rownames(udsfiles)[order(udsfiles$mtime)][nrow(udsfiles)])
  udsdb=udsdb[c("LEADS_ID","LEADS_SCREENING_VISIT","CDRGLOB","CDRSUM","MOCATOTS")] #isolate what's needed
  udsdb=subset(udsdb, udsdb$LEADS_SCREENING_VISIT=="sc") #exclude follow-up scores
  udsdb$LEADS_ID=gsub("lds","LDS",udsdb$LEADS_ID) #IDs now conform to other databases
  udsdb$MOCATOTS=ifelse(udsdb$MOCATOTS>30,NA,
                        ifelse(udsdb$MOCATOTS<0,NA,udsdb$MOCATOTS)) # clean possible weird values

  ### read APOE

  apoefiles=file.info(list.files("service/",pattern="Biospecimen_Analysis_Results*", full.names = T))
  apoedb=read.csv(rownames(apoefiles)[order(apoefiles$mtime)][nrow(apoefiles)])
  apoedb=subset(apoedb, apoedb$TESTNAME=="APOE Genotype")
  apoedb=apoedb[c("SUBJECT_CODE","TESTVALUE")] #isolate what's needed
  apoedb$apoe_pos=ifelse(grepl("4",apoedb$TESTVALUE),"1","0")

  ### read Demographics

  demfiles=file.info(list.files("service/",pattern="LEADS_PTDEMOG*", full.names = T))
  demdb=read.csv(rownames(demfiles)[order(demfiles$mtime)][nrow(demfiles)])
  demdb=demdb[c("subject_code","ptgender","ptdob","pteducat")] #isolate what's needed
  demdb$subject_code=gsub("lds","LDS",demdb$subject_code) #IDs now conform to other databases
  demdb$ptdob=as.Date(demdb$ptdob)

  colnames(mmsedb)=c("ID","visit","MMSE")
  colnames(udsdb)=c("ID","visit","CDRtot","CDRsb","MOCA")
  colnames(apoedb)=c("ID","APOE genotype","APOEpos")
  colnames(demdb)=c("ID","sex","dob","edu")

  apoedb=apoedb[!duplicated(apoedb$ID),] # There were some duplicated in the APOE database

  #### Report: Merging baseline datasets with clinical and demographic data ####

  fbbdf_bl=merge(fbbdf_bl,demdb, by="ID", all.x=TRUE)
  fbbdf_bl=merge(fbbdf_bl,apoedb, by="ID", all.x=TRUE)
  fbbdf_bl=merge(fbbdf_bl,mmsedb, by="ID", all.x=TRUE)
  fbbdf_bl=merge(fbbdf_bl,udsdb, by="ID", all.x=TRUE)

  ftpdf_bl=merge(ftpdf_bl,demdb, by="ID", all.x=TRUE)
  ftpdf_bl=merge(ftpdf_bl,apoedb, by="ID", all.x=TRUE)
  ftpdf_bl=merge(ftpdf_bl,mmsedb, by="ID", all.x=TRUE)
  ftpdf_bl=merge(ftpdf_bl,udsdb, by="ID", all.x=TRUE)

  fbbftpdf_bl=merge(fbbftpdf_bl,demdb, by="ID", all.x=TRUE)
  fbbftpdf_bl=merge(fbbftpdf_bl,apoedb, by="ID", all.x=TRUE)
  fbbftpdf_bl=merge(fbbftpdf_bl,mmsedb, by="ID", all.x=TRUE)
  fbbftpdf_bl=merge(fbbftpdf_bl,udsdb, by="ID", all.x=TRUE)

  fbbdf_bl$sex=ifelse(fbbdf_bl$sex=="1","Male",
                      ifelse(is.na(fbbdf_bl$sex),NA,"Female"))

  ftpdf_bl$sex=ifelse(ftpdf_bl$sex=="1","Male",
                      ifelse(is.na(ftpdf_bl$sex),NA,"Female"))

  fbbftpdf_bl$sex=ifelse(fbbftpdf_bl$sex=="1","Male",
                         ifelse(is.na(fbbftpdf_bl$sex),NA,"Female"))

  fbbdf_bl$age=time_length(difftime(fbbdf_bl$FBBPET_Date, fbbdf_bl$dob), "years")
  ftpdf_bl$age=time_length(difftime(ftpdf_bl$FTPPET_Date, ftpdf_bl$dob), "years")
  fbbftpdf_bl$age=time_length(difftime(fbbftpdf_bl$FBBPET_Date, fbbftpdf_bl$dob), "years")

  #### Report: Generate summary tables ####
  
  mytab1=easytable(qualvar=c("sex","APOE genotype","APOEpos","CDRtot"),
            quantvar=c("age","edu","MMSE","MOCA","CDRsb","MRIBASED_Composite_SUVR_Type2","MRIBASED_Composite_Centiloids"),
            labqual=c("Sex","APOE genotype","APOE e4 carrier","CDR total"),
            labquant=c("Age","Education (years)","MMSE","MOCA","CDR sum of boxes","FBB-PET Neocortical SUVR","FBB-PET Centiloids"),
            ndec=4,
            idgroup=9,
            dataset=fbbdf_bl,
            pcorrmethod="fdr",
            method="p",
            label="reports/Summary_FBB_Baseline_",
            export=TRUE)

  mytab2=easytable(qualvar=c("sex","APOE genotype","APOEpos","CDRtot"),
            quantvar=c("age","edu","MMSE","MOCA","CDRsb","MetaROI_MRIBASED_SUVR","Braak_1_MRIBASED_SUVR","Braak_12_MRIBASED_SUVR","Braak_34_MRIBASED_SUVR","Braak_56_MRIBASED_SUVR"),
            labqual=c("Sex","APOE genotype","APOE e4 carrier","CDR total"),
            labquant=c("Age","Education (years)","MMSE","MOCA","CDR sum of boxes","FTP-PET metaROI SUVR","FTP-PET Braak I SUVR","FTP-PET Braak I/II SUVR","FTP-PET Braak III/IV SUVR","FTP-PET Braak V/VI SUVR"),
            ndec=4,
            idgroup=4,
            dataset=ftpdf_bl,
            pcorrmethod="fdr",
            method="p",
            label="reports/Summary_FTP_Baseline_",
            export=TRUE)

  mytab3=easytable(qualvar=c("sex","APOE genotype","APOEpos","CDRtot"),
            quantvar=c("age","edu","MMSE","MOCA","CDRsb","MRIBASED_Composite_SUVR_Type2","MRIBASED_Composite_Centiloids","MetaROI_MRIBASED_SUVR","Braak_1_MRIBASED_SUVR","Braak_12_MRIBASED_SUVR","Braak_34_MRIBASED_SUVR","Braak_56_MRIBASED_SUVR"),
            labqual=c("Sex","APOE genotype","APOE e4 carrier","CDR total"),
            labquant=c("Age","Education (years)","MMSE","MOCA","CDR sum of boxes","FBB-PET Neocortical SUVR","FBB-PET Centiloids","FTP-PET metaROI SUVR","FTP-PET Braak I SUVR","FTP-PET Braak I/II SUVR","FTP-PET Braak III/IV SUVR","FTP-PET Braak V/VI SUVR"),
            ndec=4,
            idgroup=9,
            dataset=fbbftpdf_bl,
            pcorrmethod="fdr",
            method="p",
            label="reports/Summary_FBBFTP_Baseline_",
            export=TRUE)

  #### Saving Baseline merged data for possible use ####

  fname=paste("ready_datasets/FBB_Baseline_ClinDem_Merged_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(fbbdf_bl,fname, row.names = F, na="")

  fname=paste("ready_datasets/FTP_Baseline_ClinDem_Merged_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(ftpdf_bl,fname, row.names = F, na="")

  fname=paste("ready_datasets/FBBFTP_Baseline_ClinDem_Merged_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(fbbftpdf_bl,fname, row.names = F, na="")

  #### Report: Adding the FDG-PET section ####
  
  fdgdf=merge(fdgdf,demdb, by="ID", all.x=TRUE)
  fdgdf=merge(fdgdf,apoedb, by="ID", all.x=TRUE)
  fdgdf=merge(fdgdf,mmsedb, by="ID", all.x=TRUE)
  fdgdf=merge(fdgdf,udsdb, by="ID", all.x=TRUE)
  
  fdgdf$sex=ifelse(fdgdf$sex=="1","Male",
                   ifelse(is.na(fdgdf$sex),NA,"Female"))
  
  fdgdf$age=time_length(difftime(fdgdf$FDGPET_Date, fdgdf$dob), "years")
  
  #### Report: Generate summary tables ####
  
  fdgdf=subset(fdgdf, (fdgdf$CohortAssgn=="EOnonAD" | fdgdf$CohortAssgn=="CN"))
  
  mytab4=easytable(qualvar=c("sex","APOE genotype","APOEpos","CDRtot"),
            quantvar=c("age","edu","MMSE","MOCA","CDRsb"),
            labqual=c("Sex","APOE genotype","APOE e4 carrier","CDR total"),
            labquant=c("Age","Education (years)","MMSE","MOCA","CDR sum of boxes"),
            ndec=4,
            idgroup=4,
            dataset=fdgdf,
            pcorrmethod="fdr",
            method="p",
            label="reports/Summary_FDG_Baseline_",
            export=TRUE)
  
  #### Saving Baseline merged data for possible use ####
  
  fname=paste("ready_datasets/FDG_Baseline_ClinDem_Merged_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(fdgdf,fname, row.names = F, na="")
  
  ### Saving another Rdata file including also merged baseline databases ###
  
  rdatafname=paste("logs/Datasets_Generation_LONI_Step2_",format(Sys.time(), "%Y-%m-%d"),".Rdata",sep = "")
  save.image(file=rdatafname)
  
  #### Report: QC MRI-based vs. PET-based global amyloid ####

  fbbdf_bl_sym=subset(fbbdf_bl, !fbbdf_bl$CohortAssgn=="CN")
  fbbdf_bl_sym$CohortAssgn=droplevels(fbbdf_bl_sym$CohortAssgn)
  uplim=max(c(max(fbbdf_bl_sym$MRIBASED_Composite_SUVR_Type2),max(fbbdf_bl_sym$Screening_PETONLY_Composite_SUVR)))+0.10
  lolim=min(c(min(fbbdf_bl_sym$MRIBASED_Composite_SUVR_Type2),min(fbbdf_bl_sym$Screening_PETONLY_Composite_SUVR)))-0.10

  p=ggplot(data=fbbdf_bl_sym, aes(x=MRIBASED_Composite_SUVR_Type2, y=Screening_PETONLY_Composite_SUVR, color=CohortAssgn)) +geom_abline(slope=1,intercept = 0, lty="dotted") +
    geom_point(aes(x=MRIBASED_Composite_SUVR_Type2, y=Screening_PETONLY_Composite_SUVR, fill=CohortAssgn), size=2.5, alpha=0.7) +
    ylab("Neocortical PET-only FBB-PET SUVR") +
    xlab("Neocortical MRI-based FBB-PET SUVR") +
    theme_minimal() + labs(col="CohortAssgn", fill="CohortAssgn")  +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
    geom_smooth(aes(x = MRIBASED_Composite_SUVR_Type2, y = Screening_PETONLY_Composite_SUVR), data=fbbdf_bl_sym, method = "loess", se = TRUE, cex=1, col="red", lty="dashed") +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF")) + scale_x_continuous(limits=c(lolim,uplim), breaks = (seq(0.75, 2.5, 0.25))) +
    scale_y_continuous(limits=c(lolim,uplim), breaks = (seq(0.75, 2.5, 0.25))) +geom_hline(yintercept=1.18,lty="dashed",col="purple") +geom_vline(xintercept=1.20, lty="dashed",col="purple")
  p_qcpetonlyvsmribasedamy=p
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBB_PETonlyvsMRI.jpg",sep = "")
  ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

  #### Report: ADNI-style braak ROIs plot (LIANA's request) ####

  adnistdf=fbbftpdf_bl[c("ID","CohortAssgn.x","MRIBASED_Composite_SUVR_Type2","MetaROI_MRIBASED_SUVR","Braak_1_MRIBASED_SUVR","Braak_34_MRIBASED_SUVR","Braak_56_MRIBASED_SUVR")]
  adnistdf$AmyStatus=ifelse(adnistdf$MRIBASED_Composite_SUVR_Type2>1.20,"A+","A-")
  adnistdf$AmyStatus=factor(adnistdf$AmyStatus, levels=c("A-","A+"))
  adnistdf$Cohort=ifelse(adnistdf$CohortAssgn.x=="CN","CN","CI")
  adnistdf$Cohort=factor(adnistdf$Cohort, levels=c("CN","CI"))
  ggadnistdf=melt(adnistdf, id.vars = c(1,2,3,8,9), variable.name = "ROI", value.name = "SUVR")
  ggadnistdf$ROI=gsub("_MRIBASED_SUVR","",ggadnistdf$ROI)
  ggadnistdf$ROI=factor(ggadnistdf$ROI, levels=c("MetaROI","Braak_1","Braak_34","Braak_56"))
  ggadnistdf$CohortAssgn.x=factor(ggadnistdf$CohortAssgn.x, levels=c("CN","EOnonAD","EOAD"))

  allrois=levels(ggadnistdf$ROI)

  for (i in 1:length(allrois)) {

    test=allrois[i]
    tempggadnistdf=subset(ggadnistdf, ggadnistdf$ROI %in% test)
    tempggadnistdf_cnaneg=subset(tempggadnistdf, (tempggadnistdf$Cohort=="CN" & tempggadnistdf$AmyStatus=="A-"))
    tempggadnistdf_ciapos=subset(tempggadnistdf, (tempggadnistdf$Cohort!="CN" & tempggadnistdf$AmyStatus=="A+"))
    up90th=quantile(tempggadnistdf_cnaneg$SUVR, prob = seq(0, 1, length = 21))[19]

    prop=paste(round(length(which(tempggadnistdf_ciapos$SUVR<up90th))/nrow(tempggadnistdf_ciapos)*100,1),"%",sep="")

  p=ggplot(data=tempggadnistdf) + geom_boxplot(aes(x=Cohort, y=SUVR, col=AmyStatus), outlier.shape = NA) +
    geom_point(aes(x=Cohort, y=SUVR, fill=AmyStatus), position=position_jitterdodge(jitter.width=0.4), pch=21, size=2.5, alpha=0.7, colour="black", stroke=1) + ylab("FTP-PET SUVR") + xlab("") +
    theme_minimal() + labs(title=test, col="Amyloid Status", fill="Amyloid Status") + scale_color_manual(values=c("black","black")) + scale_fill_manual(values=c("white","black")) +
    theme(legend.title=element_text(size=20), legend.text = element_text(size=20), plot.title = element_text(size=20), legend.position = "bottom", axis.title.x=element_text(size=20), axis.text.x = element_text(size=20, face="bold"), axis.text.y=element_text(size=20, face="bold"), axis.title.y=element_text(size=20)) +
    geom_hline(yintercept=up90th,lty="dashed",col="red", lwd=1) + scale_y_continuous(limits=c(0.5,3.5)) +
    annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 0.8, ymax = up90th,
             alpha = .2, fill="red") + annotate("text",x=2.2,y=0.6,label=prop,size=10,col="red",fontface=2)

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_ADNIStyle_FTPPETSUVR_",test,".jpg",sep = "")
  ggsave(dbname, plot=p, width = 12, height = 12, units = "cm")


  }

  #### Report: FTP-PET metaROI vs. global amyloid ####

  p=ggplot(data=fbbftpdf_bl, aes(x=MetaROI_MRIBASED_SUVR, y=MRIBASED_Composite_SUVR_Type2, color=CohortAssgn.x)) +
    geom_point(aes(x=MetaROI_MRIBASED_SUVR, y=MRIBASED_Composite_SUVR_Type2, fill=CohortAssgn.x), size=2.5, alpha=0.7) +
    ylab("Neocortical FBB-PET SUVR") +
    xlab("Temporal MetaROI FTP-PET SUVR") +
    theme_minimal() + labs(col="CohortAssgn", fill="CohortAssgn")  +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
    geom_smooth(aes(x = MetaROI_MRIBASED_SUVR, y = MRIBASED_Composite_SUVR_Type2), data=fbbftpdf_bl, method = "loess", se = TRUE, cex=1, col="red", lty="dashed") +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

  pb=ggplot(data=fbbftpdf_bl, aes(x=MetaROI_MRIBASED_SUVR, y=MRIBASED_Composite_SUVR_Type2, color=CohortAssgn.x)) +
    geom_point(aes(x=MetaROI_MRIBASED_SUVR, y=MRIBASED_Composite_SUVR_Type2, fill=CohortAssgn.x), size=2.5, alpha=0.7) +
    ylab("Neocortical FBB-PET SUVR") +
    xlab("Temporal MetaROI FTP-PET SUVR") +
    theme_minimal() + labs(col="CohortAssgn", fill="CohortAssgn")  +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
    geom_smooth(aes(x = MetaROI_MRIBASED_SUVR, y = MRIBASED_Composite_SUVR_Type2), data=fbbftpdf_bl, method = "loess", se = FALSE, cex=1, col="red", lty="dashed") + geom_ribbon(aes(x = MetaROI_MRIBASED_SUVR, y = MRIBASED_Composite_SUVR_Type2, group=1), data=fbbftpdf_bl, stat = "smooth", method = "loess", alpha = .15) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))
  p_metaroiamy_all=pb

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_GlobalAmyvsMetaFTP.jpg",sep = "")
  ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

  p=ggplot(data=subset(fbbftpdf_bl, !fbbftpdf_bl$CohortAssgn.x=="CN"), aes(x=MetaROI_MRIBASED_SUVR, y=Screening_PETONLY_Composite_SUVR, color=CohortAssgn.x)) +
    geom_point(aes(x=MetaROI_MRIBASED_SUVR, y=Screening_PETONLY_Composite_SUVR, fill=CohortAssgn.x), size=2.5, alpha=0.7) +
    ylab("Neocortical FBB-PET SUVR (PET-based)") +
    xlab("Temporal MetaROI FTP-PET SUVR") +
    theme_minimal() + labs(col="CohortAssgn", fill="CohortAssgn")  +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
    geom_smooth(aes(x = MetaROI_MRIBASED_SUVR, y = Screening_PETONLY_Composite_SUVR), data=subset(fbbftpdf_bl, !fbbftpdf_bl$CohortAssgn.x=="CN"), method = "loess", se = TRUE, cex=1, col="red", lty="dashed") +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))
  p_metaroiamy_nocn=p
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_GlobalAmyPETPTONLYvsMetaFTP.jpg",sep = "")
  ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

  p=ggplot(data=fbbftpdf_bl, aes(x=MetaROI_MRIBASED_SUVR, y=MRIBASED_Composite_Centiloids, color=CohortAssgn.x)) +
    geom_point(aes(x=MetaROI_MRIBASED_SUVR, y=MRIBASED_Composite_Centiloids, fill=CohortAssgn.x), size=2.5, alpha=0.7) +
    ylab("Neocortical FBB-PET Centiloids") +
    xlab("Temporal MetaROI FTP-PET SUVR") +
    theme_minimal() + labs(col="CohortAssgn", fill="CohortAssgn")  +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
    geom_smooth(aes(x = MetaROI_MRIBASED_SUVR, y = MRIBASED_Composite_Centiloids), data=fbbftpdf_bl, method = "loess", se = TRUE, cex=1, col="red", lty="dashed") +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

  pb=ggplot(data=fbbftpdf_bl, aes(x=MetaROI_MRIBASED_SUVR, y=MRIBASED_Composite_Centiloids, color=CohortAssgn.x)) +
    geom_point(aes(x=MetaROI_MRIBASED_SUVR, y=MRIBASED_Composite_Centiloids, fill=CohortAssgn.x), size=2.5, alpha=0.7) +
    ylab("Neocortical FBB-PET Centiloids") +
    xlab("Temporal MetaROI FTP-PET SUVR") +
    theme_minimal() + labs(col="CohortAssgn", fill="CohortAssgn")  +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
    geom_smooth(aes(x = MetaROI_MRIBASED_SUVR, y = MRIBASED_Composite_Centiloids), data=fbbftpdf_bl, method = "loess", se = FALSE, cex=1, col="red", lty="dashed") + geom_ribbon(aes(x = MetaROI_MRIBASED_SUVR, y = MRIBASED_Composite_Centiloids, group=1), data=fbbftpdf_bl, stat = "smooth", method = "loess", alpha = .15) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

  p_metaroiamyctls_all=pb
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_GlobalAmyCTLvsMetaFTP.jpg",sep = "")
  ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

  #### Report: FTP-PET SUVRs from metaROI and Braak regions ####

  ggdf=fbbftpdf_bl[c("ID","CohortAssgn.x","MetaROI_MRIBASED_SUVR","Braak_1_MRIBASED_SUVR","Braak_34_MRIBASED_SUVR","Braak_56_MRIBASED_SUVR")]
  ggdfmelt=melt(ggdf, id.var=c("ID","CohortAssgn.x"))

  p=ggplot(data=ggdfmelt, aes(x=variable, y=value, color=CohortAssgn.x)) +
    geom_point(aes(x=variable, y=value, fill=CohortAssgn.x), position=position_jitterdodge(), size=2.5, alpha=0.7) +
    ylab("18F-Flortaucipir SUVR") +
    xlab("") +
    theme_minimal() + labs(col="Subgroup", fill="Subgroup") +
    theme(axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold")) +
    # add mean lines + datapoints
    geom_point(aes(group=CohortAssgn.x),stat="summary", fun="median", size=4, position=position_dodge(0.75), color="white", shape=17) +
    # add mean lines + datapoints
    geom_point(aes(group=CohortAssgn.x),stat="summary", fun="median", size=3, position=position_dodge(0.75), color="#F26D04", shape=17) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF")) +
    scale_x_discrete(labels=c("MetaROI_MRIBASED_SUVR" = "AD metaROI", "Braak_1_MRIBASED_SUVR" = "ERC",
                              "Braak_34_MRIBASED_SUVR" = "Braak III/IV","Braak_56_MRIBASED_SUVR" = "Braak V/VI"), expand=c(-0.5, 0))

  p_ftppet_metabraak=p

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Dotplot_FTPROIs.jpg",sep = "")
  ggsave(dbname, plot=p, width = 15, height = 10, units = "cm")

  #### Report: FBB-PET Longitudinal data in macroROIs: handling and plotting - Composite WM version ####
  
  # recalculate SUVRs #
  
  fbbcwmextr=fbbdf[c(12:13,15:133)]
  
  fbbcwmextr_new=data.frame(matrix(nrow=nrow(fbbcwmextr)))
  
  for (c in 1:ncol(fbbcwmextr)) {
    fbbcwmextr_new=cbind(fbbcwmextr_new, fbbcwmextr[c]/fbbcwmextr[8])
  }
  
  fbbcwmextr_new=fbbcwmextr_new[c(2:ncol(fbbcwmextr_new))]
  
  # We have the recalculated SUVRs
  
  fbbdf_cwm=fbbdf
  
  # fill this dataframe with the new SUVRs
  
  fbbdf_cwm[c(12:13,15:133)]=fbbcwmextr_new
  
  # Now replace Centiloids with the new formula from Royse et al. 2021 #
  
  fbbdf_cwm$MRIBASED_Composite_Centiloids=(244.20*fbbdf_cwm$MRIBASED_Composite_SUVR_Type2)-170.80
  
  # we are ready to apply the same code to the "new" dataset
  
  fbbdf_lg_cwm=data.frame()
  allids=unique(fbbdf_cwm$ID)
  
  for (i in 1:length(unique(fbbdf_cwm$ID))) {
    
    test=allids[i]
    tempfbbdf_cwm=subset(fbbdf_cwm, ID %in% test)
    tempfbbdf_cwm$FBBPET_Date=as.Date.factor(tempfbbdf_cwm$FBBPET_Date)
    tempfbbdf_cwm <- tempfbbdf_cwm[order(tempfbbdf_cwm$FBBPET_Date),]
    
    if (nrow(tempfbbdf_cwm)>1) {
      
      ## Adding small module to compute longitudinal spaghetti plots on macroROIs
      
      wfro=vector()
      wpar=vector()
      wtem=vector()
      wocc=vector()
      wcin=vector()
      
      for (lll in 1:nrow(tempfbbdf_cwm)) {
        wfro[lll]=weighted.mean(tempfbbdf_cwm[lll,c('ctx_lh_caudalmiddlefrontal_MRIBASED_SUVR',
                                                    'ctx_lh_lateralorbitofrontal_MRIBASED_SUVR',
                                                    'ctx_lh_medialorbitofrontal_MRIBASED_SUVR',
                                                    'ctx_lh_rostralmiddlefrontal_MRIBASED_SUVR',
                                                    'ctx_lh_superiorfrontal_MRIBASED_SUVR',
                                                    'ctx_lh_frontalpole_MRIBASED_SUVR',
                                                    'ctx_rh_caudalmiddlefrontal_MRIBASED_SUVR',
                                                    'ctx_rh_lateralorbitofrontal_MRIBASED_SUVR',
                                                    'ctx_rh_medialorbitofrontal_MRIBASED_SUVR',
                                                    'ctx_rh_rostralmiddlefrontal_MRIBASED_SUVR',
                                                    'ctx_rh_superiorfrontal_MRIBASED_SUVR',
                                                    'ctx_rh_frontalpole_MRIBASED_SUVR',
                                                    'ctx_lh_parsopercularis_MRIBASED_SUVR',
                                                    'ctx_lh_parsorbitalis_MRIBASED_SUVR',
                                                    'ctx_lh_parstriangularis_MRIBASED_SUVR',
                                                    'ctx_rh_parsopercularis_MRIBASED_SUVR',
                                                    'ctx_rh_parsorbitalis_MRIBASED_SUVR',
                                                    'ctx_rh_parstriangularis_MRIBASED_SUVR',
                                                    'ctx_lh_paracentral_MRIBASED_SUVR',
                                                    'ctx_lh_precentral_MRIBASED_SUVR',
                                                    'ctx_rh_paracentral_MRIBASED_SUVR',
                                                    'ctx_rh_precentral_MRIBASED_SUVR')],tempfbbdf_cwm[lll,c('ctx_lh_caudalmiddlefrontal_ClustSize',
                                                                                                            'ctx_lh_lateralorbitofrontal_ClustSize',
                                                                                                            'ctx_lh_medialorbitofrontal_ClustSize',
                                                                                                            'ctx_lh_rostralmiddlefrontal_ClustSize',
                                                                                                            'ctx_lh_superiorfrontal_ClustSize',
                                                                                                            'ctx_lh_frontalpole_ClustSize',
                                                                                                            'ctx_rh_caudalmiddlefrontal_ClustSize',
                                                                                                            'ctx_rh_lateralorbitofrontal_ClustSize',
                                                                                                            'ctx_rh_medialorbitofrontal_ClustSize',
                                                                                                            'ctx_rh_rostralmiddlefrontal_ClustSize',
                                                                                                            'ctx_rh_superiorfrontal_ClustSize',
                                                                                                            'ctx_rh_frontalpole_ClustSize',
                                                                                                            'ctx_lh_parsopercularis_ClustSize',
                                                                                                            'ctx_lh_parsorbitalis_ClustSize',
                                                                                                            'ctx_lh_parstriangularis_ClustSize',
                                                                                                            'ctx_rh_parsopercularis_ClustSize',
                                                                                                            'ctx_rh_parsorbitalis_ClustSize',
                                                                                                            'ctx_rh_parstriangularis_ClustSize',
                                                                                                            'ctx_lh_paracentral_ClustSize',
                                                                                                            'ctx_lh_precentral_ClustSize',
                                                                                                            'ctx_rh_paracentral_ClustSize',
                                                                                                            'ctx_rh_precentral_ClustSize')])
        wpar[lll]=weighted.mean(tempfbbdf_cwm[lll,c('ctx_lh_inferiorparietal_MRIBASED_SUVR',
                                                    'ctx_lh_superiorparietal_MRIBASED_SUVR',
                                                    'ctx_rh_inferiorparietal_MRIBASED_SUVR',
                                                    'ctx_rh_superiorparietal_MRIBASED_SUVR',
                                                    'ctx_lh_supramarginal_MRIBASED_SUVR',
                                                    'ctx_rh_supramarginal_MRIBASED_SUVR',
                                                    'ctx_lh_postcentral_MRIBASED_SUVR',
                                                    'ctx_rh_postcentral_MRIBASED_SUVR',
                                                    'ctx_lh_precuneus_MRIBASED_SUVR',
                                                    'ctx_rh_precuneus_MRIBASED_SUVR')],tempfbbdf_cwm[lll,c('ctx_lh_inferiorparietal_ClustSize',
                                                                                                           'ctx_lh_superiorparietal_ClustSize',
                                                                                                           'ctx_rh_inferiorparietal_ClustSize',
                                                                                                           'ctx_rh_superiorparietal_ClustSize',
                                                                                                           'ctx_lh_supramarginal_ClustSize',
                                                                                                           'ctx_rh_supramarginal_ClustSize',
                                                                                                           'ctx_lh_postcentral_ClustSize',
                                                                                                           'ctx_rh_postcentral_ClustSize',
                                                                                                           'ctx_lh_precuneus_ClustSize',
                                                                                                           'ctx_rh_precuneus_ClustSize')])
        wtem[lll]=weighted.mean(tempfbbdf_cwm[lll,c('ctx_lh_inferiortemporal_MRIBASED_SUVR',
                                                    'ctx_lh_middletemporal_MRIBASED_SUVR',
                                                    'ctx_lh_superiortemporal_MRIBASED_SUVR',
                                                    'ctx_lh_temporalpole_MRIBASED_SUVR',
                                                    'ctx_lh_transversetemporal_MRIBASED_SUVR',
                                                    'ctx_rh_inferiortemporal_MRIBASED_SUVR',
                                                    'ctx_rh_middletemporal_MRIBASED_SUVR',
                                                    'ctx_rh_superiortemporal_MRIBASED_SUVR',
                                                    'ctx_rh_temporalpole_MRIBASED_SUVR',
                                                    'ctx_rh_transversetemporal_MRIBASED_SUVR',
                                                    'ctx_lh_bankssts_MRIBASED_SUVR',
                                                    'ctx_rh_bankssts_MRIBASED_SUVR',
                                                    'ctx_lh_entorhinal_MRIBASED_SUVR',
                                                    'ctx_rh_entorhinal_MRIBASED_SUVR',
                                                    'ctx_lh_parahippocampal_MRIBASED_SUVR',
                                                    'ctx_rh_parahippocampal_MRIBASED_SUVR',
                                                    'ctx_lh_fusiform_MRIBASED_SUVR',
                                                    'ctx_rh_fusiform_MRIBASED_SUVR')],tempfbbdf_cwm[lll,c('ctx_lh_inferiortemporal_ClustSize',
                                                                                                          'ctx_lh_middletemporal_ClustSize',
                                                                                                          'ctx_lh_superiortemporal_ClustSize',
                                                                                                          'ctx_lh_temporalpole_ClustSize',
                                                                                                          'ctx_lh_transversetemporal_ClustSize',
                                                                                                          'ctx_rh_inferiortemporal_ClustSize',
                                                                                                          'ctx_rh_middletemporal_ClustSize',
                                                                                                          'ctx_rh_superiortemporal_ClustSize',
                                                                                                          'ctx_rh_temporalpole_ClustSize',
                                                                                                          'ctx_rh_transversetemporal_ClustSize',
                                                                                                          'ctx_lh_bankssts_ClustSize',
                                                                                                          'ctx_rh_bankssts_ClustSize',
                                                                                                          'ctx_lh_entorhinal_ClustSize',
                                                                                                          'ctx_rh_entorhinal_ClustSize',
                                                                                                          'ctx_lh_parahippocampal_ClustSize',
                                                                                                          'ctx_rh_parahippocampal_ClustSize',
                                                                                                          'ctx_lh_fusiform_ClustSize',
                                                                                                          'ctx_rh_fusiform_ClustSize')])
        wocc[lll]=weighted.mean(tempfbbdf_cwm[lll,c('ctx_lh_lateraloccipital_MRIBASED_SUVR',
                                                    'ctx_rh_lateraloccipital_MRIBASED_SUVR',
                                                    'ctx_lh_lingual_MRIBASED_SUVR',
                                                    'ctx_rh_lingual_MRIBASED_SUVR',
                                                    'ctx_lh_cuneus_MRIBASED_SUVR',
                                                    'ctx_rh_cuneus_MRIBASED_SUVR',
                                                    'ctx_lh_pericalcarine_MRIBASED_SUVR',
                                                    'ctx_rh_pericalcarine_MRIBASED_SUVR')],tempfbbdf_cwm[lll,c('ctx_lh_lateraloccipital_ClustSize',
                                                                                                               'ctx_rh_lateraloccipital_ClustSize',
                                                                                                               'ctx_lh_lingual_ClustSize',
                                                                                                               'ctx_rh_lingual_ClustSize',
                                                                                                               'ctx_lh_cuneus_ClustSize',
                                                                                                               'ctx_rh_cuneus_ClustSize',
                                                                                                               'ctx_lh_pericalcarine_ClustSize',
                                                                                                               'ctx_rh_pericalcarine_ClustSize')])
        wcin[lll]=weighted.mean(tempfbbdf_cwm[lll,c('ctx_lh_caudalanteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_lh_isthmuscingulate_MRIBASED_SUVR',
                                                    'ctx_lh_posteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_lh_rostralanteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_rh_caudalanteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_rh_isthmuscingulate_MRIBASED_SUVR',
                                                    'ctx_rh_posteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_rh_rostralanteriorcingulate_MRIBASED_SUVR')],tempfbbdf_cwm[lll,c('ctx_lh_caudalanteriorcingulate_ClustSize',
                                                                                                                          'ctx_lh_isthmuscingulate_ClustSize',
                                                                                                                          'ctx_lh_posteriorcingulate_ClustSize',
                                                                                                                          'ctx_lh_rostralanteriorcingulate_ClustSize',
                                                                                                                          'ctx_rh_caudalanteriorcingulate_ClustSize',
                                                                                                                          'ctx_rh_isthmuscingulate_ClustSize',
                                                                                                                          'ctx_rh_posteriorcingulate_ClustSize',
                                                                                                                          'ctx_rh_rostralanteriorcingulate_ClustSize')])
        
        
      }
      
      tempfbbdf_cwm$Frontal=wfro
      tempfbbdf_cwm$Parietal=wpar
      tempfbbdf_cwm$Temporal=wtem
      tempfbbdf_cwm$Occipital=wocc
      tempfbbdf_cwm$Cingulate=wcin
      
      ### END of MODULE, now we can save also those macroROIs in the next step
      
      tempfbbdf_cwm=tempfbbdf_cwm[c("ID","CohortAssgn","FBBPET_Date","ScalingFactor_CompositeWM","MRIBASED_Composite_SUVR_Type2","MRIBASED_Composite_Centiloids","Frontal","Parietal","Temporal","Occipital","Cingulate")]
      tempfbbdf_cwm$Timepoint=seq(0,(nrow(tempfbbdf_cwm)-1),1)
      fbbdf_lg_cwm=rbind(fbbdf_lg_cwm, tempfbbdf_cwm)
      
    }
    
  }
  
  fbbdf_lg_cwm$Timepoint=factor(fbbdf_lg_cwm$Timepoint, levels=unique(fbbdf_lg_cwm$Timepoint))
  
  p=ggplot(fbbdf_lg_cwm, aes(x=Timepoint, y=ScalingFactor_CompositeWM, group=ID, colour = CohortAssgn)) +
    geom_line(alpha=0.3) + stat_summary(aes(group = CohortAssgn), geom = "line", fun = "median") + stat_summary(aes(group = CohortAssgn), geom = "point", fun = "median", shape = 17, size = 3)+
    ylab("Composite WM Scaling Factor") +xlab("Timepoint") + theme_minimal() +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF")) +
    theme(panel.grid.major.x = element_blank() ,panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + scale_x_discrete(expand=c(0.1,0.1))
  
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongFBB_CompWM_ScalingFactor.jpg",sep = "")
  ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")
  
  #### Report: Save longitudinal FBB-PET database for convenience ####
  
  alltps=unique(fbbdf_lg_cwm$Timepoint)
  dftps=list()
  for (ntp in 1:length(alltps)) {
    temptp=alltps[ntp]
    temptpdf=subset(fbbdf_lg_cwm, fbbdf_lg_cwm$Timepoint %in% temptp)  ## create and add new data frame
    colnames(temptpdf)[3:ncol(temptpdf)]=paste("TP",temptp,"_",colnames(temptpdf)[3:ncol(temptpdf)], sep="")
    dftps[[ntp]] <- temptpdf
  }
  fbbdf_lg_cwm_exp=Reduce(function(x,y)merge(x,y,by=c("ID","CohortAssgn"),all.x=TRUE), dftps)
  
  ### quick module to remove unwanted subjects for now ###
  
  # delids=c("LDS0110021","LDS0370163","LDS0730044","LDS0730055","LDS0730150","LDS1770107","LDS1770112","LDS1770181","LDS3600187")
  # fbbdf_lg_cwm=subset(fbbdf_lg_cwm, !fbbdf_lg_cwm$ID %in% delids)
  
  ###
  
  dbname=paste("ready_datasets/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBB_LongitudinalDatabaseCompWM.xlsx",sep = "")
  write.xlsx(fbbdf_lg_cwm, dbname, row.names=FALSE, col.names = TRUE)
  
  ### Added module to get other longitudinal plot type ###
  
  iddobfbb=fbbdf_bl[c(1,256,259,263)]
  fbbdf_lg_cwm_plot2=merge(fbbdf_lg_cwm,iddobfbb,by=c("ID"))
  fbbdf_lg_cwm_plot2$age=time_length(difftime(fbbdf_lg_cwm_plot2$FBBPET_Date, fbbdf_lg_cwm_plot2$dob), "years")
  fbbdf_lg_cwm_plot2$CDRtot=ifelse(is.na(fbbdf_lg_cwm_plot2$CDRtot),"NaN",as.character(fbbdf_lg_cwm_plot2$CDRtot))
  fbbdf_lg_cwm_plot2$CDRtot=factor(fbbdf_lg_cwm_plot2$CDRtot, levels=c("0.5","1","NaN"))
  
  fbblg_plot2=ggplot(subset(fbbdf_lg_cwm_plot2, fbbdf_lg_cwm_plot2$CohortAssgn=="EOAD"), aes(x=age, y=MRIBASED_Composite_Centiloids, group=ID, colour = CohortAssgn)) +
    geom_line(alpha=1, size=1) +
    ylab("18F-Florbetaben-PET Centiloids (Ref: CompWM)") +xlab("Age") + theme_minimal() +
    theme(axis.title.x=element_text(size=16), axis.text.x = element_text(size=14), axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=14)) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF")) + scale_y_continuous(limits=c(0,200),breaks=seq(0,200,25)) + labs(title="Longitudinal FBB-PET Centiloids") +
    theme(panel.grid.minor = element_blank(), legend.position = "none")
  
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongFBBCompWM_Plot2.jpg",sep = "")
  ggsave(dbname, plot=fbblg_plot2, width = 18, height = 15, units = "cm")
  
  fbblg_plot2_apoe=ggplot(subset(fbbdf_lg_cwm_plot2, fbbdf_lg_cwm_plot2$CohortAssgn=="EOAD"), aes(x=age, y=MRIBASED_Composite_Centiloids, group=ID, colour = APOEpos)) +
    geom_line(alpha=1, size=1) +
    ylab("18F-Florbetaben-PET Centiloids (Ref: CompWM)") +xlab("Age") + theme_minimal() + scale_color_wsj() +
    theme(axis.title.x=element_text(size=16), axis.text.x = element_text(size=14), axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=14)) + scale_y_continuous(limits=c(0,200),breaks=seq(0,200,25)) + labs(title="Longitudinal FBB-PET Centiloids") +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongFBBCompWM_Plot2_APOE.jpg",sep = "")
  ggsave(dbname, plot=fbblg_plot2_apoe, width = 18, height = 15, units = "cm")
  
  
  fbblg_plot2_cdr=ggplot(subset(fbbdf_lg_cwm_plot2, fbbdf_lg_cwm_plot2$CohortAssgn=="EOAD"), aes(x=age, y=MRIBASED_Composite_Centiloids, group=ID, colour = CDRtot)) +
    geom_line(alpha=1, size=1) +
    ylab("18F-Florbetaben-PET Centiloids (Ref: CompWM)") +xlab("Age") + theme_minimal() + scale_color_manual(values=c("#716FB2","#EB093C","#B4B9Bf")) +
    theme(axis.title.x=element_text(size=16), axis.text.x = element_text(size=14), axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=14)) + scale_y_continuous(limits=c(0,200),breaks=seq(0,200,25)) + labs(title="Longitudinal FBB-PET Centiloids") +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongFBBCompWM_Plot2_CDR.jpg",sep = "")
  ggsave(dbname, plot=fbblg_plot2_cdr, width = 18, height = 15, units = "cm")
  
  
  #### Report: FTP-PET Longitudinal data in macroROIs: handling and plotting - Eroded WM version ####
  
  # recalculate SUVRs #
  
  ftpewmextr=ftpdf[c(8:126)]
  
  ftpewmextr_new=data.frame(matrix(nrow=nrow(ftpewmextr)))
  
  for (c in 1:ncol(ftpewmextr)) {
    ftpewmextr_new=cbind(ftpewmextr_new, ftpewmextr[c]/ftpewmextr[6])
  }
  
  ftpewmextr_new=ftpewmextr_new[c(2:ncol(ftpewmextr_new))]
  
  # We have the recalculated SUVRs
  
  ftpdf_ewm=ftpdf
  
  # fill this dataframe with the new SUVRs
  
  ftpdf_ewm[c(8:126)]=ftpewmextr_new
  
  ftpdf_ewm_lg=data.frame()
  allids=unique(ftpdf_ewm$ID)
  
  for (i in 1:length(unique(ftpdf_ewm$ID))) {
    
    test=allids[i]
    tempftpdf_ewm=subset(ftpdf_ewm, ID %in% test)
    tempftpdf_ewm$FTPPET_Date=as.Date.factor(tempftpdf_ewm$FTPPET_Date)
    tempftpdf_ewm <- tempftpdf_ewm[order(tempftpdf_ewm$FTPPET_Date),]
    
    if (nrow(tempftpdf_ewm)>1) {
      
      ## Adding small module to compute longitudinal spaghetti plots on macroROIs
      
      wfro=vector()
      wpar=vector()
      wtem=vector()
      wocc=vector()
      wcin=vector()
      
      for (lll in 1:nrow(tempftpdf_ewm)) {
        wfro[lll]=weighted.mean(tempftpdf_ewm[lll,c('ctx_lh_caudalmiddlefrontal_MRIBASED_SUVR',
                                                    'ctx_lh_lateralorbitofrontal_MRIBASED_SUVR',
                                                    'ctx_lh_medialorbitofrontal_MRIBASED_SUVR',
                                                    'ctx_lh_rostralmiddlefrontal_MRIBASED_SUVR',
                                                    'ctx_lh_superiorfrontal_MRIBASED_SUVR',
                                                    'ctx_lh_frontalpole_MRIBASED_SUVR',
                                                    'ctx_rh_caudalmiddlefrontal_MRIBASED_SUVR',
                                                    'ctx_rh_lateralorbitofrontal_MRIBASED_SUVR',
                                                    'ctx_rh_medialorbitofrontal_MRIBASED_SUVR',
                                                    'ctx_rh_rostralmiddlefrontal_MRIBASED_SUVR',
                                                    'ctx_rh_superiorfrontal_MRIBASED_SUVR',
                                                    'ctx_rh_frontalpole_MRIBASED_SUVR',
                                                    'ctx_lh_parsopercularis_MRIBASED_SUVR',
                                                    'ctx_lh_parsorbitalis_MRIBASED_SUVR',
                                                    'ctx_lh_parstriangularis_MRIBASED_SUVR',
                                                    'ctx_rh_parsopercularis_MRIBASED_SUVR',
                                                    'ctx_rh_parsorbitalis_MRIBASED_SUVR',
                                                    'ctx_rh_parstriangularis_MRIBASED_SUVR',
                                                    'ctx_lh_paracentral_MRIBASED_SUVR',
                                                    'ctx_lh_precentral_MRIBASED_SUVR',
                                                    'ctx_rh_paracentral_MRIBASED_SUVR',
                                                    'ctx_rh_precentral_MRIBASED_SUVR')],tempftpdf_ewm[lll,c('ctx_lh_caudalmiddlefrontal_ClustSize',
                                                                                                            'ctx_lh_lateralorbitofrontal_ClustSize',
                                                                                                            'ctx_lh_medialorbitofrontal_ClustSize',
                                                                                                            'ctx_lh_rostralmiddlefrontal_ClustSize',
                                                                                                            'ctx_lh_superiorfrontal_ClustSize',
                                                                                                            'ctx_lh_frontalpole_ClustSize',
                                                                                                            'ctx_rh_caudalmiddlefrontal_ClustSize',
                                                                                                            'ctx_rh_lateralorbitofrontal_ClustSize',
                                                                                                            'ctx_rh_medialorbitofrontal_ClustSize',
                                                                                                            'ctx_rh_rostralmiddlefrontal_ClustSize',
                                                                                                            'ctx_rh_superiorfrontal_ClustSize',
                                                                                                            'ctx_rh_frontalpole_ClustSize',
                                                                                                            'ctx_lh_parsopercularis_ClustSize',
                                                                                                            'ctx_lh_parsorbitalis_ClustSize',
                                                                                                            'ctx_lh_parstriangularis_ClustSize',
                                                                                                            'ctx_rh_parsopercularis_ClustSize',
                                                                                                            'ctx_rh_parsorbitalis_ClustSize',
                                                                                                            'ctx_rh_parstriangularis_ClustSize',
                                                                                                            'ctx_lh_paracentral_ClustSize',
                                                                                                            'ctx_lh_precentral_ClustSize',
                                                                                                            'ctx_rh_paracentral_ClustSize',
                                                                                                            'ctx_rh_precentral_ClustSize')])
        wpar[lll]=weighted.mean(tempftpdf_ewm[lll,c('ctx_lh_inferiorparietal_MRIBASED_SUVR',
                                                    'ctx_lh_superiorparietal_MRIBASED_SUVR',
                                                    'ctx_rh_inferiorparietal_MRIBASED_SUVR',
                                                    'ctx_rh_superiorparietal_MRIBASED_SUVR',
                                                    'ctx_lh_supramarginal_MRIBASED_SUVR',
                                                    'ctx_rh_supramarginal_MRIBASED_SUVR',
                                                    'ctx_lh_postcentral_MRIBASED_SUVR',
                                                    'ctx_rh_postcentral_MRIBASED_SUVR',
                                                    'ctx_lh_precuneus_MRIBASED_SUVR',
                                                    'ctx_rh_precuneus_MRIBASED_SUVR')],tempftpdf_ewm[lll,c('ctx_lh_inferiorparietal_ClustSize',
                                                                                                           'ctx_lh_superiorparietal_ClustSize',
                                                                                                           'ctx_rh_inferiorparietal_ClustSize',
                                                                                                           'ctx_rh_superiorparietal_ClustSize',
                                                                                                           'ctx_lh_supramarginal_ClustSize',
                                                                                                           'ctx_rh_supramarginal_ClustSize',
                                                                                                           'ctx_lh_postcentral_ClustSize',
                                                                                                           'ctx_rh_postcentral_ClustSize',
                                                                                                           'ctx_lh_precuneus_ClustSize',
                                                                                                           'ctx_rh_precuneus_ClustSize')])
        wtem[lll]=weighted.mean(tempftpdf_ewm[lll,c('ctx_lh_inferiortemporal_MRIBASED_SUVR',
                                                    'ctx_lh_middletemporal_MRIBASED_SUVR',
                                                    'ctx_lh_superiortemporal_MRIBASED_SUVR',
                                                    'ctx_lh_temporalpole_MRIBASED_SUVR',
                                                    'ctx_lh_transversetemporal_MRIBASED_SUVR',
                                                    'ctx_rh_inferiortemporal_MRIBASED_SUVR',
                                                    'ctx_rh_middletemporal_MRIBASED_SUVR',
                                                    'ctx_rh_superiortemporal_MRIBASED_SUVR',
                                                    'ctx_rh_temporalpole_MRIBASED_SUVR',
                                                    'ctx_rh_transversetemporal_MRIBASED_SUVR',
                                                    'ctx_lh_bankssts_MRIBASED_SUVR',
                                                    'ctx_rh_bankssts_MRIBASED_SUVR',
                                                    'ctx_lh_entorhinal_MRIBASED_SUVR',
                                                    'ctx_rh_entorhinal_MRIBASED_SUVR',
                                                    'ctx_lh_parahippocampal_MRIBASED_SUVR',
                                                    'ctx_rh_parahippocampal_MRIBASED_SUVR',
                                                    'ctx_lh_fusiform_MRIBASED_SUVR',
                                                    'ctx_rh_fusiform_MRIBASED_SUVR')],tempftpdf_ewm[lll,c('ctx_lh_inferiortemporal_ClustSize',
                                                                                                          'ctx_lh_middletemporal_ClustSize',
                                                                                                          'ctx_lh_superiortemporal_ClustSize',
                                                                                                          'ctx_lh_temporalpole_ClustSize',
                                                                                                          'ctx_lh_transversetemporal_ClustSize',
                                                                                                          'ctx_rh_inferiortemporal_ClustSize',
                                                                                                          'ctx_rh_middletemporal_ClustSize',
                                                                                                          'ctx_rh_superiortemporal_ClustSize',
                                                                                                          'ctx_rh_temporalpole_ClustSize',
                                                                                                          'ctx_rh_transversetemporal_ClustSize',
                                                                                                          'ctx_lh_bankssts_ClustSize',
                                                                                                          'ctx_rh_bankssts_ClustSize',
                                                                                                          'ctx_lh_entorhinal_ClustSize',
                                                                                                          'ctx_rh_entorhinal_ClustSize',
                                                                                                          'ctx_lh_parahippocampal_ClustSize',
                                                                                                          'ctx_rh_parahippocampal_ClustSize',
                                                                                                          'ctx_lh_fusiform_ClustSize',
                                                                                                          'ctx_rh_fusiform_ClustSize')])
        wocc[lll]=weighted.mean(tempftpdf_ewm[lll,c('ctx_lh_lateraloccipital_MRIBASED_SUVR',
                                                    'ctx_rh_lateraloccipital_MRIBASED_SUVR',
                                                    'ctx_lh_lingual_MRIBASED_SUVR',
                                                    'ctx_rh_lingual_MRIBASED_SUVR',
                                                    'ctx_lh_cuneus_MRIBASED_SUVR',
                                                    'ctx_rh_cuneus_MRIBASED_SUVR',
                                                    'ctx_lh_pericalcarine_MRIBASED_SUVR',
                                                    'ctx_rh_pericalcarine_MRIBASED_SUVR')],tempftpdf_ewm[lll,c('ctx_lh_lateraloccipital_ClustSize',
                                                                                                               'ctx_rh_lateraloccipital_ClustSize',
                                                                                                               'ctx_lh_lingual_ClustSize',
                                                                                                               'ctx_rh_lingual_ClustSize',
                                                                                                               'ctx_lh_cuneus_ClustSize',
                                                                                                               'ctx_rh_cuneus_ClustSize',
                                                                                                               'ctx_lh_pericalcarine_ClustSize',
                                                                                                               'ctx_rh_pericalcarine_ClustSize')])
        wcin[lll]=weighted.mean(tempftpdf_ewm[lll,c('ctx_lh_caudalanteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_lh_isthmuscingulate_MRIBASED_SUVR',
                                                    'ctx_lh_posteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_lh_rostralanteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_rh_caudalanteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_rh_isthmuscingulate_MRIBASED_SUVR',
                                                    'ctx_rh_posteriorcingulate_MRIBASED_SUVR',
                                                    'ctx_rh_rostralanteriorcingulate_MRIBASED_SUVR')],tempftpdf_ewm[lll,c('ctx_lh_caudalanteriorcingulate_ClustSize',
                                                                                                                          'ctx_lh_isthmuscingulate_ClustSize',
                                                                                                                          'ctx_lh_posteriorcingulate_ClustSize',
                                                                                                                          'ctx_lh_rostralanteriorcingulate_ClustSize',
                                                                                                                          'ctx_rh_caudalanteriorcingulate_ClustSize',
                                                                                                                          'ctx_rh_isthmuscingulate_ClustSize',
                                                                                                                          'ctx_rh_posteriorcingulate_ClustSize',
                                                                                                                          'ctx_rh_rostralanteriorcingulate_ClustSize')])
        
        
      }
      
      tempftpdf_ewm$Frontal=wfro
      tempftpdf_ewm$Parietal=wpar
      tempftpdf_ewm$Temporal=wtem
      tempftpdf_ewm$Occipital=wocc
      tempftpdf_ewm$Cingulate=wcin
      
      ### END of MODULE, now we can save also those macroROIs in the next step
      
      tempftpdf_ewm=tempftpdf_ewm[c("ID","CohortAssgn","FTPPET_Date","ScalingFactor_ErodedWM","MetaROI_MRIBASED_SUVR","Frontal","Parietal","Temporal","Occipital","Cingulate")]
      tempftpdf_ewm$Timepoint=seq(0,(nrow(tempftpdf_ewm)-1),1)
      ftpdf_ewm_lg=rbind(ftpdf_ewm_lg, tempftpdf_ewm)
      
    }
    
  }
  
  ftpdf_ewm_lg$Timepoint=factor(ftpdf_ewm_lg$Timepoint, levels=unique(ftpdf_ewm_lg$Timepoint))
  
  p=ggplot(ftpdf_ewm_lg, aes(x=Timepoint, y=ScalingFactor_ErodedWM, group=ID, colour = CohortAssgn)) +
    geom_line(alpha=0.3) + stat_summary(aes(group = CohortAssgn), geom = "line", fun= "median") + stat_summary(aes(group = CohortAssgn), geom = "point", fun = "median", shape = 17, size = 3)+
    ylab("Inf Cerebellar GM Scaling Factor") +xlab("Timepoint") + theme_minimal() +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF")) +
    theme(panel.grid.major.x = element_blank() ,panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + scale_x_discrete(expand=c(0.1,0.1))
  
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongFTP_EWM_ErodedWM_ScalingFactor.jpg",sep = "")
  ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")
  
  
  #### Report: Save longitudinal FTP-PET database for convenience ####
  
  alltps=unique(ftpdf_ewm_lg$Timepoint)
  dftps=list()
  for (ntp in 1:length(alltps)) {
    temptp=alltps[ntp]
    temptpdf=subset(ftpdf_ewm_lg, ftpdf_ewm_lg$Timepoint %in% temptp)  ## create and add new data frame
    colnames(temptpdf)[3:ncol(temptpdf)]=paste("TP",temptp,"_",colnames(temptpdf)[3:ncol(temptpdf)], sep="")
    dftps[[ntp]] <- temptpdf
  }
  ftpdf_ewm_lg_exp=Reduce(function(x,y)merge(x,y,by=c("ID","CohortAssgn"),all.x=TRUE), dftps)
  
  ### Removing errors for now 
  
  # delids=c("LDS0110021","LDS0370163","LDS0730044","LDS0730055","LDS0730150","LDS1770107","LDS1770112","LDS1770181","LDS3600187")
  # ftpdf_ewm_lg=subset(ftpdf_ewm_lg, !ftpdf_ewm_lg$ID %in% delids)
  
  
  dbname=paste("ready_datasets/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FTP_LongitudinalDatabaseErodedWM.xlsx",sep = "")
  write.xlsx(ftpdf_ewm_lg, dbname, row.names=FALSE, col.names = TRUE)
  
  ### Added module to get other longitudinal plot type ###
  
  iddobftp=ftpdf_bl[c(1,249,252,256)]
  ftpdf_ewm_lg_plot2=merge(ftpdf_ewm_lg,iddobftp,by=c("ID"))
  ftpdf_ewm_lg_plot2$age=time_length(difftime(ftpdf_ewm_lg_plot2$FTPPET_Date, ftpdf_ewm_lg_plot2$dob), "years")
  ftpdf_ewm_lg_plot2$CDRtot=ifelse(is.na(ftpdf_ewm_lg_plot2$CDRtot),"NaN",as.character(ftpdf_ewm_lg_plot2$CDRtot))
  ftpdf_ewm_lg_plot2$CDRtot=factor(ftpdf_ewm_lg_plot2$CDRtot, levels=c("0.5","1","NaN"))
  
  ftplg_plot2=ggplot(subset(ftpdf_ewm_lg_plot2,ftpdf_ewm_lg_plot2$CohortAssgn=="EOAD"), aes(x=age, y=MetaROI_MRIBASED_SUVR, group=ID, colour = CohortAssgn)) +
    geom_line(alpha=1, size=1) +
    ylab("18F-Flortaucipir-PET metaROI SUVR (Ref: ErodedWM)") +xlab("Age") + theme_minimal() +
    theme(axis.title.x=element_text(size=16), axis.text.x = element_text(size=14), axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=14)) +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF")) + labs(title="Longitudinal FTP-PET metaROI SUVR") +
    theme(panel.grid.minor = element_blank(), legend.position = "none") + scale_y_continuous(limits=c(0.75,2.5))
  
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongFTP_EWM_Plot2.jpg",sep = "")
  ggsave(dbname, plot=ftplg_plot2, width = 18, height = 15, units = "cm")
  
  ftplg_plot2_apoe=ggplot(subset(ftpdf_ewm_lg_plot2, ftpdf_ewm_lg_plot2$CohortAssgn=="EOAD"), aes(x=age, y=MetaROI_MRIBASED_SUVR, group=ID, colour = APOEpos)) +
    geom_line(alpha=1, size=1) +
    ylab("18F-Flortaucipir-PET metaROI SUVR (Ref: ErodedWM)") +xlab("Age") + theme_minimal() + scale_color_wsj() +
    theme(axis.title.x=element_text(size=16), axis.text.x = element_text(size=14), axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=14)) + labs(title="Longitudinal FTP-PET metaROI SUVR") +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom") + scale_y_continuous(limits=c(0.75,2.5))
  
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongFTP_EWM_Plot2_APOE.jpg",sep = "")
  ggsave(dbname, plot=ftplg_plot2_apoe, width = 18, height = 15, units = "cm")
  
  
  ftplg_plot2_cdr=ggplot(subset(ftpdf_ewm_lg_plot2, ftpdf_ewm_lg_plot2$CohortAssgn=="EOAD"), aes(x=age, y=MetaROI_MRIBASED_SUVR, group=ID, colour = CDRtot)) +
    geom_line(alpha=1, size=1) +
    ylab("18F-Flortaucipir-PET metaROI SUVR (Ref: ErodedWM)") +xlab("Age") + theme_minimal() + scale_color_manual(values=c("#716FB2","#EB093C","#B4B9Bf")) +
    theme(axis.title.x=element_text(size=16), axis.text.x = element_text(size=14), axis.text.y=element_text(size=14, face="bold"), axis.title.y=element_text(size=14)) + labs(title="Longitudinal FTP-PET metaROI SUVR") +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom") + scale_y_continuous(limits=c(0.75,2.5))
  
  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongFTP_EWM_Plot2_CDR.jpg",sep = "")
  ggsave(dbname, plot=fbblg_plot2_cdr, width = 18, height = 15, units = "cm")
  
  
  
  
  
  
  
  
  
  #### Report: ggbrain summary images ####

  fbbdf_bl_eoad=subset(fbbdf_bl, fbbdf_bl$CohortAssgn=="EOAD")
  ggseg_fbb_eoad=fbbdf_bl_eoad[65:133]
  ggseg_fbb_eoad=describe(ggseg_fbb_eoad)
  ggseg_fbb_eoad=data.frame(row.names(ggseg_fbb_eoad),ggseg_fbb_eoad$mean)
  colnames(ggseg_fbb_eoad)=c("ROI","em")
  ggseg_fbb_eoad$hemi=ifelse(grepl("lh",ggseg_fbb_eoad$ROI)==TRUE,"left","right")
  ggseg_fbb_eoad$ROI=gsub("ctx_lh_","",ggseg_fbb_eoad$ROI)
  ggseg_fbb_eoad$ROI=gsub("ctx_rh_","",ggseg_fbb_eoad$ROI)
  ggseg_fbb_eoad$ROI=gsub("_MRIBASED_SUVR","",ggseg_fbb_eoad$ROI)

  ggseg_fbb_eoad$region=ifelse(ggseg_fbb_eoad$ROI=="caudalanteriorcingulate","caudal anterior cingulate",
                               ifelse(ggseg_fbb_eoad$ROI=="caudalmiddlefrontal","caudal middle frontal",
                                      ifelse(ggseg_fbb_eoad$ROI=="inferiorparietal","inferior parietal",
                                             ifelse(ggseg_fbb_eoad$ROI=="inferiortemporal","inferior temporal",
                                                    ifelse(ggseg_fbb_eoad$ROI=="isthmuscingulate","isthmus cingulate",
                                                           ifelse(ggseg_fbb_eoad$ROI=="lateraloccipital","lateral occipital",
                                                                  ifelse(ggseg_fbb_eoad$ROI=="lateralorbitofrontal","lateral orbitofrontal",
                                                                         ifelse(ggseg_fbb_eoad$ROI=="medialorbitofrontal","medial orbitofrontal",
                                                                                ifelse(ggseg_fbb_eoad$ROI=="middletemporal","middle temporal",
                                                                                       ifelse(ggseg_fbb_eoad$ROI=="parsopercularis","pars opercularis",
                                                                                              ifelse(ggseg_fbb_eoad$ROI=="parsorbitalis","pars orbitalis",
                                                                                                     ifelse(ggseg_fbb_eoad$ROI=="parstriangularis","pars triangularis",
                                                                                                            ifelse(ggseg_fbb_eoad$ROI=="posteriorcingulate","posterior cingulate",
                                                                                                                   ifelse(ggseg_fbb_eoad$ROI=="rostralanteriorcingulate","rostral anterior cingulate",
                                                                                                                          ifelse(ggseg_fbb_eoad$ROI=="rostralmiddlefrontal","rostral middle frontal",
                                                                                                                                 ifelse(ggseg_fbb_eoad$ROI=="superiorfrontal","superior frontal",
                                                                                                                                        ifelse(ggseg_fbb_eoad$ROI=="superiorparietal","superior parietal",
                                                                                                                                               ifelse(ggseg_fbb_eoad$ROI=="superiortemporal","superior temporal",
                                                                                                                                                      ifelse(ggseg_fbb_eoad$ROI=="temporalpole","temporal pole",
                                                                                                                                                             ifelse(ggseg_fbb_eoad$ROI=="frontalpole","frontal pole",
                                                                                                                                                                    ifelse(ggseg_fbb_eoad$ROI=="transversetemporal","transverse temporal",as.character(ggseg_fbb_eoad$ROI))))))))))))))))))))))

  ggseg_fbb_eoad=ggseg_fbb_eoad[c(4,2,3)]
  ggseg_fbb_eoad=subset(ggseg_fbb_eoad, !ggseg_fbb_eoad$region=="unknown")

  p1=ggseg_fbb_eoad %>%
    ggseg(col="black",size=.1, mapping=aes(fill=em), position="stacked")+ scale_fill_viridis_c(na.value="grey95", limits=c(0.5,2)) +
    theme(legend.position="right") + labs(title="FBB-PET in EOAD", fill="SUVR") +theme(text=element_text(size=12),axis.text.x = element_blank(), axis.title.x = element_blank(),
                                                                                       axis.text.y=element_text(size=10)) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))
  p_fbb_eoad=p1

  fbbdf_bl_eononad=subset(fbbdf_bl, fbbdf_bl$CohortAssgn=="EOnonAD")
  ggseg_fbb_eononad=fbbdf_bl_eononad[65:133]
  ggseg_fbb_eononad=describe(ggseg_fbb_eononad)
  ggseg_fbb_eononad=data.frame(row.names(ggseg_fbb_eononad),ggseg_fbb_eononad$mean)
  colnames(ggseg_fbb_eononad)=c("ROI","em")
  ggseg_fbb_eononad$hemi=ifelse(grepl("lh",ggseg_fbb_eononad$ROI)==TRUE,"left","right")
  ggseg_fbb_eononad$ROI=gsub("ctx_lh_","",ggseg_fbb_eononad$ROI)
  ggseg_fbb_eononad$ROI=gsub("ctx_rh_","",ggseg_fbb_eononad$ROI)
  ggseg_fbb_eononad$ROI=gsub("_MRIBASED_SUVR","",ggseg_fbb_eononad$ROI)

  ggseg_fbb_eononad$region=ifelse(ggseg_fbb_eononad$ROI=="caudalanteriorcingulate","caudal anterior cingulate",
                                  ifelse(ggseg_fbb_eononad$ROI=="caudalmiddlefrontal","caudal middle frontal",
                                         ifelse(ggseg_fbb_eononad$ROI=="inferiorparietal","inferior parietal",
                                                ifelse(ggseg_fbb_eononad$ROI=="inferiortemporal","inferior temporal",
                                                       ifelse(ggseg_fbb_eononad$ROI=="isthmuscingulate","isthmus cingulate",
                                                              ifelse(ggseg_fbb_eononad$ROI=="lateraloccipital","lateral occipital",
                                                                     ifelse(ggseg_fbb_eononad$ROI=="lateralorbitofrontal","lateral orbitofrontal",
                                                                            ifelse(ggseg_fbb_eononad$ROI=="medialorbitofrontal","medial orbitofrontal",
                                                                                   ifelse(ggseg_fbb_eononad$ROI=="middletemporal","middle temporal",
                                                                                          ifelse(ggseg_fbb_eononad$ROI=="parsopercularis","pars opercularis",
                                                                                                 ifelse(ggseg_fbb_eononad$ROI=="parsorbitalis","pars orbitalis",
                                                                                                        ifelse(ggseg_fbb_eononad$ROI=="parstriangularis","pars triangularis",
                                                                                                               ifelse(ggseg_fbb_eononad$ROI=="posteriorcingulate","posterior cingulate",
                                                                                                                      ifelse(ggseg_fbb_eononad$ROI=="rostralanteriorcingulate","rostral anterior cingulate",
                                                                                                                             ifelse(ggseg_fbb_eononad$ROI=="rostralmiddlefrontal","rostral middle frontal",
                                                                                                                                    ifelse(ggseg_fbb_eononad$ROI=="superiorfrontal","superior frontal",
                                                                                                                                           ifelse(ggseg_fbb_eononad$ROI=="superiorparietal","superior parietal",
                                                                                                                                                  ifelse(ggseg_fbb_eononad$ROI=="superiortemporal","superior temporal",
                                                                                                                                                         ifelse(ggseg_fbb_eononad$ROI=="temporalpole","temporal pole",
                                                                                                                                                                ifelse(ggseg_fbb_eononad$ROI=="frontalpole","frontal pole",
                                                                                                                                                                       ifelse(ggseg_fbb_eononad$ROI=="transversetemporal","transverse temporal",as.character(ggseg_fbb_eononad$ROI))))))))))))))))))))))

  ggseg_fbb_eononad=ggseg_fbb_eononad[c(4,2,3)]
  ggseg_fbb_eononad=subset(ggseg_fbb_eononad, !ggseg_fbb_eononad$region=="unknown")


  p1=ggseg_fbb_eononad %>%
    ggseg(col="black",size=.1, mapping=aes(fill=em), position="stacked")+ scale_fill_viridis_c(na.value="grey95", limits=c(0.5,2)) +
    theme(legend.position="right") + labs(title="FBB-PET in EOnonAD", fill="SUVR") +theme(text=element_text(size=12),axis.text.x = element_blank(), axis.title.x = element_blank(),
                                                                                          axis.text.y=element_text(size=10)) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))

  p_fbb_eononad=p1

  fbbdf_bl_cn=subset(fbbdf_bl, fbbdf_bl$CohortAssgn=="CN")
  ggseg_fbb_cn=fbbdf_bl_cn[65:133]
  ggseg_fbb_cn=describe(ggseg_fbb_cn)
  ggseg_fbb_cn=data.frame(row.names(ggseg_fbb_cn),ggseg_fbb_cn$mean)
  colnames(ggseg_fbb_cn)=c("ROI","em")
  ggseg_fbb_cn$hemi=ifelse(grepl("lh",ggseg_fbb_cn$ROI)==TRUE,"left","right")
  ggseg_fbb_cn$ROI=gsub("ctx_lh_","",ggseg_fbb_cn$ROI)
  ggseg_fbb_cn$ROI=gsub("ctx_rh_","",ggseg_fbb_cn$ROI)
  ggseg_fbb_cn$ROI=gsub("_MRIBASED_SUVR","",ggseg_fbb_cn$ROI)

  ggseg_fbb_cn$region=ifelse(ggseg_fbb_cn$ROI=="caudalanteriorcingulate","caudal anterior cingulate",
                             ifelse(ggseg_fbb_cn$ROI=="caudalmiddlefrontal","caudal middle frontal",
                                    ifelse(ggseg_fbb_cn$ROI=="inferiorparietal","inferior parietal",
                                           ifelse(ggseg_fbb_cn$ROI=="inferiortemporal","inferior temporal",
                                                  ifelse(ggseg_fbb_cn$ROI=="isthmuscingulate","isthmus cingulate",
                                                         ifelse(ggseg_fbb_cn$ROI=="lateraloccipital","lateral occipital",
                                                                ifelse(ggseg_fbb_cn$ROI=="lateralorbitofrontal","lateral orbitofrontal",
                                                                       ifelse(ggseg_fbb_cn$ROI=="medialorbitofrontal","medial orbitofrontal",
                                                                              ifelse(ggseg_fbb_cn$ROI=="middletemporal","middle temporal",
                                                                                     ifelse(ggseg_fbb_cn$ROI=="parsopercularis","pars opercularis",
                                                                                            ifelse(ggseg_fbb_cn$ROI=="parsorbitalis","pars orbitalis",
                                                                                                   ifelse(ggseg_fbb_cn$ROI=="parstriangularis","pars triangularis",
                                                                                                          ifelse(ggseg_fbb_cn$ROI=="posteriorcingulate","posterior cingulate",
                                                                                                                 ifelse(ggseg_fbb_cn$ROI=="rostralanteriorcingulate","rostral anterior cingulate",
                                                                                                                        ifelse(ggseg_fbb_cn$ROI=="rostralmiddlefrontal","rostral middle frontal",
                                                                                                                               ifelse(ggseg_fbb_cn$ROI=="superiorfrontal","superior frontal",
                                                                                                                                      ifelse(ggseg_fbb_cn$ROI=="superiorparietal","superior parietal",
                                                                                                                                             ifelse(ggseg_fbb_cn$ROI=="superiortemporal","superior temporal",
                                                                                                                                                    ifelse(ggseg_fbb_cn$ROI=="temporalpole","temporal pole",
                                                                                                                                                           ifelse(ggseg_fbb_cn$ROI=="frontalpole","frontal pole",
                                                                                                                                                                  ifelse(ggseg_fbb_cn$ROI=="transversetemporal","transverse temporal",as.character(ggseg_fbb_cn$ROI))))))))))))))))))))))

  ggseg_fbb_cn=ggseg_fbb_cn[c(4,2,3)]
  ggseg_fbb_cn=subset(ggseg_fbb_cn, !ggseg_fbb_cn$region=="unknown")


  p1=ggseg_fbb_cn %>%
    ggseg(col="black",size=.1, mapping=aes(fill=em), position="stacked")+ scale_fill_viridis_c(na.value="grey95", limits=c(0.5,2)) +
    theme(legend.position="right") + labs(title="FBB-PET in CN", fill="SUVR") +theme(text=element_text(size=12),axis.text.x = element_blank(), axis.title.x = element_blank(),
                                                                                     axis.text.y=element_text(size=10)) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))

  p_fbb_cn=p1

  ftpdf_bl_eoad=subset(ftpdf_bl, ftpdf_bl$CohortAssgn=="EOAD")
  ggseg_ftp_eoad=ftpdf_bl_eoad[58:126]
  ggseg_ftp_eoad=describe(ggseg_ftp_eoad)
  ggseg_ftp_eoad=data.frame(row.names(ggseg_ftp_eoad),ggseg_ftp_eoad$mean)
  colnames(ggseg_ftp_eoad)=c("ROI","em")
  ggseg_ftp_eoad$hemi=ifelse(grepl("lh",ggseg_ftp_eoad$ROI)==TRUE,"left","right")
  ggseg_ftp_eoad$ROI=gsub("ctx_lh_","",ggseg_ftp_eoad$ROI)
  ggseg_ftp_eoad$ROI=gsub("ctx_rh_","",ggseg_ftp_eoad$ROI)
  ggseg_ftp_eoad$ROI=gsub("_MRIBASED_SUVR","",ggseg_ftp_eoad$ROI)

  ggseg_ftp_eoad$region=ifelse(ggseg_ftp_eoad$ROI=="caudalanteriorcingulate","caudal anterior cingulate",
                               ifelse(ggseg_ftp_eoad$ROI=="caudalmiddlefrontal","caudal middle frontal",
                                      ifelse(ggseg_ftp_eoad$ROI=="inferiorparietal","inferior parietal",
                                             ifelse(ggseg_ftp_eoad$ROI=="inferiortemporal","inferior temporal",
                                                    ifelse(ggseg_ftp_eoad$ROI=="isthmuscingulate","isthmus cingulate",
                                                           ifelse(ggseg_ftp_eoad$ROI=="lateraloccipital","lateral occipital",
                                                                  ifelse(ggseg_ftp_eoad$ROI=="lateralorbitofrontal","lateral orbitofrontal",
                                                                         ifelse(ggseg_ftp_eoad$ROI=="medialorbitofrontal","medial orbitofrontal",
                                                                                ifelse(ggseg_ftp_eoad$ROI=="middletemporal","middle temporal",
                                                                                       ifelse(ggseg_ftp_eoad$ROI=="parsopercularis","pars opercularis",
                                                                                              ifelse(ggseg_ftp_eoad$ROI=="parsorbitalis","pars orbitalis",
                                                                                                     ifelse(ggseg_ftp_eoad$ROI=="parstriangularis","pars triangularis",
                                                                                                            ifelse(ggseg_ftp_eoad$ROI=="posteriorcingulate","posterior cingulate",
                                                                                                                   ifelse(ggseg_ftp_eoad$ROI=="rostralanteriorcingulate","rostral anterior cingulate",
                                                                                                                          ifelse(ggseg_ftp_eoad$ROI=="rostralmiddlefrontal","rostral middle frontal",
                                                                                                                                 ifelse(ggseg_ftp_eoad$ROI=="superiorfrontal","superior frontal",
                                                                                                                                        ifelse(ggseg_ftp_eoad$ROI=="superiorparietal","superior parietal",
                                                                                                                                               ifelse(ggseg_ftp_eoad$ROI=="superiortemporal","superior temporal",
                                                                                                                                                      ifelse(ggseg_ftp_eoad$ROI=="temporalpole","temporal pole",
                                                                                                                                                             ifelse(ggseg_ftp_eoad$ROI=="frontalpole","frontal pole",
                                                                                                                                                                    ifelse(ggseg_ftp_eoad$ROI=="transversetemporal","transverse temporal",as.character(ggseg_ftp_eoad$ROI))))))))))))))))))))))

  ggseg_ftp_eoad=ggseg_ftp_eoad[c(4,2,3)]
  ggseg_ftp_eoad=subset(ggseg_ftp_eoad, !ggseg_ftp_eoad$region=="unknown")

  p1=ggseg_ftp_eoad %>%
    ggseg(col="black",size=.1, mapping=aes(fill=em), position="stacked")+ scale_fill_viridis_c(option="inferno",na.value="grey95", limits=c(0.5,2.8)) +
    theme(legend.position="right") + labs(title="FTP-PET in EOAD", fill="SUVR") +theme(text=element_text(size=12),axis.text.x = element_blank(), axis.title.x = element_blank(),
                                                                                       axis.text.y=element_text(size=10)) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))
  p_ftp_eoad=p1

  ftpdf_bl_eononad=subset(ftpdf_bl, ftpdf_bl$CohortAssgn=="EOnonAD")
  ggseg_ftp_eononad=ftpdf_bl_eononad[58:126]
  ggseg_ftp_eononad=describe(ggseg_ftp_eononad)
  ggseg_ftp_eononad=data.frame(row.names(ggseg_ftp_eononad),ggseg_ftp_eononad$mean)
  colnames(ggseg_ftp_eononad)=c("ROI","em")
  ggseg_ftp_eononad$hemi=ifelse(grepl("lh",ggseg_ftp_eononad$ROI)==TRUE,"left","right")
  ggseg_ftp_eononad$ROI=gsub("ctx_lh_","",ggseg_ftp_eononad$ROI)
  ggseg_ftp_eononad$ROI=gsub("ctx_rh_","",ggseg_ftp_eononad$ROI)
  ggseg_ftp_eononad$ROI=gsub("_MRIBASED_SUVR","",ggseg_ftp_eononad$ROI)

  ggseg_ftp_eononad$region=ifelse(ggseg_ftp_eononad$ROI=="caudalanteriorcingulate","caudal anterior cingulate",
                                  ifelse(ggseg_ftp_eononad$ROI=="caudalmiddlefrontal","caudal middle frontal",
                                         ifelse(ggseg_ftp_eononad$ROI=="inferiorparietal","inferior parietal",
                                                ifelse(ggseg_ftp_eononad$ROI=="inferiortemporal","inferior temporal",
                                                       ifelse(ggseg_ftp_eononad$ROI=="isthmuscingulate","isthmus cingulate",
                                                              ifelse(ggseg_ftp_eononad$ROI=="lateraloccipital","lateral occipital",
                                                                     ifelse(ggseg_ftp_eononad$ROI=="lateralorbitofrontal","lateral orbitofrontal",
                                                                            ifelse(ggseg_ftp_eononad$ROI=="medialorbitofrontal","medial orbitofrontal",
                                                                                   ifelse(ggseg_ftp_eononad$ROI=="middletemporal","middle temporal",
                                                                                          ifelse(ggseg_ftp_eononad$ROI=="parsopercularis","pars opercularis",
                                                                                                 ifelse(ggseg_ftp_eononad$ROI=="parsorbitalis","pars orbitalis",
                                                                                                        ifelse(ggseg_ftp_eononad$ROI=="parstriangularis","pars triangularis",
                                                                                                               ifelse(ggseg_ftp_eononad$ROI=="posteriorcingulate","posterior cingulate",
                                                                                                                      ifelse(ggseg_ftp_eononad$ROI=="rostralanteriorcingulate","rostral anterior cingulate",
                                                                                                                             ifelse(ggseg_ftp_eononad$ROI=="rostralmiddlefrontal","rostral middle frontal",
                                                                                                                                    ifelse(ggseg_ftp_eononad$ROI=="superiorfrontal","superior frontal",
                                                                                                                                           ifelse(ggseg_ftp_eononad$ROI=="superiorparietal","superior parietal",
                                                                                                                                                  ifelse(ggseg_ftp_eononad$ROI=="superiortemporal","superior temporal",
                                                                                                                                                         ifelse(ggseg_ftp_eononad$ROI=="temporalpole","temporal pole",
                                                                                                                                                                ifelse(ggseg_ftp_eononad$ROI=="frontalpole","frontal pole",
                                                                                                                                                                       ifelse(ggseg_ftp_eononad$ROI=="transversetemporal","transverse temporal",as.character(ggseg_ftp_eononad$ROI))))))))))))))))))))))

  ggseg_ftp_eononad=ggseg_ftp_eononad[c(4,2,3)]
  ggseg_ftp_eononad=subset(ggseg_ftp_eononad, !ggseg_ftp_eononad$region=="unknown")


  p1=ggseg_ftp_eononad %>%
    ggseg(col="black",size=.1, mapping=aes(fill=em), position="stacked")+ scale_fill_viridis_c(option="inferno",na.value="grey95", limits=c(0.5,2.8)) +
    theme(legend.position="right") + labs(title="FTP-PET in EOnonAD", fill="SUVR") +theme(text=element_text(size=12),axis.text.x = element_blank(), axis.title.x = element_blank(),
                                                                                          axis.text.y=element_text(size=10)) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))

  p_ftp_eononad=p1

  ftpdf_bl_cn=subset(ftpdf_bl, ftpdf_bl$CohortAssgn=="CN")
  ggseg_ftp_cn=ftpdf_bl_cn[58:126]
  ggseg_ftp_cn=describe(ggseg_ftp_cn)
  ggseg_ftp_cn=data.frame(row.names(ggseg_ftp_cn),ggseg_ftp_cn$mean)
  colnames(ggseg_ftp_cn)=c("ROI","em")
  ggseg_ftp_cn$hemi=ifelse(grepl("lh",ggseg_ftp_cn$ROI)==TRUE,"left","right")
  ggseg_ftp_cn$ROI=gsub("ctx_lh_","",ggseg_ftp_cn$ROI)
  ggseg_ftp_cn$ROI=gsub("ctx_rh_","",ggseg_ftp_cn$ROI)
  ggseg_ftp_cn$ROI=gsub("_MRIBASED_SUVR","",ggseg_ftp_cn$ROI)

  ggseg_ftp_cn$region=ifelse(ggseg_ftp_cn$ROI=="caudalanteriorcingulate","caudal anterior cingulate",
                             ifelse(ggseg_ftp_cn$ROI=="caudalmiddlefrontal","caudal middle frontal",
                                    ifelse(ggseg_ftp_cn$ROI=="inferiorparietal","inferior parietal",
                                           ifelse(ggseg_ftp_cn$ROI=="inferiortemporal","inferior temporal",
                                                  ifelse(ggseg_ftp_cn$ROI=="isthmuscingulate","isthmus cingulate",
                                                         ifelse(ggseg_ftp_cn$ROI=="lateraloccipital","lateral occipital",
                                                                ifelse(ggseg_ftp_cn$ROI=="lateralorbitofrontal","lateral orbitofrontal",
                                                                       ifelse(ggseg_ftp_cn$ROI=="medialorbitofrontal","medial orbitofrontal",
                                                                              ifelse(ggseg_ftp_cn$ROI=="middletemporal","middle temporal",
                                                                                     ifelse(ggseg_ftp_cn$ROI=="parsopercularis","pars opercularis",
                                                                                            ifelse(ggseg_ftp_cn$ROI=="parsorbitalis","pars orbitalis",
                                                                                                   ifelse(ggseg_ftp_cn$ROI=="parstriangularis","pars triangularis",
                                                                                                          ifelse(ggseg_ftp_cn$ROI=="posteriorcingulate","posterior cingulate",
                                                                                                                 ifelse(ggseg_ftp_cn$ROI=="rostralanteriorcingulate","rostral anterior cingulate",
                                                                                                                        ifelse(ggseg_ftp_cn$ROI=="rostralmiddlefrontal","rostral middle frontal",
                                                                                                                               ifelse(ggseg_ftp_cn$ROI=="superiorfrontal","superior frontal",
                                                                                                                                      ifelse(ggseg_ftp_cn$ROI=="superiorparietal","superior parietal",
                                                                                                                                             ifelse(ggseg_ftp_cn$ROI=="superiortemporal","superior temporal",
                                                                                                                                                    ifelse(ggseg_ftp_cn$ROI=="temporalpole","temporal pole",
                                                                                                                                                           ifelse(ggseg_ftp_cn$ROI=="frontalpole","frontal pole",
                                                                                                                                                                  ifelse(ggseg_ftp_cn$ROI=="transversetemporal","transverse temporal",as.character(ggseg_ftp_cn$ROI))))))))))))))))))))))

  ggseg_ftp_cn=ggseg_ftp_cn[c(4,2,3)]
  ggseg_ftp_cn=subset(ggseg_ftp_cn, !ggseg_ftp_cn$region=="unknown")


  p1=ggseg_ftp_cn %>%
    ggseg(col="black",size=.1, mapping=aes(fill=em), position="stacked")+ scale_fill_viridis_c(option="inferno",na.value="grey95", limits=c(0.5,2.8)) +
    theme(legend.position="right") + labs(title="FTP-PET in CN", fill="SUVR") +theme(text=element_text(size=12),axis.text.x = element_blank(), axis.title.x = element_blank(),
                                                                                     axis.text.y=element_text(size=10)) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))

  p_ftp_cn=p1


  #### Report: timetrends ####

  ## FBB-PET

  timetrends_fbb_bl=fbbdf_bl[c(1,2,9)]
  timetrends_fbb_bl_cn=subset(timetrends_fbb_bl, timetrends_fbb_bl$CohortAssgn=="CN")
  timetrends_fbb_bl_eoad=subset(timetrends_fbb_bl, timetrends_fbb_bl$CohortAssgn=="EOAD")
  timetrends_fbb_bl_eononad=subset(timetrends_fbb_bl, timetrends_fbb_bl$CohortAssgn=="EOnonAD")

  timetrends_fbb_bl_eononad=as.data.frame(table(timetrends_fbb_bl_eononad$FBBPET_Date))
  timetrends_fbb_bl_eononad$cumsum=cumsum(timetrends_fbb_bl_eononad$Freq)
  timetrends_fbb_bl_eononad$Cohort=rep("EOnonAD",nrow(timetrends_fbb_bl_eononad))
  timetrends_fbb_bl_eononad=timetrends_fbb_bl_eononad[order(timetrends_fbb_bl_eononad$Var1),]

  timetrends_fbb_bl_eoad=as.data.frame(table(timetrends_fbb_bl_eoad$FBBPET_Date))
  timetrends_fbb_bl_eoad$cumsum=cumsum(timetrends_fbb_bl_eoad$Freq)
  timetrends_fbb_bl_eoad$Cohort=rep("EOAD",nrow(timetrends_fbb_bl_eoad))
  timetrends_fbb_bl_eoad=timetrends_fbb_bl_eoad[order(timetrends_fbb_bl_eoad$Var1),]

  timetrends_fbb_bl_cn=as.data.frame(table(timetrends_fbb_bl_cn$FBBPET_Date))
  timetrends_fbb_bl_cn$cumsum=cumsum(timetrends_fbb_bl_cn$Freq)
  timetrends_fbb_bl_cn$Cohort=rep("CN",nrow(timetrends_fbb_bl_cn))
  timetrends_fbb_bl_cn=timetrends_fbb_bl_cn[order(timetrends_fbb_bl_cn$Var1),]

  timetrends_fbb_final=rbind(timetrends_fbb_bl_cn, timetrends_fbb_bl_eoad, timetrends_fbb_bl_eononad)
  timetrends_fbb_final$Var1=as.Date(as.character(timetrends_fbb_final$Var1, format = "%m/%d/%Y"))
  timetrends_fbb_final=timetrends_fbb_final[order(timetrends_fbb_final$Var1),]
  timetrends_fbb_final$Cohort=factor(timetrends_fbb_final$Cohort, levels=c("EOAD","EOnonAD","CN"))

  labsttdf_fbb=rbind(timetrends_fbb_bl_cn[nrow(timetrends_fbb_bl_cn),],timetrends_fbb_bl_eoad[nrow(timetrends_fbb_bl_eoad),],timetrends_fbb_bl_eononad[nrow(timetrends_fbb_bl_eononad),])
  labsttdf_fbb$Cohort=factor(labsttdf_fbb$Cohort, levels=c("EOAD","EOnonAD","CN"))
  labsttdf_fbb$Var1=as.Date(as.character(labsttdf_fbb$Var1, format = "%m/%d/%Y"))


  p_timetrend_fbb=ggplot(data=timetrends_fbb_final, aes(x=Var1, y=cumsum, group=Cohort,color=Cohort)) + geom_line(size=1.5) +
    ylab("Number of baseline FBB-PET scans processed") +
    xlab("Date") + geom_text(data = labsttdf_fbb,
                             aes(label = cumsum),
                             nudge_x =60, nudge_y=5, size=4) +
    theme_minimal() +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10))  +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

  ## FTP-PET

  timetrends_ftp_bl=ftpdf_bl[c(1,2,4)]
  timetrends_ftp_bl_cn=subset(timetrends_ftp_bl, timetrends_ftp_bl$CohortAssgn=="CN")
  timetrends_ftp_bl_eoad=subset(timetrends_ftp_bl, timetrends_ftp_bl$CohortAssgn=="EOAD")
  timetrends_ftp_bl_eononad=subset(timetrends_ftp_bl, timetrends_ftp_bl$CohortAssgn=="EOnonAD")

  timetrends_ftp_bl_eononad=as.data.frame(table(timetrends_ftp_bl_eononad$FTPPET_Date))
  timetrends_ftp_bl_eononad$cumsum=cumsum(timetrends_ftp_bl_eononad$Freq)
  timetrends_ftp_bl_eononad$Cohort=rep("EOnonAD",nrow(timetrends_ftp_bl_eononad))
  timetrends_ftp_bl_eononad=timetrends_ftp_bl_eononad[order(timetrends_ftp_bl_eononad$Var1),]

  timetrends_ftp_bl_eoad=as.data.frame(table(timetrends_ftp_bl_eoad$FTPPET_Date))
  timetrends_ftp_bl_eoad$cumsum=cumsum(timetrends_ftp_bl_eoad$Freq)
  timetrends_ftp_bl_eoad$Cohort=rep("EOAD",nrow(timetrends_ftp_bl_eoad))
  timetrends_ftp_bl_eoad=timetrends_ftp_bl_eoad[order(timetrends_ftp_bl_eoad$Var1),]

  timetrends_ftp_bl_cn=as.data.frame(table(timetrends_ftp_bl_cn$FTPPET_Date))
  timetrends_ftp_bl_cn$cumsum=cumsum(timetrends_ftp_bl_cn$Freq)
  timetrends_ftp_bl_cn$Cohort=rep("CN",nrow(timetrends_ftp_bl_cn))
  timetrends_ftp_bl_cn=timetrends_ftp_bl_cn[order(timetrends_ftp_bl_cn$Var1),]

  timetrends_ftp_final=rbind(timetrends_ftp_bl_cn, timetrends_ftp_bl_eoad, timetrends_ftp_bl_eononad)
  timetrends_ftp_final$Var1=as.Date(as.character(timetrends_ftp_final$Var1, format = "%m/%d/%Y"))
  timetrends_ftp_final=timetrends_ftp_final[order(timetrends_ftp_final$Var1),]
  timetrends_ftp_final$Cohort=factor(timetrends_ftp_final$Cohort, levels=c("EOAD","EOnonAD","CN"))

  labsttdf_ftp=rbind(timetrends_ftp_bl_cn[nrow(timetrends_ftp_bl_cn),],timetrends_ftp_bl_eoad[nrow(timetrends_ftp_bl_eoad),],timetrends_ftp_bl_eononad[nrow(timetrends_ftp_bl_eononad),])
  labsttdf_ftp$Cohort=factor(labsttdf_ftp$Cohort, levels=c("EOAD","EOnonAD","CN"))
  labsttdf_ftp$Var1=as.Date(as.character(labsttdf_ftp$Var1, format = "%m/%d/%Y"))


  p_timetrend_ftp=ggplot(data=timetrends_ftp_final, aes(x=Var1, y=cumsum, group=Cohort,color=Cohort)) + geom_line(size=1.5) +
    ylab("Number of baseline FTP-PET scans processed") +
    xlab("Date") + geom_text(data = labsttdf_ftp,
                             aes(label = cumsum),
                             nudge_x =60, nudge_y=5, size=4) +
    theme_minimal() +
    theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10))  +
    scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

## FDG-PET

timetrends_fdg_bl=fdgdf[c(1,2,4)]
timetrends_fdg_bl_cn=subset(timetrends_fdg_bl, timetrends_fdg_bl$CohortAssgn=="CN")
timetrends_fdg_bl_eoad=subset(timetrends_fdg_bl, timetrends_fdg_bl$CohortAssgn=="EOAD")
timetrends_fdg_bl_eononad=subset(timetrends_fdg_bl, timetrends_fdg_bl$CohortAssgn=="EOnonAD")

timetrends_fdg_bl_eononad=as.data.frame(table(timetrends_fdg_bl_eononad$FDGPET_Date))
timetrends_fdg_bl_eononad$cumsum=cumsum(timetrends_fdg_bl_eononad$Freq)
timetrends_fdg_bl_eononad$Cohort=rep("EOnonAD",nrow(timetrends_fdg_bl_eononad))
timetrends_fdg_bl_eononad=timetrends_fdg_bl_eononad[order(timetrends_fdg_bl_eononad$Var1),]

timetrends_fdg_bl_cn=as.data.frame(table(timetrends_fdg_bl_cn$FDGPET_Date))
timetrends_fdg_bl_cn$cumsum=cumsum(timetrends_fdg_bl_cn$Freq)
timetrends_fdg_bl_cn$Cohort=rep("CN",nrow(timetrends_fdg_bl_cn))
timetrends_fdg_bl_cn=timetrends_fdg_bl_cn[order(timetrends_fdg_bl_cn$Var1),]

timetrends_fdg_final=rbind(timetrends_fdg_bl_cn, timetrends_fdg_bl_eoad, timetrends_fdg_bl_eononad)
timetrends_fdg_final$Var1=as.Date(as.character(timetrends_fdg_final$Var1, format = "%m/%d/%Y"))
timetrends_fdg_final=timetrends_fdg_final[order(timetrends_fdg_final$Var1),]
timetrends_fdg_final$Cohort=factor(timetrends_fdg_final$Cohort, levels=c("EOAD","EOnonAD","CN"))

labsttdf_fdg=rbind(timetrends_fdg_bl_cn[nrow(timetrends_fdg_bl_cn),],timetrends_fdg_bl_eoad[nrow(timetrends_fdg_bl_eoad),],timetrends_fdg_bl_eononad[nrow(timetrends_fdg_bl_eononad),])
labsttdf_fdg$Cohort=factor(labsttdf_fdg$Cohort, levels=c("EOAD","EOnonAD","CN"))
labsttdf_fdg$Var1=as.Date(as.character(labsttdf_fdg$Var1, format = "%m/%d/%Y"))

p_timetrend_fdg=ggplot(data=timetrends_fdg_final, aes(x=Var1, y=cumsum, group=Cohort,color=Cohort)) + geom_line(size=1.5) +
  ylab("Number of FDG-PET scans processed") +
  xlab("Date") + geom_text(data = labsttdf_fdg,
                           aes(label = cumsum),
                           nudge_x =60, nudge_y=5, size=4)  + expand_limits(x=as.Date(as.character("2018-05-01", format = "%m/%d/%Y")), y=165) +
  theme_minimal() +
  theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10))  +
  scale_color_manual(values=c("#007CBE","#B4B9BF","#052049"))

#### Liana's Request: FBB-PET Plotting baseline vs. change automatically grabbing any Timepoint available ####

# Centiloids

fbblgnames=colnames(fbbdf_lg_cwm_exp)
fbbintids=grepl("Centiloids",fbblgnames)
tempfbbch=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbintids])
fbbdatids=grepl("Date",fbblgnames)
tempfbbchdat=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbdatids])
tempfbbch$Centiloids_Change=ifelse(is.na(tempfbbch$TP2_MRIBASED_Composite_Centiloids), 
                                   (tempfbbch$TP1_MRIBASED_Composite_Centiloids-tempfbbch$TP0_MRIBASED_Composite_Centiloids)/time_length(difftime(tempfbbchdat$TP1_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"), 
                                   (tempfbbch$TP2_MRIBASED_Composite_Centiloids-tempfbbch$TP0_MRIBASED_Composite_Centiloids)/time_length(difftime(tempfbbchdat$TP2_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"))

# tempfbbch$Centiloids_Change=ifelse(is.na(tempfbbch$TP2_MRIBASED_Composite_Centiloids), tempfbbch$TP1_MRIBASED_Composite_Centiloids-tempfbbch$TP0_MRIBASED_Composite_Centiloids, tempfbbch$TP2_MRIBASED_Composite_Centiloids-tempfbbch$TP0_MRIBASED_Composite_Centiloids)
tempfbbch$Centiloids_NTimepoints=ifelse(is.na(tempfbbch$TP2_MRIBASED_Composite_Centiloids),2,3)
tempfbbch=tempfbbch[c(1,ncol(tempfbbch)-1,ncol(tempfbbch))]
fbbdf_lg_cwm_exp=merge(fbbdf_lg_cwm_exp,tempfbbch,by="ID")

# Composite SUVR

fbbintids=grepl("Composite_SUVR",fbblgnames)
tempfbbch=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbintids])

fbbdatids=grepl("Date",fbblgnames)
tempfbbchdat=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbdatids])
tempfbbch$SUVR_Change=ifelse(is.na(tempfbbch$TP2_MRIBASED_Composite_SUVR_Type2), 
                             (tempfbbch$TP1_MRIBASED_Composite_SUVR_Type2-tempfbbch$TP0_MRIBASED_Composite_SUVR_Type2)/time_length(difftime(tempfbbchdat$TP1_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"), 
                             (tempfbbch$TP2_MRIBASED_Composite_SUVR_Type2-tempfbbch$TP0_MRIBASED_Composite_SUVR_Type2)/time_length(difftime(tempfbbchdat$TP2_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"))

# tempfbbch$SUVR_Change=ifelse(is.na(tempfbbch$TP2_MRIBASED_Composite_SUVR_Type2), tempfbbch$TP1_MRIBASED_Composite_SUVR_Type2-tempfbbch$TP0_MRIBASED_Composite_SUVR_Type2, tempfbbch$TP2_MRIBASED_Composite_SUVR_Type2-tempfbbch$TP0_MRIBASED_Composite_SUVR_Type2)
tempfbbch$SUVR_NTimepoints=ifelse(is.na(tempfbbch$TP2_MRIBASED_Composite_SUVR_Type2),2,3)
tempfbbch=tempfbbch[c(1,ncol(tempfbbch)-1,ncol(tempfbbch))]
fbbdf_lg_cwm_exp=merge(fbbdf_lg_cwm_exp,tempfbbch,by="ID")

# Frontal

fbblgnames=colnames(fbbdf_lg_cwm_exp)
fbbintids=grepl("Frontal",fbblgnames)
tempfbbch=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbintids])

fbbdatids=grepl("Date",fbblgnames)
tempfbbchdat=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbdatids])
tempfbbch$Frontal_Change=ifelse(is.na(tempfbbch$TP2_Frontal), 
                                (tempfbbch$TP1_Frontal-tempfbbch$TP0_Frontal)/time_length(difftime(tempfbbchdat$TP1_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"), 
                                (tempfbbch$TP2_Frontal-tempfbbch$TP0_Frontal)/time_length(difftime(tempfbbchdat$TP2_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"))


# tempfbbch$Frontal_Change=ifelse(is.na(tempfbbch$TP2_Frontal), tempfbbch$TP1_Frontal-tempfbbch$TP0_Frontal, tempfbbch$TP2_Frontal-tempfbbch$TP0_Frontal)
tempfbbch$Frontal_NTimepoints=ifelse(is.na(tempfbbch$TP2_Frontal),2,3)
tempfbbch=tempfbbch[c(1,ncol(tempfbbch)-1,ncol(tempfbbch))]
fbbdf_lg_cwm_exp=merge(fbbdf_lg_cwm_exp,tempfbbch,by="ID")

# Parietal

fbblgnames=colnames(fbbdf_lg_cwm_exp)
fbbintids=grepl("Parietal",fbblgnames)
tempfbbch=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbintids])

fbbdatids=grepl("Date",fbblgnames)
tempfbbchdat=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbdatids])
tempfbbch$Parietal_Change=ifelse(is.na(tempfbbch$TP2_Parietal), 
                                 (tempfbbch$TP1_Parietal-tempfbbch$TP0_Parietal)/time_length(difftime(tempfbbchdat$TP1_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"), 
                                 (tempfbbch$TP2_Parietal-tempfbbch$TP0_Parietal)/time_length(difftime(tempfbbchdat$TP2_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"))

# tempfbbch$Parietal_Change=ifelse(is.na(tempfbbch$TP2_Parietal), tempfbbch$TP1_Parietal-tempfbbch$TP0_Parietal, tempfbbch$TP2_Parietal-tempfbbch$TP0_Parietal)
tempfbbch$Parietal_NTimepoints=ifelse(is.na(tempfbbch$TP2_Parietal),2,3)
tempfbbch=tempfbbch[c(1,ncol(tempfbbch)-1,ncol(tempfbbch))]
fbbdf_lg_cwm_exp=merge(fbbdf_lg_cwm_exp,tempfbbch,by="ID")

# Occipital

fbblgnames=colnames(fbbdf_lg_cwm_exp)
fbbintids=grepl("Occipital",fbblgnames)
tempfbbch=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbintids])

fbbdatids=grepl("Date",fbblgnames)
tempfbbchdat=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbdatids])
tempfbbch$Occipital_Change=ifelse(is.na(tempfbbch$TP2_Occipital), 
                                  (tempfbbch$TP1_Occipital-tempfbbch$TP0_Occipital)/time_length(difftime(tempfbbchdat$TP1_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"), 
                                  (tempfbbch$TP2_Occipital-tempfbbch$TP0_Occipital)/time_length(difftime(tempfbbchdat$TP2_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"))


# tempfbbch$Occipital_Change=ifelse(is.na(tempfbbch$TP2_Occipital), tempfbbch$TP1_Occipital-tempfbbch$TP0_Occipital, tempfbbch$TP2_Occipital-tempfbbch$TP0_Occipital)
tempfbbch$Occipital_NTimepoints=ifelse(is.na(tempfbbch$TP2_Occipital),2,3)
tempfbbch=tempfbbch[c(1,ncol(tempfbbch)-1,ncol(tempfbbch))]
fbbdf_lg_cwm_exp=merge(fbbdf_lg_cwm_exp,tempfbbch,by="ID")

# Temporal

fbblgnames=colnames(fbbdf_lg_cwm_exp)
fbbintids=grepl("Temporal",fbblgnames)
tempfbbch=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbintids])

fbbdatids=grepl("Date",fbblgnames)
tempfbbchdat=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbdatids])
tempfbbch$Temporal_Change=ifelse(is.na(tempfbbch$TP2_Temporal), 
                                 (tempfbbch$TP1_Temporal-tempfbbch$TP0_Temporal)/time_length(difftime(tempfbbchdat$TP1_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"), 
                                 (tempfbbch$TP2_Temporal-tempfbbch$TP0_Temporal)/time_length(difftime(tempfbbchdat$TP2_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"))


# tempfbbch$Temporal_Change=ifelse(is.na(tempfbbch$TP2_Temporal), tempfbbch$TP1_Temporal-tempfbbch$TP0_Temporal, tempfbbch$TP2_Temporal-tempfbbch$TP0_Temporal)
tempfbbch$Temporal_NTimepoints=ifelse(is.na(tempfbbch$TP2_Temporal),2,3)
tempfbbch=tempfbbch[c(1,ncol(tempfbbch)-1,ncol(tempfbbch))]
fbbdf_lg_cwm_exp=merge(fbbdf_lg_cwm_exp,tempfbbch,by="ID")

# Cingulate

fbblgnames=colnames(fbbdf_lg_cwm_exp)
fbbintids=grepl("Cingulate",fbblgnames)
tempfbbch=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbintids])

fbbdatids=grepl("Date",fbblgnames)
tempfbbchdat=data.frame(fbbdf_lg_cwm_exp[1],fbbdf_lg_cwm_exp[,fbbdatids])
tempfbbch$Cingulate_Change=ifelse(is.na(tempfbbch$TP2_Cingulate), 
                                  (tempfbbch$TP1_Cingulate-tempfbbch$TP0_Cingulate)/time_length(difftime(tempfbbchdat$TP1_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"), 
                                  (tempfbbch$TP2_Cingulate-tempfbbch$TP0_Cingulate)/time_length(difftime(tempfbbchdat$TP2_FBBPET_Date, tempfbbchdat$TP0_FBBPET_Date), "years"))


# tempfbbch$Cingulate_Change=ifelse(is.na(tempfbbch$TP2_Cingulate), tempfbbch$TP1_Cingulate-tempfbbch$TP0_Cingulate, tempfbbch$TP2_Cingulate-tempfbbch$TP0_Cingulate)
tempfbbch$Cingulate_NTimepoints=ifelse(is.na(tempfbbch$TP2_Cingulate),2,3)
tempfbbch=tempfbbch[c(1,ncol(tempfbbch)-1,ncol(tempfbbch))]
fbbdf_lg_cwm_exp=merge(fbbdf_lg_cwm_exp,tempfbbch,by="ID")

# Now we can plot baseline vs. change, color coding by number of timepoints?

fbbdf_lg_cwm_exp$Centiloids_NTimepoints=factor(fbbdf_lg_cwm_exp$Centiloids_NTimepoints, levels=unique(fbbdf_lg_cwm_exp$Centiloids_NTimepoints))
fbbdf_lg_cwm_exp$SUVR_NTimepoints=factor(fbbdf_lg_cwm_exp$SUVR_NTimepoints, levels=unique(fbbdf_lg_cwm_exp$SUVR_NTimepoints))
fbbdf_lg_cwm_exp$Frontal_NTimepoints=factor(fbbdf_lg_cwm_exp$Frontal_NTimepoints, levels=unique(fbbdf_lg_cwm_exp$Frontal_NTimepoints))
fbbdf_lg_cwm_exp$Occipital_NTimepoints=factor(fbbdf_lg_cwm_exp$Occipital_NTimepoints, levels=unique(fbbdf_lg_cwm_exp$Occipital_NTimepoints))
fbbdf_lg_cwm_exp$Parietal_NTimepoints=factor(fbbdf_lg_cwm_exp$Parietal_NTimepoints, levels=unique(fbbdf_lg_cwm_exp$Parietal_NTimepoints))
fbbdf_lg_cwm_exp$Temporal_NTimepoints=factor(fbbdf_lg_cwm_exp$Temporal_NTimepoints, levels=unique(fbbdf_lg_cwm_exp$Temporal_NTimepoints))
fbbdf_lg_cwm_exp$Cingulate_NTimepoints=factor(fbbdf_lg_cwm_exp$Cingulate_NTimepoints, levels=unique(fbbdf_lg_cwm_exp$Cingulate_NTimepoints))

### Small module to remove cases that need to be removed from the longitudinal analysis, if any

# delids=c("")
# fbbdf_lg_cwm_exp=subset(fbbdf_lg_cwm_exp, !fbbdf_lg_cwm_exp$ID %in% delids)

###

p=ggplot(data=subset(fbbdf_lg_cwm_exp, fbbdf_lg_cwm_exp$CohortAssgn=="EOAD"), aes(x=TP0_MRIBASED_Composite_Centiloids, y=Centiloids_Change)) +
  geom_point(aes(x=TP0_MRIBASED_Composite_Centiloids, y=Centiloids_Change, col=CohortAssgn, shape=Centiloids_NTimepoints), size=2, alpha=0.7) +
  ylab("Centiloids (Ann. change)") +
  xlab("Centiloids (baseline)") +
  theme_minimal() + labs(title="FBB-PET Centiloids",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_MRIBASED_Composite_Centiloids, y = Centiloids_Change), data=fbbdf_lg_cwm_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_MRIBASED_Composite_Centiloids, y = Centiloids_Change, group=1), data=fbbdf_lg_cwm_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_amyctls_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Centiloids_BLvsCHANGE_CompWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=subset(fbbdf_lg_cwm_exp, fbbdf_lg_cwm_exp$CohortAssgn=="EOAD"), aes(x=TP0_MRIBASED_Composite_SUVR_Type2, y=SUVR_Change)) +
  geom_point(aes(x=TP0_MRIBASED_Composite_SUVR_Type2, y=SUVR_Change, col=CohortAssgn, shape=SUVR_NTimepoints), size=2, alpha=0.7) +
  ylab("FBB-PET SUVR (Ann. change)") +
  xlab("FBB-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FBB-PET Neocortical SUVR", col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_MRIBASED_Composite_SUVR_Type2, y = SUVR_Change), data=fbbdf_lg_cwm_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_MRIBASED_Composite_SUVR_Type2, y = SUVR_Change, group=1), data=fbbdf_lg_cwm_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_amysuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBBSUVR_BLvsCHANGE_CompWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

###

p=ggplot(data=subset(fbbdf_lg_cwm_exp, fbbdf_lg_cwm_exp$CohortAssgn=="EOAD"), aes(x=TP0_Frontal, y=Frontal_Change)) +
  geom_point(aes(x=TP0_Frontal, y=Frontal_Change, col=CohortAssgn, shape=Frontal_NTimepoints), size=2, alpha=0.7) +
  ylab("FBB-PET SUVR (Ann. change)") +
  xlab("FBB-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FBB-PET Frontal SUVR",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Frontal, y = Frontal_Change), data=fbbdf_lg_cwm_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Frontal, y = Frontal_Change, group=1), data=fbbdf_lg_cwm_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_amyfrontalsuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBB_Frontal_SUVR_BLvsCHANGE_CompWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=subset(fbbdf_lg_cwm_exp, fbbdf_lg_cwm_exp$CohortAssgn=="EOAD"), aes(x=TP0_Occipital, y=Occipital_Change)) +
  geom_point(aes(x=TP0_Occipital, y=Occipital_Change, col=CohortAssgn, shape=Occipital_NTimepoints), size=2, alpha=0.7) +
  ylab("FBB-PET SUVR (Ann. change)") +
  xlab("FBB-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FBB-PET Occipital SUVR",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Occipital, y = Occipital_Change), data=fbbdf_lg_cwm_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Occipital, y = Occipital_Change, group=1), data=fbbdf_lg_cwm_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_amyoccipitalsuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBB_Occipital_SUVR_BLvsCHANGE_CompWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=subset(fbbdf_lg_cwm_exp, fbbdf_lg_cwm_exp$CohortAssgn=="EOAD"), aes(x=TP0_Temporal, y=Temporal_Change)) +
  geom_point(aes(x=TP0_Temporal, y=Temporal_Change, col=CohortAssgn, shape=Temporal_NTimepoints), size=2, alpha=0.7) +
  ylab("FBB-PET SUVR (Ann. change)") +
  xlab("FBB-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FBB-PET Temporal SUVR",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none", axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Temporal, y = Temporal_Change), data=fbbdf_lg_cwm_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Temporal, y = Temporal_Change, group=1), data=fbbdf_lg_cwm_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_amytemporalsuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBB_Temporal_SUVR_BLvsCHANGE_CompWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=subset(fbbdf_lg_cwm_exp, fbbdf_lg_cwm_exp$CohortAssgn=="EOAD"), aes(x=TP0_Cingulate, y=Cingulate_Change)) +
  geom_point(aes(x=TP0_Cingulate, y=Cingulate_Change, col=CohortAssgn, shape=Cingulate_NTimepoints), size=2, alpha=0.7) +
  ylab("FBB-PET SUVR (Ann. change)") +
  xlab("FBB-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FBB-PET Cingulate SUVR", col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Cingulate, y = Cingulate_Change), data=fbbdf_lg_cwm_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Cingulate, y = Cingulate_Change, group=1), data=fbbdf_lg_cwm_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_amycingulatesuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBB_Cingulate_SUVR_BLvsCHANGE_CompWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=subset(fbbdf_lg_cwm_exp, fbbdf_lg_cwm_exp$CohortAssgn=="EOAD"), aes(x=TP0_Parietal, y=Parietal_Change)) +
  geom_point(aes(x=TP0_Parietal, y=Parietal_Change, col=CohortAssgn, shape=Parietal_NTimepoints), size=2, alpha=0.7) +
  ylab("FBB-PET SUVR (Ann. change)") +
  xlab("FBB-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FBB-PET Parietal SUVR",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none", axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Parietal, y = Parietal_Change), data=fbbdf_lg_cwm_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Parietal, y = Parietal_Change, group=1), data=fbbdf_lg_cwm_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_amyparietalsuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBB_Parietal_SUVR_BLvsCHANGE_CompWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

dbname=paste("ready_datasets/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBB_LongitudinalDatabaseCompWM_exp.xlsx",sep = "")
write.xlsx(fbbdf_lg_cwm_exp, dbname, row.names=FALSE, col.names = TRUE)

#### Liana's Request: FTP-PET Plotting baseline vs. change automatically grabbing any Timepoint available ####

# MetaROI

ftplgnames=colnames(ftpdf_ewm_lg_exp)
ftpintids=grepl("MetaROI",ftplgnames)
tempftpch=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpintids])

ftpdatids=grepl("Date",ftplgnames)
tempftpchdat=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpdatids])
tempftpch$MetaROI_Change=ifelse(is.na(tempftpch$TP2_MetaROI_MRIBASED_SUVR), 
                                (tempftpch$TP1_MetaROI_MRIBASED_SUVR-tempftpch$TP0_MetaROI_MRIBASED_SUVR)/time_length(difftime(tempftpchdat$TP1_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"), 
                                (tempftpch$TP2_MetaROI_MRIBASED_SUVR-tempftpch$TP0_MetaROI_MRIBASED_SUVR)/time_length(difftime(tempftpchdat$TP2_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"))

# tempftpch$MetaROI_Change=ifelse(is.na(tempftpch$TP2_MetaROI_MRIBASED_SUVR), tempftpch$TP1_MetaROI_MRIBASED_SUVR-tempftpch$TP0_MetaROI_MRIBASED_SUVR, tempftpch$TP2_MetaROI_MRIBASED_SUVR-tempftpch$TP0_MetaROI_MRIBASED_SUVR)
tempftpch$MetaROI_NTimepoints=ifelse(is.na(tempftpch$TP2_MetaROI_MRIBASED_SUVR),2,3)
tempftpch=tempftpch[c(1,ncol(tempftpch)-1,ncol(tempftpch))]
ftpdf_ewm_lg_exp=merge(ftpdf_ewm_lg_exp,tempftpch,by="ID")

# Frontal

ftplgnames=colnames(ftpdf_ewm_lg_exp)
ftpintids=grepl("Frontal",ftplgnames)
tempftpch=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpintids])

ftpdatids=grepl("Date",ftplgnames)
tempftpchdat=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpdatids])
tempftpch$Frontal_Change=ifelse(is.na(tempftpch$TP2_Frontal), 
                                (tempftpch$TP1_Frontal-tempftpch$TP0_Frontal)/time_length(difftime(tempftpchdat$TP1_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"), 
                                (tempftpch$TP2_Frontal-tempftpch$TP0_Frontal)/time_length(difftime(tempftpchdat$TP2_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"))



# tempftpch$Frontal_Change=ifelse(is.na(tempftpch$TP2_Frontal), tempftpch$TP1_Frontal-tempftpch$TP0_Frontal, tempftpch$TP2_Frontal-tempftpch$TP0_Frontal)
tempftpch$Frontal_NTimepoints=ifelse(is.na(tempftpch$TP2_Frontal),2,3)
tempftpch=tempftpch[c(1,ncol(tempftpch)-1,ncol(tempftpch))]
ftpdf_ewm_lg_exp=merge(ftpdf_ewm_lg_exp,tempftpch,by="ID")

# Parietal

ftplgnames=colnames(ftpdf_ewm_lg_exp)
ftpintids=grepl("Parietal",ftplgnames)
tempftpch=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpintids])

ftpdatids=grepl("Date",ftplgnames)
tempftpchdat=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpdatids])
tempftpch$Parietal_Change=ifelse(is.na(tempftpch$TP2_Parietal), 
                                 (tempftpch$TP1_Parietal-tempftpch$TP0_Parietal)/time_length(difftime(tempftpchdat$TP1_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"), 
                                 (tempftpch$TP2_Parietal-tempftpch$TP0_Parietal)/time_length(difftime(tempftpchdat$TP2_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"))


# tempftpch$Parietal_Change=ifelse(is.na(tempftpch$TP2_Parietal), tempftpch$TP1_Parietal-tempftpch$TP0_Parietal, tempftpch$TP2_Parietal-tempftpch$TP0_Parietal)
tempftpch$Parietal_NTimepoints=ifelse(is.na(tempftpch$TP2_Parietal),2,3)
tempftpch=tempftpch[c(1,ncol(tempftpch)-1,ncol(tempftpch))]
ftpdf_ewm_lg_exp=merge(ftpdf_ewm_lg_exp,tempftpch,by="ID")

# Occipital

ftplgnames=colnames(ftpdf_ewm_lg_exp)
ftpintids=grepl("Occipital",ftplgnames)
tempftpch=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpintids])

ftpdatids=grepl("Date",ftplgnames)
tempftpchdat=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpdatids])
tempftpch$Occipital_Change=ifelse(is.na(tempftpch$TP2_Occipital), 
                                  (tempftpch$TP1_Occipital-tempftpch$TP0_Occipital)/time_length(difftime(tempftpchdat$TP1_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"), 
                                  (tempftpch$TP2_Occipital-tempftpch$TP0_Occipital)/time_length(difftime(tempftpchdat$TP2_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"))


# tempftpch$Occipital_Change=ifelse(is.na(tempftpch$TP2_Occipital), tempftpch$TP1_Occipital-tempftpch$TP0_Occipital, tempftpch$TP2_Occipital-tempftpch$TP0_Occipital)
tempftpch$Occipital_NTimepoints=ifelse(is.na(tempftpch$TP2_Occipital),2,3)
tempftpch=tempftpch[c(1,ncol(tempftpch)-1,ncol(tempftpch))]
ftpdf_ewm_lg_exp=merge(ftpdf_ewm_lg_exp,tempftpch,by="ID")

# Temporal

ftplgnames=colnames(ftpdf_ewm_lg_exp)
ftpintids=grepl("Temporal",ftplgnames)
tempftpch=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpintids])

ftpdatids=grepl("Date",ftplgnames)
tempftpchdat=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpdatids])
tempftpch$Temporal_Change=ifelse(is.na(tempftpch$TP2_Temporal), 
                                 (tempftpch$TP1_Temporal-tempftpch$TP0_Temporal)/time_length(difftime(tempftpchdat$TP1_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"), 
                                 (tempftpch$TP2_Temporal-tempftpch$TP0_Temporal)/time_length(difftime(tempftpchdat$TP2_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"))


#tempftpch$Temporal_Change=ifelse(is.na(tempftpch$TP2_Temporal), tempftpch$TP1_Temporal-tempftpch$TP0_Temporal, tempftpch$TP2_Temporal-tempftpch$TP0_Temporal)
tempftpch$Temporal_NTimepoints=ifelse(is.na(tempftpch$TP2_Temporal),2,3)
tempftpch=tempftpch[c(1,ncol(tempftpch)-1,ncol(tempftpch))]
ftpdf_ewm_lg_exp=merge(ftpdf_ewm_lg_exp,tempftpch,by="ID")

# Cingulate

ftplgnames=colnames(ftpdf_ewm_lg_exp)
ftpintids=grepl("Cingulate",ftplgnames)
tempftpch=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpintids])

ftpdatids=grepl("Date",ftplgnames)
tempftpchdat=data.frame(ftpdf_ewm_lg_exp[1],ftpdf_ewm_lg_exp[,ftpdatids])
tempftpch$Cingulate_Change=ifelse(is.na(tempftpch$TP2_Cingulate), 
                                  (tempftpch$TP1_Cingulate-tempftpch$TP0_Cingulate)/time_length(difftime(tempftpchdat$TP1_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"), 
                                  (tempftpch$TP2_Cingulate-tempftpch$TP0_Cingulate)/time_length(difftime(tempftpchdat$TP2_FTPPET_Date, tempftpchdat$TP0_FTPPET_Date), "years"))


# tempftpch$Cingulate_Change=ifelse(is.na(tempftpch$TP2_Cingulate), tempftpch$TP1_Cingulate-tempftpch$TP0_Cingulate, tempftpch$TP2_Cingulate-tempftpch$TP0_Cingulate)
tempftpch$Cingulate_NTimepoints=ifelse(is.na(tempftpch$TP2_Cingulate),2,3)
tempftpch=tempftpch[c(1,ncol(tempftpch)-1,ncol(tempftpch))]
ftpdf_ewm_lg_exp=merge(ftpdf_ewm_lg_exp,tempftpch,by="ID")

###

# Now we can plot baseline vs. change, color coding by number of timepoints?

ftpdf_ewm_lg_exp$MetaROI_NTimepoints=factor(ftpdf_ewm_lg_exp$MetaROI_NTimepoints, levels=unique(ftpdf_ewm_lg_exp$MetaROI_NTimepoints))
ftpdf_ewm_lg_exp$Frontal_NTimepoints=factor(ftpdf_ewm_lg_exp$Frontal_NTimepoints, levels=unique(ftpdf_ewm_lg_exp$Frontal_NTimepoints))
ftpdf_ewm_lg_exp$Occipital_NTimepoints=factor(ftpdf_ewm_lg_exp$Occipital_NTimepoints, levels=unique(ftpdf_ewm_lg_exp$Occipital_NTimepoints))
ftpdf_ewm_lg_exp$Parietal_NTimepoints=factor(ftpdf_ewm_lg_exp$Parietal_NTimepoints, levels=unique(ftpdf_ewm_lg_exp$Parietal_NTimepoints))
ftpdf_ewm_lg_exp$Temporal_NTimepoints=factor(ftpdf_ewm_lg_exp$Temporal_NTimepoints, levels=unique(ftpdf_ewm_lg_exp$Temporal_NTimepoints))
ftpdf_ewm_lg_exp$Cingulate_NTimepoints=factor(ftpdf_ewm_lg_exp$Cingulate_NTimepoints, levels=unique(ftpdf_ewm_lg_exp$Cingulate_NTimepoints))


### Small module to remove cases that need to be removed from the longitudinal analysis, if any

# delids=c("")
# ftpdf_ewm_lg_exp=subset(ftpdf_ewm_lg_exp, !ftpdf_ewm_lg_exp$ID %in% delids)

p=ggplot(data=ftpdf_ewm_lg_exp, aes(x=TP0_MetaROI_MRIBASED_SUVR, y=MetaROI_Change)) +
  geom_point(aes(x=TP0_MetaROI_MRIBASED_SUVR, y=MetaROI_Change, col=CohortAssgn, shape=MetaROI_NTimepoints), size=2, alpha=0.7) +
  ylab("FTP-PET SUVR (Ann. change)") +
  xlab("FTP-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FTP-PET metaROI SUVR",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_MetaROI_MRIBASED_SUVR, y = MetaROI_Change), data=ftpdf_ewm_lg_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_MetaROI_MRIBASED_SUVR, y = MetaROI_Change, group=1), data=ftpdf_ewm_lg_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_ftpmeta_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FTPmetaROI_BLvsCHANGE_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=ftpdf_ewm_lg_exp, aes(x=TP0_Frontal, y=Frontal_Change)) +
  geom_point(aes(x=TP0_Frontal, y=Frontal_Change, col=CohortAssgn, shape=Frontal_NTimepoints), size=2, alpha=0.7) +
  ylab("FTP-PET SUVR (Ann. change)") +
  xlab("FTP-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FTP-PET Frontal SUVR", col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Frontal, y = Frontal_Change), data=ftpdf_ewm_lg_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Frontal, y = Frontal_Change, group=1), data=ftpdf_ewm_lg_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_ftpfrontalsuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FTP_Frontal_SUVR_BLvsCHANGE_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=ftpdf_ewm_lg_exp, aes(x=TP0_Occipital, y=Occipital_Change)) +
  geom_point(aes(x=TP0_Occipital, y=Occipital_Change, col=CohortAssgn, shape=Occipital_NTimepoints), size=2, alpha=0.7) +
  ylab("FTP-PET SUVR (Ann. change)") +
  xlab("FTP-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FTP-PET Occipital SUVR", col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Occipital, y = Occipital_Change), data=ftpdf_ewm_lg_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Occipital, y = Occipital_Change, group=1), data=ftpdf_ewm_lg_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_ftpoccipitalsuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FTP_Occipital_SUVR_BLvsCHANGE_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=ftpdf_ewm_lg_exp, aes(x=TP0_Temporal, y=Temporal_Change)) +
  geom_point(aes(x=TP0_Temporal, y=Temporal_Change, col=CohortAssgn, shape=Temporal_NTimepoints), size=2, alpha=0.7) +
  ylab("FTP-PET SUVR (Ann. change)") +
  xlab("FTP-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FTP-PET Temporal SUVR",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Temporal, y = Temporal_Change), data=ftpdf_ewm_lg_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Temporal, y = Temporal_Change, group=1), data=ftpdf_ewm_lg_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_ftptemporalsuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FTP_Temporal_SUVR_BLvsCHANGE_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=ftpdf_ewm_lg_exp, aes(x=TP0_Cingulate, y=Cingulate_Change)) +
  geom_point(aes(x=TP0_Cingulate, y=Cingulate_Change, col=CohortAssgn, shape=Cingulate_NTimepoints), size=2, alpha=0.7) +
  ylab("FTP-PET SUVR (Ann. change)") +
  xlab("FTP-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FTP-PET Cingulate SUVR",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Cingulate, y = Cingulate_Change), data=ftpdf_ewm_lg_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Cingulate, y = Cingulate_Change, group=1), data=ftpdf_ewm_lg_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_ftpcingulatesuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FTP_Cingulate_SUVR_BLvsCHANGE_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=ftpdf_ewm_lg_exp, aes(x=TP0_Parietal, y=Parietal_Change)) +
  geom_point(aes(x=TP0_Parietal, y=Parietal_Change, col=CohortAssgn, shape=Parietal_NTimepoints), size=2, alpha=0.7) +
  ylab("FTP-PET SUVR (Ann. change)") +
  xlab("FTP-PET SUVR (baseline)") +
  theme_minimal() + labs(title="FTP-PET Parietal SUVR",col="Cohort", shape="N Timepoints")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = TP0_Parietal, y = Parietal_Change), data=ftpdf_ewm_lg_exp, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = TP0_Parietal, y = Parietal_Change, group=1), data=ftpdf_ewm_lg_exp, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_ftpparietalsuvr_change=p +geom_hline(yintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FTP_Parietal_SUVR_BLvsCHANGE_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

dbname=paste("ready_datasets/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FTP_LongitudinalDatabaseErodedWM_exp.xlsx",sep = "")
write.xlsx(ftpdf_ewm_lg_exp, dbname, row.names=FALSE, col.names = TRUE)

#### Plots Changes in FBB vs. Changes in FTP ####

df_fbblg_changes=data.frame(fbbdf_lg_cwm_exp$ID,
                            fbbdf_lg_cwm_exp$CohortAssgn,
                            fbbdf_lg_cwm_exp$Centiloids_Change,
                            fbbdf_lg_cwm_exp$Frontal_Change,
                            fbbdf_lg_cwm_exp$Occipital_Change,
                            fbbdf_lg_cwm_exp$Parietal_Change,
                            fbbdf_lg_cwm_exp$Temporal_Change,
                            fbbdf_lg_cwm_exp$Cingulate_Change)

colnames(df_fbblg_changes)=c("ID",
                             "Cohort",
                             "FBB_Centiloids_Change",
                             "FBB_Frontal_Change",
                             "FBB_Occipital_Change",
                             "FBB_Parietal_Change",
                             "FBB_Temporal_Change",
                             "FBB_Cingulate_Change")

df_ftplg_changes=data.frame(ftpdf_ewm_lg_exp$ID,
                            ftpdf_ewm_lg_exp$CohortAssgn,
                            ftpdf_ewm_lg_exp$MetaROI_Change,
                            ftpdf_ewm_lg_exp$Frontal_Change,
                            ftpdf_ewm_lg_exp$Occipital_Change,
                            ftpdf_ewm_lg_exp$Parietal_Change,
                            ftpdf_ewm_lg_exp$Temporal_Change,
                            ftpdf_ewm_lg_exp$Cingulate_Change)

colnames(df_ftplg_changes)=c("ID",
                             "Cohort",
                             "FTP_MetaROI_Change",
                             "FTP_Frontal_Change",
                             "FTP_Occipital_Change",
                             "FTP_Parietal_Change",
                             "FTP_Temporal_Change",
                             "FTP_Cingulate_Change")

df_lg_changes=merge(df_fbblg_changes, df_ftplg_changes, by=c("ID","Cohort"))

p=ggplot(data=df_lg_changes, aes(x=FBB_Centiloids_Change, y=FTP_MetaROI_Change)) +
  geom_point(aes(x=FBB_Centiloids_Change, y=FTP_MetaROI_Change, col=Cohort), size=2, alpha=0.7) +
  ylab("\u0394FTP-PET metaSUVR") +
  xlab("\u0394FBB-PET Centiloids") +
  theme_minimal() + labs(title="Composite vs. MetaROI",col="Cohort")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = FBB_Centiloids_Change, y = FTP_MetaROI_Change), data=df_lg_changes, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = FBB_Centiloids_Change, y = FTP_MetaROI_Change, group=1), data=df_lg_changes, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_centmeta_change=p +geom_hline(yintercept=0,col="black") +geom_vline(xintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_CompMeta_CHANGES_CompWM_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=df_lg_changes, aes(x=FBB_Frontal_Change, y=FTP_Frontal_Change)) +
  geom_point(aes(x=FBB_Frontal_Change, y=FTP_Frontal_Change, col=Cohort), size=2, alpha=0.7) +
  ylab("\u0394FTP-PET SUVR") +
  xlab("\u0394FBB-PET SUVR") +
  theme_minimal() + labs(title="Frontal",col="Cohort")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = FBB_Frontal_Change, y = FTP_Frontal_Change), data=df_lg_changes, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = FBB_Frontal_Change, y = FTP_Frontal_Change, group=1), data=df_lg_changes, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_frontal_change=p +geom_hline(yintercept=0,col="black") +geom_vline(xintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Frontal_CHANGES_CompWM_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=df_lg_changes, aes(x=FBB_Occipital_Change, y=FTP_Occipital_Change)) +
  geom_point(aes(x=FBB_Occipital_Change, y=FTP_Occipital_Change, col=Cohort), size=2, alpha=0.7) +
  ylab("\u0394FTP-PET SUVR") +
  xlab("\u0394FBB-PET SUVR") +
  theme_minimal() + labs(title="Occipital",col="Cohort")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = FBB_Occipital_Change, y = FTP_Occipital_Change), data=df_lg_changes, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = FBB_Occipital_Change, y = FTP_Occipital_Change, group=1), data=df_lg_changes, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_occipital_change=p +geom_hline(yintercept=0,col="black") +geom_vline(xintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Occipital_CHANGES_CompWM_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=df_lg_changes, aes(x=FBB_Parietal_Change, y=FTP_Parietal_Change)) +
  geom_point(aes(x=FBB_Parietal_Change, y=FTP_Parietal_Change, col=Cohort), size=2, alpha=0.7) +
  ylab("\u0394FTP-PET SUVR") +
  xlab("\u0394FBB-PET SUVR") +
  theme_minimal() + labs(title="Parietal",col="Cohort")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = FBB_Parietal_Change, y = FTP_Parietal_Change), data=df_lg_changes, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = FBB_Parietal_Change, y = FTP_Parietal_Change, group=1), data=df_lg_changes, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_parietal_change=p +geom_hline(yintercept=0,col="black") +geom_vline(xintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Parietal_CHANGES_CompWM_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=df_lg_changes, aes(x=FBB_Temporal_Change, y=FTP_Temporal_Change)) +
  geom_point(aes(x=FBB_Temporal_Change, y=FTP_Temporal_Change, col=Cohort), size=2, alpha=0.7) +
  ylab("\u0394FTP-PET SUVR") +
  xlab("\u0394FBB-PET SUVR") +
  theme_minimal() + labs(title="Temporal",col="Cohort")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = FBB_Temporal_Change, y = FTP_Temporal_Change), data=df_lg_changes, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = FBB_Temporal_Change, y = FTP_Temporal_Change, group=1), data=df_lg_changes, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_temporal_change=p +geom_hline(yintercept=0,col="black") +geom_vline(xintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Temporal_CHANGES_CompWM_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

p=ggplot(data=df_lg_changes, aes(x=FBB_Cingulate_Change, y=FTP_Cingulate_Change)) +
  geom_point(aes(x=FBB_Cingulate_Change, y=FTP_Cingulate_Change, col=Cohort), size=2, alpha=0.7) +
  ylab("\u0394FTP-PET SUVR") +
  xlab("\u0394FBB-PET SUVR") +
  theme_minimal() + labs(title="Cingulate",col="Cohort")  +
  theme(legend.position = "none",axis.title.x=element_text(size=10), axis.text.x = element_text(size=8), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = FBB_Cingulate_Change, y = FTP_Cingulate_Change), data=df_lg_changes, method = "lm", se = FALSE, cex=0.7, col="red", lty="dotted") + geom_ribbon(aes(x = FBB_Cingulate_Change, y = FTP_Cingulate_Change, group=1), data=df_lg_changes, stat = "smooth", method = "lm", alpha = .15) +
  scale_color_manual(values=c("#052049","#007CBE","#B4B9BF"))

p_cingulate_change=p +geom_hline(yintercept=0,col="black") +geom_vline(xintercept=0,col="black")

dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Cingulate_CHANGES_CompWM_ErodedWM.jpg",sep = "")
ggsave(dbname, plot=p, width = 15, height = 15, units = "cm")

dbname=paste("ready_datasets/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_LongitudinalDatabases_Change_Merged.xlsx",sep = "")
write.xlsx(df_lg_changes, dbname, row.names=FALSE, col.names = TRUE)

### Save workspace again ###

rdatafname=paste("logs/Datasets_Generation_LONI_Step3_",format(Sys.time(), "%Y-%m-%d"),".Rdata",sep = "")
save.image(file=rdatafname)

### Quickplot to show comparison between CBL based and CompWM based reference for global amyloid ###

quickplot=fbbdf[c(1,12,20)]
quickplot$GlobalSUVR_Compwm=quickplot$MRIBASED_Composite_SUVR_Type2/quickplot$CompositeWM_MRIBASED_SUVR

p=ggplot(data=quickplot, aes(x=MRIBASED_Composite_SUVR_Type2, y=GlobalSUVR_Compwm)) +
  geom_point(aes(x=MRIBASED_Composite_SUVR_Type2, y=GlobalSUVR_Compwm), size=2.5, alpha=0.7) +
  ylab("Neocortical MRI-based FBB-PET SUVR (CompWM reference)") +
  xlab("Neocortical MRI-based FBB-PET SUVR (WCBL reference)") +
  theme_minimal() + labs(col="CohortAssgn", fill="CohortAssgn")  +
  theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), axis.text.y=element_text(size=10, face="bold"), axis.title.y=element_text(size=10)) +
  geom_smooth(aes(x = MRIBASED_Composite_SUVR_Type2, y = GlobalSUVR_Compwm), data=quickplot, method = "loess", se = TRUE, cex=1, col="red", lty="dashed") 
pname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_FBBPET_GlobalAmy_WcblvsCompWM.jpg",sep = "")
ggsave(pname,plot=p, height = 14, width=12, units = "cm")
  
  #### Service for dashboard ####

  npt_screened=nrow(fbb_pt_reads)
  npt_discordant=length(which(fbb_pt_reads$Screening_PETONLY_Disagreement=="1"))
  rate_concordance=length(which(fbb_pt_reads$Screening_PETONLY_Disagreement=="0"))/nrow(fbb_pt_reads)*100
  rate_amypos=length(which(fbb_pt_reads$Screening_PETONLY_Final_Read=="1"))/nrow(fbb_pt_reads)*100
  rate_amyneg=100-rate_amypos
  neoad=length(which(fbb_pt_reads$Cohort=="EOAD"))
  neononad=nrow(fbb_pt_reads)-length(which(fbb_pt_reads$Cohort=="EOAD"))

  ftpdf_bl_eoad=subset(ftpdf_bl, ftpdf_bl$CohortAssgn=="EOAD")
  ftpdf_bl_eononad=subset(ftpdf_bl, ftpdf_bl$CohortAssgn=="EOnonAD")
  rate_tauposeoad=length(which(ftpdf_bl_eoad$MetaROI_MRIBASED_SUVR>1.2))/nrow(ftpdf_bl_eoad)*100
  rate_tauposeononad=length(which(ftpdf_bl_eononad$MetaROI_MRIBASED_SUVR>1.2))/nrow(ftpdf_bl_eononad)*100
  nlongftp=nrow(ftpdf_ewm_lg_exp)
  nlongfbb=nrow(fbbdf_lg_cwm_exp)
  nlongfbbftp=length(which(fbbdf_lg_cwm_exp$ID %in% ftpdf_ewm_lg_exp$ID)==TRUE)
  nptfbbftp=nrow(subset(fbbftpdf_bl, !fbbftpdf_bl$CohortAssgn.x=="CN"))
  nfbbftp=nrow(fbbftpdf_bl)
  ntotalimagesproc=nrow(fbbdf)+nrow(ftpdf)+nrow(fdgdf)

  dbname=paste("reports/",format(Sys.time(), "%Y-%m-%d"),"_","LEADS_Dashboard_data.Rdata",sep = "")

  save(p_ctlsvsvisread_poscp,
       p_ctlsbycohort_poscp,
       p_qcpetonlyvsmribasedamy,
       p_metaroiamy_all,
       p_metaroiamyctls_all,
       p_ftppet_metabraak,
       p_ftpmeta_change,
       p_amyctls_change,
       p_amysuvr_change,
       p_amyfrontalsuvr_change,
       p_amyoccipitalsuvr_change,
       p_amycingulatesuvr_change,
       p_amyparietalsuvr_change,
       p_amytemporalsuvr_change,
       p_ftpfrontalsuvr_change,
       p_ftpoccipitalsuvr_change,
       p_ftpcingulatesuvr_change,
       p_ftpparietalsuvr_change,
       p_ftptemporalsuvr_change,
       p_centmeta_change,
       p_frontal_change,
       p_occipital_change,
       p_temporal_change,
       p_parietal_change,
       p_cingulate_change,
       p_fbb_cn,
       p_fbb_eoad,
       p_fbb_eononad,
       p_ftp_cn,
       p_ftp_eoad,
       p_ftp_eononad,
       p_timetrend_fbb,
       p_timetrend_ftp,
       p_timetrend_fdg,
       npt_screened,
       npt_discordant,
       rate_amypos,
       rate_amyneg,
       neoad,
       neononad,
       rate_tauposeoad,
       rate_tauposeononad,
       nlongfbb,
       nlongftp,
       nlongfbbftp,
       rate_concordance,
       nptfbbftp,
       ntotalimagesproc,
       nfbbftp, 
       ftplg_plot2,
       ftplg_plot2_apoe,
       ftplg_plot2_cdr,
       fbblg_plot2,
       fbblg_plot2_apoe,
       fbblg_plot2_cdr,
       file = dbname)

  #### Final Greetings ####

  message("**********************")
  message(paste("Finished! I saved Plots and Tables for you to enjoy. \n********************** \nHere is a piece a wisdom:"))
  print(statquote(topic="Science"))
  message("**********************")


} else {

  message(paste("Something is wrong in the cohort column.\nIt is possible one or more symptomatic patients have been processed but not read yet. \nHere is subjects that do not match:\n"))
  print(idsedc_sbs)
  fnameqc=paste("logs/FailedCohortQC_",format(Sys.time(), "%Y-%m-%d"),".csv",sep = "")
  write.csv(idsedc_sbs,fnameqc, row.names = F, na="")
}


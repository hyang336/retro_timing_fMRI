library(rio)
library(tidyverse)

#------------------------------------------------------------------------------------------------------------
# Load files and extract subject IDs
SAMS1_dir='C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\2023-08-21\\SAMS1\\'
SAMS2_dir='C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\2023-08-21\\SAMS2\\'
AMASS_haams_dir='C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\2023-08-21\\AMASS\\haams_onset_files\\'
AMASS_yaams_dir='C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\2023-08-21\\AMASS\\yaams_onset_files\\'

SSID_pattern='sub-\\d{3,4}'

SAMS1_files=list.files(SAMS1_dir)
SAMS1_SSID=str_extract(SAMS1_files,SSID_pattern)

SAMS2_files=list.files(SAMS2_dir)
SAMS2_SSID=str_extract(SAMS2_files,SSID_pattern)

AMASS_haams_files=list.files(AMASS_haams_dir)
AMASS_haams_SSID=str_extract(AMASS_haams_files,SSID_pattern)

AMASS_yaams_files=list.files(AMASS_yaams_dir)
AMASS_yaams_SSID=str_extract(AMASS_yaams_files,SSID_pattern)

#----------------------------------------------------------------------------------------------------------
#load in each subject's data for each run, and calculate trial length based on the time stamps, need to handle different projects differently due to file 

###!!! The duration of each trial component (pregoal, goalcue, prestim, stim) are pulled out from the onset files !!!###
###!!! While the trial-to-trial ITI was calculated using the differences of pregoal_onset between successive trials !!!###

AMASS_frame=data.frame(SSID=character(),trial=numeric(),run=numeric(),PregoalITI=numeric(),GoalCue_Pregoal=numeric(),Prestim_GoalCue=numeric(),Stim_Prestim=numeric(),Stim_duration=numeric())

## AMASS_haams
for (i in c(1:length(AMASS_haams_files))){
  ss_frame=import(paste(AMASS_haams_dir,AMASS_haams_files[i],sep=''))
  runs=unique(ss_frame$run)
  for (j in c(1:length(runs))){
    run_frame<-ss_frame %>% filter(run==runs[j]) %>% mutate(PregoalITI=pregoal_onset - lag(pregoal_onset, default = NA))
    temp_frame=data.frame(SSID=character(dim(run_frame)[1]),trial=numeric(dim(run_frame)[1]),run=numeric(dim(run_frame)[1]),PregoalITI=numeric(dim(run_frame)[1]),GoalCue_Pregoal=numeric(dim(run_frame)[1]),Prestim_GoalCue=numeric(dim(run_frame)[1]),Stim_Prestim=numeric(dim(run_frame)[1]),Stim_duration=numeric(dim(run_frame)[1]))
    temp_frame$trial=run_frame$trial_index
    temp_frame$run=runs[j]
    temp_frame$PregoalITI=run_frame$PregoalITI
    temp_frame$GoalCue_Pregoal=run_frame$pregoal_duration
    temp_frame$Prestim_GoalCue=run_frame$goal_duration
    temp_frame$Stim_Prestim=run_frame$prestim_duration
    temp_frame$Stim_duration=run_frame$stim_duration
    temp_frame$SSID=AMASS_haams_SSID[i]
    AMASS_frame=rbind(AMASS_frame,temp_frame)
  }
  
}

## AMASS_yaams
for (i in c(1:length(AMASS_yaams_files))){
  ss_frame=import(paste(AMASS_yaams_dir,AMASS_yaams_files[i],sep=''))
  runs=unique(ss_frame$run)
  for (j in c(1:length(runs))){
    run_frame<-ss_frame %>% filter(run==runs[j]) %>% mutate(PregoalITI=pregoal_onset - lag(pregoal_onset, default = NA))
    temp_frame=data.frame(SSID=character(dim(run_frame)[1]),trial=numeric(dim(run_frame)[1]),run=numeric(dim(run_frame)[1]),PregoalITI=numeric(dim(run_frame)[1]),GoalCue_Pregoal=numeric(dim(run_frame)[1]),Prestim_GoalCue=numeric(dim(run_frame)[1]),Stim_Prestim=numeric(dim(run_frame)[1]),Stim_duration=numeric(dim(run_frame)[1]))
    temp_frame$trial=run_frame$trial_index
    temp_frame$run=runs[j]
    temp_frame$PregoalITI=run_frame$PregoalITI
    temp_frame$GoalCue_Pregoal=run_frame$pregoal_duration
    temp_frame$Prestim_GoalCue=run_frame$goal_duration
    temp_frame$Stim_Prestim=run_frame$prestim_duration
    temp_frame$Stim_duration=run_frame$stim_duration
    temp_frame$SSID=AMASS_yaams_SSID[i]
    AMASS_frame=rbind(AMASS_frame,temp_frame)
  }
  
}


## SAMS1
SAMS1_frame=data.frame(SSID=character(),trial=numeric(),run=numeric(),phase=character(),Trial_duration=numeric())
for (i in c(1:length(SAMS1_files))){
  ss_frame_study=import(paste(SAMS1_dir,SAMS1_files[i],'/',SAMS1_files[i],'_task-AssocMemStudy_onsets.csv',sep=''))
  ss_frame_test=import(paste(SAMS1_dir,SAMS1_files[i],'/',SAMS1_files[i],'_task-AssocMemTest_onsets.csv',sep=''))
  study_runs=unique(ss_frame_study$run_id)
  test_runs=unique(ss_frame_test$run_id)
  #study
  for (j in c(1:length(study_runs))){
    run_frame<-ss_frame_study %>% filter(run_id==study_runs[j]) %>% mutate(Trial_duration=onset - lag(onset, default = NA))
    temp_frame=data.frame(SSID=character(dim(run_frame)[1]),trial=numeric(dim(run_frame)[1]),run=numeric(dim(run_frame)[1]),Trial_duration=numeric(dim(run_frame)[1]))
    temp_frame$trial=run_frame$trial_id
    temp_frame$run=study_runs[j]
    temp_frame$Trial_duration=run_frame$Trial_duration
    temp_frame$SSID=SAMS1_SSID[i]
    temp_frame$phase='study'
    SAMS1_frame=rbind(SAMS1_frame,temp_frame)
  }
  #test
  for (j in c(1:length(test_runs))){
    run_frame<-ss_frame_test %>% filter(run_id==test_runs[j]) %>% mutate(Trial_duration=onset - lag(onset, default = NA))
    temp_frame=data.frame(SSID=character(dim(run_frame)[1]),trial=numeric(dim(run_frame)[1]),run=numeric(dim(run_frame)[1]),Trial_duration=numeric(dim(run_frame)[1]))
    temp_frame$trial=run_frame$trial_id
    temp_frame$run=test_runs[j]
    temp_frame$Trial_duration=run_frame$Trial_duration
    temp_frame$SSID=SAMS1_SSID[i]
    temp_frame$phase='test'
    SAMS1_frame=rbind(SAMS1_frame,temp_frame)
  }
  
}


## SAMS2, for some reason the file structure is inconsistent among subjects
SAMS2_frame=data.frame(SSID=character(),run=numeric(),phase=character(),Trial_duration=numeric())
for (i in c(1:length(SAMS2_files))){
  #handle files with (1) in their name
  ss_folder=paste(SAMS2_dir,SAMS2_files[i],'/behav/',sep='')
  
  ss_files=list.files(ss_folder)
  study_run_X <- numeric(0)
  test_run_X <- numeric(0)
  for (file_path in ss_files) {
    file_name <- basename(file_path)
    if (grepl("study_run", file_name)) {
      study_run_X <- c(study_run_X, as.numeric(sub(".*study_run_(\\d+)_.*", "\\1", file_name)))
    } else if (grepl("test_run", file_name)) {
      test_run_X <- c(test_run_X, as.numeric(sub(".*test_run_(\\d+)_.*", "\\1", file_name)))
    }
  }
  
  study_runs=unique(study_run_X)
  test_runs=unique(test_run_X)

  #study
  for (j in c(1:length(study_runs))){
    ss_frame_study=import(paste(SAMS2_dir,SAMS2_files[i],'/behav/',SAMS2_SSID[i],'_study_run_',j,'_data.csv',sep=''))
    if (dim(ss_frame_study)[1]>0){
      run_frame<-ss_frame_study %>% mutate(Trial_duration=onset - lag(onset, default = NA))
      temp_frame=data.frame(SSID=character(dim(run_frame)[1]),run=numeric(dim(run_frame)[1]),Trial_duration=numeric(dim(run_frame)[1]))

      temp_frame$run=study_runs[j]
      temp_frame$Trial_duration=run_frame$Trial_duration
      temp_frame$SSID=SAMS2_SSID[i]
      temp_frame$phase='study'
      SAMS2_frame=rbind(SAMS2_frame,temp_frame)
    }
  }
  #test
  for (j in c(1:length(test_runs))){
    ss_frame_test=import(paste(SAMS2_dir,SAMS2_files[i],'/behav/',SAMS2_SSID[i],'_test_run_',j,'_data.csv',sep=''))
    if (dim(ss_frame_test)[1]>0){
      run_frame<-ss_frame_test %>% mutate(Trial_duration=onset - lag(onset, default = NA))
      temp_frame=data.frame(SSID=character(dim(run_frame)[1]),run=numeric(dim(run_frame)[1]),Trial_duration=numeric(dim(run_frame)[1]))

      temp_frame$run=test_runs[j]
      temp_frame$Trial_duration=run_frame$Trial_duration
      temp_frame$SSID=SAMS2_SSID[i]
      temp_frame$phase='test'
      SAMS2_frame=rbind(SAMS2_frame,temp_frame)
    }
  }
  
}

#----------------------------------------------------------------------------------------------------------
#histogram

# AMASS
AMASS_weird=AMASS_frame[AMASS_frame$PregoalITI>20.16&!is.na(AMASS_frame$PregoalITI),]
  
AMASS_hist_ITI=ggplot(data=AMASS_frame,aes(x=PregoalITI))+geom_histogram(binwidth = 0.001,aes(y=..count../sum(..count..)))+labs(title='trial length distribution (all subjects)', x='trial length (pregoal onset to pregoal onset)', y='trial proportion')+scale_x_continuous(breaks=seq(20,20.2,by=0.01))

AMASS_hist_pregoal_duration=ggplot(data=AMASS_frame,aes(x=GoalCue_Pregoal))+geom_histogram(binwidth = 0.001,aes(y=..count../sum(..count..)))+labs(title='pregoal duration distribution (all subjects)', x='pregoal duration', y='trial proportion')+scale_x_continuous(breaks=seq(8,8.2,by=0.01))

AMASS_hist_goal_duration=ggplot(data=AMASS_frame,aes(x=Prestim_GoalCue))+geom_histogram(binwidth = 0.001,aes(y=..count../sum(..count..)))+labs(title='goal-cue duration distribution (all subjects)', x='goal-cue duration', y='trial proportion')+scale_x_continuous(breaks=seq(2,2.2,by=0.001))

AMASS_hist_prestim_duration=ggplot(data=AMASS_frame,aes(x=Stim_Prestim))+geom_histogram(binwidth = 0.001,aes(y=..count../sum(..count..)))+labs(title='prestim duration distribution (all subjects)', x='prestim duration', y='trial proportion')+scale_x_continuous(breaks=seq(6,6.2,by=0.01))

AMASS_hist_stim_duration=ggplot(data=AMASS_frame,aes(x=Stim_duration))+geom_histogram(binwidth = 0.001,aes(y=..count../sum(..count..)))+labs(title='stim duration distribution (all subjects)', x='stim duration', y='trial proportion')+scale_x_continuous(breaks=seq(4,4.2,by=0.01))

AMASS_hist_ITI
AMASS_hist_pregoal_duration
AMASS_hist_goal_duration
AMASS_hist_prestim_duration
AMASS_hist_stim_duration

# SAMS1
SAMS1_weird1=SAMS1_frame[SAMS1_frame$Trial_duration<11.99&!is.na(SAMS1_frame$Trial_duration),]
SAMS1_weird2=SAMS1_frame[SAMS1_frame$Trial_duration>12.1&!is.na(SAMS1_frame$Trial_duration),]

SAMS1_hist_ITI=ggplot(data=SAMS1_frame,aes(x=Trial_duration))+geom_histogram(binwidth = 0.001,aes(y=..count../sum(..count..)))+labs(title='trial length distribution (all subjects)', x='trial length (onset to onset)', y='trial proportion')+xlim(11.9,12.1)
SAMS1_hist_ITI


# SAMS2
SAMS2_weird1=SAMS2_frame[SAMS2_frame$Trial_duration<12&!is.na(SAMS2_frame$Trial_duration),]

SAMS2_hist_ITI=ggplot(data=SAMS2_frame,aes(x=Trial_duration))+geom_histogram(binwidth = 0.001,aes(y=..count../sum(..count..)))+labs(title='trial length distribution (all subjects)', x='trial length (onset to onset)', y='trial proportion')+xlim(12,12.5)
SAMS2_hist_ITI


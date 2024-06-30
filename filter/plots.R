library(ggplot2)
library(tidyverse)
library(MASS)

data010Hz=read_csv('C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\Naive_filter\\V1_phase_results.csv')
output_dir='C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\Naive_filter'

#initialize dataframes
col_names <- c("SSID", "V1_BOLD_corr", "V1_phase_corr")
cor_df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(cor_df) <- col_names

col_names <- c("TR", "run", "BOLD", "SSID")
BOLD_df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(BOLD_df) <- col_names

col_names <- c("TR", "run", "phase", "SSID")
phase_df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(phase_df) <- col_names

# loop over subject
subs=list.dirs('C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\Naive_filter',full.names=TRUE,recursive=FALSE)

for (i in seq_along(subs)){
  BOLD_files_V1=list.files(subs[i],pattern='V1_mean_BOLD_signal_bandpass.txt',full.names = TRUE) 
  phase_files_V1=list.files(subs[i],pattern='V1_phase_bandpass.txt',full.names = TRUE) 
  #load BOLD V1 data
  BOLD_v1_list=list()
  phase_v1_list=list()
  
  sub_match=gregexpr(pattern='sub-\\d{3,4}',subs[i],perl = TRUE)
  sub_extracted <- regmatches(subs[i], sub_match)
  sub_extracted=unlist(sub_extracted)
  
  for (j in seq_along(BOLD_files_V1)){
    run_match=gregexpr(pattern='run-\\d{2}',BOLD_files_V1[j],perl = TRUE)
    run_extracted <- regmatches(BOLD_files_V1[j], run_match)
    run_extracted=unlist(run_extracted)
    temp=read_csv(BOLD_files_V1[j],col_names = FALSE)
    if (length(unlist(temp))<215){
      temp<-temp%>%
        bind_rows(tibble(X1=rep(NA,215-length(unlist(temp)))))
    }
    temp=temp%>%
      rename(!!run_extracted:=X1)
    BOLD_v1_list[[j]]=temp
  }
  BOLD_V1_frame=as.data.frame(BOLD_v1_list)
  #remove runs that are not full length
  BOLD_V1_frame <- BOLD_V1_frame %>%
    select(where(~ !any(is.na(.))))
  # compute pairwise correlation of filtered V1 signal and phase across runs
  BOLD_V1_corr=cor(BOLD_V1_frame)
  BOLD_V1_corr_vector <- BOLD_V1_corr[upper.tri(BOLD_V1_corr)]
  
  #load phase V1 data
  phase_v1_list=list()
  for (j in seq_along(phase_files_V1)){
    run_match=gregexpr(pattern='run-\\d{2}',phase_files_V1[j],perl = TRUE)
    run_extracted <- regmatches(phase_files_V1[j], run_match)
    run_extracted=unlist(run_extracted)
    temp=read_csv(phase_files_V1[j],col_names = FALSE)
    if (length(unlist(temp))<215){
      temp<-temp%>%
        bind_rows(tibble(X1=rep(NA,215-length(unlist(temp)))))
    }
    temp=temp%>%
      rename(!!run_extracted:=X1)
    phase_v1_list[[j]]=temp
  }
  phase_V1_frame=as.data.frame(phase_v1_list)
  #remove runs that are not full length
  phase_V1_frame <- phase_V1_frame %>%
    select(where(~ !any(is.na(.))))
  # compute pairwise correlation of filtered V1 signal and phase across runs
  phase_V1_corr=cor(phase_V1_frame)
  phase_V1_corr_vector <- phase_V1_corr[upper.tri(phase_V1_corr)]
  # save to dataframe
    #save mean BOLD
  BOLD_V1_frame<-BOLD_V1_frame%>%
    mutate(TR=row_number())
  BOLD_V1_frame_long<-BOLD_V1_frame%>%
    pivot_longer(cols=-TR,
                 names_to = 'run',
                 values_to = 'BOLD')
  BOLD_V1_frame_long<-BOLD_V1_frame_long%>%
    mutate(SSID=sub_extracted)
  BOLD_df=rbind(BOLD_df,BOLD_V1_frame_long)
    #save phase
  phase_V1_frame<-phase_V1_frame%>%
    mutate(TR=row_number())
  phase_V1_frame_long<-phase_V1_frame%>%
    pivot_longer(cols=-TR,
                 names_to = 'run',
                 values_to = 'phase')
  phase_V1_frame_long<-phase_V1_frame_long%>%
    mutate(SSID=sub_extracted)
  phase_df=rbind(phase_df,phase_V1_frame_long)
    #save correlations
  BOLD_V1_corr_frame=as.data.frame(BOLD_V1_corr_vector)
  colnames(BOLD_V1_corr_frame)='V1_BOLD_corr'
  phase_V1_corr_frame=as.data.frame(phase_V1_corr_vector)
  colnames(phase_V1_corr_frame)='V1_phase_corr'
  corr_frame=cbind(BOLD_V1_corr_frame,phase_V1_corr_frame)
  corr_frame<-corr_frame%>%
    mutate(SSID=sub_extracted)
  cor_df=rbind(cor_df,corr_frame)
}

#save dataframes to CSV
write.csv(BOLD_df,file=paste0(output_dir,'\\BOLD_df.csv'))
write.csv(phase_df,file=paste0(output_dir,'\\phase_df.csv'))
###############################plots#############################################
dir.create(paste0(output_dir,'\\histograms'))
# Ripley's L function to quantify spread
Ripley.K <- function(x, scale) {
  # Arguments:
  # x is an array of data.
  # scale (not actually used) is an option to rescale the data.
  #
  # Return value:
  # A function that calculates Ripley's K for any value between 0 and 1 (or `scale`).
  #
  x.pairs <- outer(x, x, function(a,b) abs(a-b))  # All pairwise distances
  x.pairs <- x.pairs[lower.tri(x.pairs)]          # Distances between distinct pairs
  if(missing(scale)) scale <- diff(range(x.pairs))# Rescale distances to [0,1]
  x.pairs <- x.pairs / scale
  #
  # The built-in `ecdf` function returns the proportion of values in `x.pairs` that
  # are less than or equal to its argument.
  #
  return (ecdf(x.pairs))
}
#
# The one-dimensional L function.
# It merely subtracts 1 - (1-y)^2 from `Ripley.K(x)(y)`.  
# Its argument `x` is an array of data values.
#
Ripley.L <- function(x) {function(y) Ripley.K(x)(y) - 1 + (1-y)^2}
# the location of the first peak of Ripley's L function denotes the width of the tightest peak and its concentration


TRs=unique(phase_df$TR)
for (i in TRs){
  png(paste0(output_dir,'\\histograms\\phase_hist_',i,'.png'))
  hist(phase_df$phase[phase_df$TR==i],nclass=100)
  dev.off()
}# The middle TRs have pretty spread-out distributions of phase, which is expected since subjects likely diverged in their behaviors (e.g., eye-closing pattern) during a run.


################quantify how much error we will make if we estimate a subject's run-onsets with the initial phase of all other subjects#######
for (sub in unique(phase_df$SSID)){
  #leaving the current sub out
  phase_df_remain=phase_df[phase_df$SSID!=sub,]
  
  #calculate the peaks of phase histogram (using density()) in the remaining subjects at TR 1
  density_est=density(phase_df_remain$phase[phase_df_remain$TR==1])
  local_maxima <- which(diff(sign(diff(density_est$y))) == -2) + 1#detect changes in the sign of first derivative
  peak_locations <- density_est$x[local_maxima]
  
  #convert phase from [-pi, pi] to [0, 2pi]
  peak_locations=peak_locations+pi
  phase_df_remain$phase=phase_df_remain$phase+pi
  #extract the first 6 TRs
  first_6TR=phase_df[phase_df$TR<7&phase_df$SSID==sub,]
  first_6TR$phase=first_6TR$phase+pi
  
  #fit gaussians to the histogram at peak locations to extract fwhm, after removing everything left or right of the mid point between two peaks
  phase_remain_TR1=phase_df_remain$phase[phase_df_remain$TR==1]
  phase_remain_TR1_left=phase_remain_TR1[phase_remain_TR1<=mean(peak_locations)]
  phase_remain_TR1_right=phase_remain_TR1[phase_remain_TR1>mean(peak_locations)]
  gaus_fit_left=fitdistr(phase_remain_TR1_left,'normal')
  gaus_fit_right=fitdistr(phase_remain_TR1_right,'normal')
  
  #first we go over each run and look at TR1, and classify it in terms of the 3 categories (ealier than peak 1, close to peak 1, or close to peak 2)
  for (run in unique(first_6TR$run)){
    
  }
  #first we go from the first TR to the 6th TR and discard ones with a phase earlier than the earlier peak in the histogram
  #note that this has to be done in sequence since for a run with a first TR aligned to the later phase peak, the subsequent TRs will wrap around
  for (TR in c(1,2,3,4,5,6)){
    
  }
  
  first_5TR_earlypeak_runs=first_5TR$run[abs(first_5TR$phase[first_5TR$TR==1]+pi-peak_locations[1])<abs(first_5TR$phase[first_5TR$TR==1]+pi-peak_locations[2])]
  #each run needs to be considered independently, therefor we only have 1 peak to work with
  phase_leftout_runs=phase_df$phase[phase_df$SSID==sub&phase_df$TR==1]
  
  for (p in phase_leftout_runs){
    phase_differences=abs(peak_locations-p)
    
  }
  #because these are all good subjects, any shift from the above scheme is error, convert it to seconds
}
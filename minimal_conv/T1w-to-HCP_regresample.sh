#!/bin/bash

#load necessary modules
ml biology freesurfer/7.4.1

#change freesurfer environmental variable
export SUBJECTS_DIR=/oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/sourcedata/freesurfer/

#convert ANTS registration file to lta (freesurfer) T1-to-fsnative (no error, not sure if correct though...)
lta_convert --invox /oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/sub-015/anat/sub-015_from-T1w_to-fsnative_mode-image_xfm.txt --outlta /scratch/users/hyang336/sub-015_Run-01_T1-to-fsnative.lta --outreg /scratch/users/hyang336/sub-015_Run-01_T1-to-fsnative.dat --src /oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/sub-015/func/sub-015_task-GoalAttnMemTest_dir-PA_run-01_space-T1w_boldref.nii.gz --trg /oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/sourcedata/freesurfer/sub-015/mri/T1.mgz

#use mri_vol2surf (freesurfer) to register and resample GLM results in T1 to fsnative 
mri_vol2surf --src /oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/sub-015/func/sub-015_task-GoalAttnMemTest_dir-PA_run-01_space-T1w_boldref.nii.gz --out /scratch/users/hyang336/mriVol2Surf_test.mgh --srcreg /scratch/users/hyang336/sub-015_Run-01_T1-to-fsnative.lta --regheader sub-015 --hemi lh 


mri_vol2surf --src <T1 GLM output.nii> --out <fsnative GLM output.mgz> --regheader sub-015 --hemi lh
mri_vol2surf --src <T1 GLM output.nii> --out <fsnative GLM output.mgz> --regheader sub-015 --hemi rh

#convert to GIFTI
mris_convert /scratch/users/hyang336/mriVol2Surf_test.mgh /scratch/users/hyang336/mriVol2Surf_test.gii #(doesn't work yet, the result from above has 0 vertices)

/scratch/users/hyang336/retro_timing_fMRI_data/yaam/sub-015_Run_1_ResMS_fold-1/ResMS-0.1.nii
#register fsnative to fs_LR

#resample "metric" (i.e. scalar, e.g. residual map) data from fsnative to fs_LR


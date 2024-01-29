#!/bin/bash

#load necessary modules
ml biology freesurfer/7.4.1
ml biology workbench

#change freesurfer environmental variable
export SUBJECTS_DIR=/oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/sourcedata/freesurfer/

#convert ANTS registration file to lta (freesurfer) T1-to-fsnative (no need according to https://neurostars.org/t/volume-to-surface-mapping-mri-vol2surf-using-fmriprep-outputs/4079)
#lta_convert --invox /oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/sub-015/anat/sub-015_from-T1w_to-fsnative_mode-image_xfm.txt --outlta /scratch/users/hyang336/sub-015_Run-01_T1-to-fsnative.lta --outreg /scratch/users/hyang336/sub-015_Run-01_T1-to-fsnative.dat --src /scratch/users/hyang336/retro_timing_fMRI_data/yaam/sub-015_Run_1_ResMS_fold-1/ResMS-0.1.nii --trg /oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/sourcedata/freesurfer/sub-015/mri/T1.mgz

#use mri_vol2surf (freesurfer) to sample whole-brain volumatric data (in native T1) to fsnative white-matter surface 
mri_vol2surf --src /scratch/users/hyang336/retro_timing_fMRI_data/yaam/sub-015_Run_1_ResMS_fold-1/ResMS-0.1.nii --out /scratch/users/hyang336/mriVol2Surf_test_lh.func.gii --projfrac 0.5 --hemi lh --regheader sub-015
mri_vol2surf --src /scratch/users/hyang336/retro_timing_fMRI_data/yaam/sub-015_Run_1_ResMS_fold-1/ResMS-0.1.nii --out /scratch/users/hyang336/mriVol2Surf_test_rh.func.gii --projfrac 0.5 --hemi rh --regheader sub-015

#register fsnative to fs_LR

#resample "metric" (i.e. scalar, e.g. residual map) data from fsnative to fs_LR


sub='sub-031';
run=1;
fold=1; 
TR=2;
fmriprep_dir='/oak/stanford/groups/awagner/yaams-haams/fmri/yaams/bids_trimmed/derivatives/fmriprep-23.0.1/';
derivative_dir=fmriprep_dir;
behav_dir='/oak/stanford/groups/awagner/yaams-haams/fmri/yaams/onset_files/';
output='~/scratch/debug';
mask=strcat(fmriprep_dir,'/',sub,'/anat/',sub,'_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz');

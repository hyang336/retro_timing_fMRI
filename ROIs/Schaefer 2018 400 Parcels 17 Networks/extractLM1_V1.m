%extract left M1 and bilateral V1 from Schaefer2018-400Parcels17Networks in
%MNINLin2009cAsym space (from template flow)

templ=niftiread('tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-400Parcels17Networks_dseg.nii.gz');
header=niftiinfo("tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-400Parcels17Networks_dseg.nii.gz");
%according to
%https://raw.githubusercontent.com/templateflow/tpl-MNI152NLin2009cAsym/master/tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-400Parcels17Networks_dseg.tsv
%and https://github.com/ThomasYeoLab/CBIG/blob/master/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/Yeo2011_17networks_N1000.split_components.glossary.csv
% SomMotA and VisCent_ExStr are the ones roughly correspond to M1 and V1
lM1_label=[25:1:43];
V1_label=[[1:1:12],[201:1:212]];
lM1_bV1=ismember(templ,[lM1_label,V1_label]);

filename='MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-400Parcels17Networks_lSomMotA_bVisCent_ExStr.nii';
header.Filename=[pwd,'\',filename];

niftiwrite(int16(lM1_bV1),filename,header);
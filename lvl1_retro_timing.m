%% Retroactively figure out the timing between task start and scan start
% Takes in post-fmriprep functional data and fit GLMs with different
% time-offsets, then determine the correct timing based on overall R^2

% Prototypical script that only looks at one run from one subject

% For now doesn't care about subject responses, just code the goal cue and
% the stimulus as two conditions
function lvl1_retro_timing(fmriprep_dir,behav_dir,sub,output,TR,age)

%% step 1 generate alltrial regressor and noise regressor
%set up time shifting bins
% since we know that the scanner was started sometime between the first 10
% seconds (as tracked in the behaviour file), we align the behavioural file
% between [0, 10], in other words, the first pregoal fixation either
% happened at 10 seconds after the beginning of the first TR, or at 0
% second of the first TR, respectively (If you started scan after 10
% seconds, then the onset of the first pregoal fixation is 0 as all time in
% GLM is measured with respect to the start of the scan) We shift the onset
% (i.e. condition onsets with 100 ms per step). Note that the confound
% regressors need not to be shifted since those were estimated from the
% imaging data

time_to_TR1=[0:0.1:10];

%hard-coded cue duration and stim duration
cue_dr=2;
stim_dr=4;

%subject-specific output dir
if ~exist(strcat(output,'/',sub),'dir')
    mkdir (output,sub);
end
if ~exist(strcat(output,'/',sub,'_ResMS'),'dir')
    mkdir (output,[sub,'_ResMS']);
end
temp_dir=strcat(output,'/',sub,'/');
ResMS_dir=strcat(output,'/',sub,'_ResMS/');

%% Had to figure out the timing run-by-run because the delay varied across runs
runkey=fullfile(strcat(fmriprep_dir,'/',sub,'/func/'),'*GoalAttnMemTest*run-01_space-T1w*preproc_bold.nii.gz');

runfile=dir(runkey);
substr=struct();%put everythin in a struct for easy organize
substr.run=extractfield(runfile,'name');
substr.id=sub;

%unzip the nii.gz files into the temp directory
gunzip(strcat(fmriprep_dir,'/',sub,'/func/',substr.run),temp_dir);
%gunzip(maskfile,temp_dir); load the nii files, primarily to get the number
%of time points
substr.runexp=spm_vol(strcat(temp_dir,erase(substr.run,'.gz')));

%call smooth function, which is in analyses/pilot/ smooth the unzipped .nii
%files, return smoothed .nii as 1-by-run cells to a field in substr
%substr.runsmooth=crapsmoothspm(temp_dir,erase(substr.run,'.gz'),[4 4 4]);

load('retro_timing_matlabbatch_lvl1.mat');%initialized matlabbatch template MUST HAVE ALL THE NECESSARY FIELDS

%record which condtions each run has, useful for specifying design matrix
%at the end runbycond=cell(length(substr.run),2);%maximam 2 condtions per
%run

%map different naming schemes
sub_beh=sub;
sub_beh(4)='_';
sub_beh=[sub_beh(1:4),age,sub_beh(5:end)];

%load behavioural file
raw=readtable(strcat(behav_dir,'/','amass_fmri_',sub_beh,'_session_1_block_1_phase_ret_data.csv'));
%extract the two relevant timing info (goal cue and stim)
substr.runevent=[raw.goal_onset,raw.stim_onset];

%confounds
conf_name=strcat(fmriprep_dir,'/',sub,'/func/',sub,'_','*run-01_desc-confound*.tsv');%use task{1} and run{1} since it's iteratively defined
confstruct=dir(conf_name);
conffile=strcat(confstruct.folder,'/',confstruct.name);
substr.runconf=tdfread(conffile,'tab');

%build the cell structure for loading each TR of the image into matlabbatch
slice=(1:length(substr.runexp{1}));
slice=cellstr(num2str(slice'));
slice=cellfun(@strtrim,slice,'UniformOutput',false);%get rid of the white spaces
comma=repmat(',',(length(substr.runexp{1})-1+1),1);
comma=cellstr(comma);
prefix=cell(length(slice),1);
prefix(:)={substr.runexp{1}.fname};%should be the same unique run name repeated for # of TRs
%prefix=prefix';
sliceinfo=cellfun(@strcat,prefix,comma,slice,'UniformOutput',false);
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = sliceinfo;


%6 acompcor for WM and CSF, this part assumes that the confounds in the
%json file are ordered in terms of variance explained
% Modified to be compatible with fMRIprep 23.0.1, now no longer need to
% look into the json file
fn=fieldnames(substr.runconf);

%find the index for the first 6 occurance of WM and CSF in masks, also
%account for the posssibility that there might be fewer than 6
WM_ind=find(cellfun(@(x)contains(x,'w_comp_cor'),fn));
if length(WM_ind)>=6
    WM_1='w_comp_cor_00';
    WM_2='w_comp_cor_01';
    WM_3='w_comp_cor_02';
    WM_4='w_comp_cor_03';
    WM_5='w_comp_cor_04';
    WM_6='w_comp_cor_05';

    %add these components as regressors into the GLM
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'acomp_WM1';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = substr.runconf.(WM_1)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'acomp_WM2';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = substr.runconf.(WM_2)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'acomp_WM3';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = substr.runconf.(WM_3)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'acomp_WM4';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = substr.runconf.(WM_4)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).name = 'acomp_WM5';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).val = substr.runconf.(WM_5)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).name = 'acomp_WM6';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).val = substr.runconf.(WM_6)(1:end);

    w=6;%how many WM regressors we have
else
    for w=1:length(WM_ind)
        WM=fn{WM_ind(w)};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w).name = strcat('acomp_WM',num2str(w));
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w).val = substr.runconf.(WM)(1:end);
    end
end

CSF_ind=find(cellfun(@(x)contains(x,'c_comp_cor'),fn));
if length(CSF_ind)>=6
    CSF_1='c_comp_cor_00';
    CSF_2='c_comp_cor_01';
    CSF_3='c_comp_cor_02';
    CSF_4='c_comp_cor_03';
    CSF_5='c_comp_cor_04';
    CSF_6='c_comp_cor_05';

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+1).name = 'acomp_CSF1';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+1).val = substr.runconf.(CSF_1)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+2).name = 'acomp_CSF2';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+2).val = substr.runconf.(CSF_2)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+3).name = 'acomp_CSF3';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+3).val = substr.runconf.(CSF_3)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+4).name = 'acomp_CSF4';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+4).val = substr.runconf.(CSF_4)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+5).name = 'acomp_CSF5';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+5).val = substr.runconf.(CSF_5)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+6).name = 'acomp_CSF6';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+6).val = substr.runconf.(CSF_6)(1:end);
else
    for c=1:length(CSF_ind)
        CSF=fn{CSF_ind(c)};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+c).name = strcat('acomp_CSF',num2str(c));
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+c).val = substr.runconf.(CSF)(1:end);
    end
end



%specify run-agnostic fields
matlabbatch{1}.spm.stats.fmri_spec.dir = {temp_dir};%all runs are combined into one
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;%remember to change this according to actual TR in second
%matlabbatch{1}.spm.stats.fmri_spec.mask = {maskfile};%specify explicit
%mask, using avg whole-brain mask specify the to-be-output SPM.mat dir
matlabbatch{2}.spm.stats.fmri_est.spmmat = {strcat(temp_dir,'SPM.mat')};

%condition names and boxcar duration
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'goalcue';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = cue_dr;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'stim';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = stim_dr;

%gotta fill these fields too
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;

%initil setup for SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Dumb for-loop over all the possible time shift. WARNING, this fits 100+
% GLM per run (depending on the time resolution of your choice), go grab a
% coffee or rub one out
for shift=1:length(time_to_TR1)
    goal_onset=substr.runevent(:,1)-time_to_TR1(shift);
    stim_onset=substr.runevent(:,2)-time_to_TR1(shift);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = goal_onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = stim_onset;
    %run here to generate SPM.mat
    spm_jobman('run',matlabbatch(1:2));

    %rename ResMS.nii, which contains the mean square of the residual time
    %series after fitting the GLM, so it doesn't get overwritten
    resms=[temp_dir,'/ResMS.nii'];
    newname=['ResMS',num2str(time_to_TR1(shift))];
    movefile([temp_dir,'/ResMS.nii'],[ResMS_dir,'/',newname,'.nii']);
    
    delete([temp_dir,'/SPM.mat']);
end


%% plot some figures
% In a perfect world every voxel will have a single min value in terms of
% the ResMS as a function of time shifting. We may need to first filter out
% the informative voxels (e.g. those having only one minimum and the peak
% is tall), then using those voxels to select the correct timing.

end
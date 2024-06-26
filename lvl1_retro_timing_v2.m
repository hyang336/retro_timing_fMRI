%% no rand_v, otherwise the same as lvl1_retro_timing

function lvl1_retro_timing_v2(fmriprep_dir,derivative_dir,behav_dir,sub,run,output,TR,fold,mask,varargin)

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

%% input parser
% Optional input:
% 'bin_num':
% 'start_time':
% 'end_time':
p=inputParser;

%microtime resolution and bin number for matlabbatch
default.mt_res=20;
default.mt_0=10;

default.bin_num=101; %number of endpoints between the starting and the end time points (bin number + 1)

default.start_time=-5;
default.end_time=5;

%these isint functions are from the Princeton MVPA toolbox
addParameter(p,'mt_res',default.mt_res,@isint);
addParameter(p,'mt_0',default.mt_0,@isint);
addParameter(p,'bin_num',default.bin_num,@isint);
addParameter(p,'start_time',default.start_time,@isnumeric);
addParameter(p,'end_time',default.end_time,@isnumeric);

parse(p,varargin{:});
%%
time_to_TR1=linspace(p.Results.start_time,p.Results.end_time,p.Results.bin_num);

%hard-coded cue duration and stim duration
cue_dr=2;
stim_dr=4;

%subject,run-specific output dir
if ~exist(strcat(output,'/',sub,'_Run_',num2str(run)),'dir')
    mkdir (output,[sub,'_Run_',num2str(run)]);
end
if ~exist(strcat(output,'/',sub,'_Run_',num2str(run),'_ResMS_fold-',num2str(fold)),'dir')
    mkdir (output,[sub,'_Run_',num2str(run),'_ResMS_fold-',num2str(fold)]);
end

temp_dir=strcat(output,'/',sub,'_Run_',num2str(run),'/');
ResMS_dir=strcat(output,'/',sub,'_Run_',num2str(run),'_ResMS_fold-',num2str(fold));

%% Had to figure out the timing run-by-run because the delay varied across runs
%runkey=fullfile(strcat(fmriprep_dir,'/',sub,'/func/'),'*GoalAttnMemTest*run-01_space-T1w*preproc_bold.nii.gz');
runkey=fullfile(strcat(derivative_dir,'/',sub,'/func/'),strcat('*GoalAttnMemTest*run-0',num2str(run),'_space-T1w*preproc_bold.nii.gz'));

runfile=dir(runkey);
substr=struct();%put everythin in a struct for easy organize
substr.run=extractfield(runfile,'name');
substr.id=sub;

%unzip the nii.gz files into the temp directory
gunzip(strcat(fmriprep_dir,'/',sub,'/func/',substr.run),temp_dir);
%also need to unzip the mask
pattern='\.gz$';
match = regexp(mask, pattern, 'once');
if ~isempty(match)
    gunzip(mask,temp_dir);
end
%load the nii files, primarily to get the number
%of time points
substr.runexp=spm_vol(strcat(temp_dir,erase(substr.run,'.gz')));
maskfiles = dir(fullfile(temp_dir, '*brain_mask*'));
maskfile=[temp_dir,maskfiles.name];
%call smooth function, which is in analyses/pilot/ smooth the unzipped .nii
%files, return smoothed .nii as 1-by-run cells to a field in substr
%substr.runsmooth=crapsmoothspm(temp_dir,erase(substr.run,'.gz'),[4 4 4]);

load('retro_timing_matlabbatch_lvl1.mat');%initialized matlabbatch template MUST HAVE ALL THE NECESSARY FIELDS

%change microtime resolution and ref bin
matlabbatch{1, 1}.spm.stats.fmri_spec.timing.fmri_t=p.Results.mt_res;
matlabbatch{1, 1}.spm.stats.fmri_spec.timing.fmri_t0=p.Results.mt_0;

%record which condtions each run has, useful for specifying design matrix
%at the end runbycond=cell(length(substr.run),2);%maximam 2 condtions per
%run

%map different naming schemes
% sub_beh=sub;
% sub_beh(4)='_';
% sub_beh=[sub_beh(1:4),age,sub_beh(5:end)];

%load behavioural file
raw=readtable(strcat(behav_dir,'/',sub,'_onsets.csv'));
% In older version of MATLAB, columns with non-numerical data will be converted to cell, which will fuck up the later code
% Remove noresp trials and convert RT to double
if strcmp(class(raw.respRT),'cell')
    noresp_trials=find(strcmp(raw.resp,'None'));
    raw(noresp_trials,:)=[];
    raw.respRT=str2double(raw.respRT);
end
current_run=raw.run==run;
raw_run=raw(current_run,:);
%extract relevant timing info (goal cue, stim onset, resp, and respRT)
substr.runevent=[num2cell(raw_run.goal_onset),num2cell(raw_run.stim_onset),raw_run.resp,num2cell(raw_run.respRT)];
%separate key presses and calculate onsets, excluding no response
substr.right_index_RT=substr.runevent(find(strcmp(substr.runevent(:,3),'y')),4);%this is a cell
%substr.right_index_onsets=cellfun(@str2double,substr.runevent(find(strcmp(substr.runevent(:,3),'y')),4))+cell2mat(substr.runevent(find(strcmp(substr.runevent(:,3),'y')),2));%this is a mat
substr.right_index_onsets=cellfun(@(x,y) x+y,substr.runevent(find(strcmp(substr.runevent(:,3),'y')),4),substr.runevent(find(strcmp(substr.runevent(:,3),'y')),2));%this is a mat
substr.right_middle_RT=substr.runevent(find(strcmp(substr.runevent(:,3),'g')),4);%this is a cell
%substr.right_middle_onsets=cellfun(@str2double,substr.runevent(find(strcmp(substr.runevent(:,3),'g')),4))+cell2mat(substr.runevent(find(strcmp(substr.runevent(:,3),'g')),2));%this is a mat
substr.right_middle_onsets=cellfun(@(x,y) x+y,substr.runevent(find(strcmp(substr.runevent(:,3),'g')),4),substr.runevent(find(strcmp(substr.runevent(:,3),'g')),2));%this is a mat

%confounds
conf_name=strcat(fmriprep_dir,'/',sub,'/func/',sub,'_',['*run-0',num2str(run),'_desc-confound*.tsv']);%use task{1} and run{1} since it's iteratively defined
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

%6 motion regressors since we did not run ICA on these data
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'trans_x';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = substr.runconf.('trans_x')(1:end);
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'trans_y';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = substr.runconf.('trans_y')(1:end);
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'trans_z';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = substr.runconf.('trans_z')(1:end);
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'rot_x';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = substr.runconf.('rot_x')(1:end);
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).name = 'rot_y';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).val = substr.runconf.('rot_y')(1:end);
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).name = 'rot_z';
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).val = substr.runconf.('rot_z')(1:end);

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
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(7).name = 'acomp_WM1';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(7).val = substr.runconf.(WM_1)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(8).name = 'acomp_WM2';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(8).val = substr.runconf.(WM_2)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(9).name = 'acomp_WM3';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(9).val = substr.runconf.(WM_3)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(10).name = 'acomp_WM4';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(10).val = substr.runconf.(WM_4)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(11).name = 'acomp_WM5';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(11).val = substr.runconf.(WM_5)(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(12).name = 'acomp_WM6';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(12).val = substr.runconf.(WM_6)(1:end);

    w=6+6;%how many regressors we have so far
else
    for w=1:length(WM_ind)
        WM=fn{WM_ind(w)};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+6).name = strcat('acomp_WM',num2str(w));
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(w+6).val = substr.runconf.(WM)(1:end);
    end
    w=w+6;%account for the 6 motion regressor
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
matlabbatch{1}.spm.stats.fmri_spec.mask = {maskfile};%specify explicit mask, using subject-specific MNI mask by default

%condition names and boxcar duration
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'goalcue';
%matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = cue_dr;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'stim';
%matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = stim_dr;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
%2023-Dec: also add in button presses and parametrically modulate it by RT as in the Mumford paper, this wouldn't make much sense
%if we are only looking at V1
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'right_index';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'right_middle';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 0;

%gotta fill these fields too
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
%matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
%matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', 'RInd_RT', 'param', cell2mat(substr.right_index_RT), 'poly', {});
%matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
%matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
%matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', 'RMid_RT', 'param', cell2mat(substr.right_middle_RT), 'poly', {});
%matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth = 1;

% Changing the filter may not be very useful, although there may be a drift
% in the BOLD signal overtime, it is unlikely to be related to the task
% onset, moreover, the regressors in the GLM are all stationary so the
% model fit won't reflect the degree to which the regressors capture this
% drift.

%matlabbatch{1}.spm.stats.fmri_spec.sess.hpf=1280;%lower the cut-off frequency
%for high-pass filter, units in Second

%initil setup for SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Dumb for-loop over all the possible time shift. WARNING, this fits 100+
% GLM per run (depending on the time resolution of your choice), go grab a
% coffee or rub one out
for shift=1:length(time_to_TR1)
    goal_onset=cell2mat(substr.runevent(:,1))-time_to_TR1(shift);
    stim_onset=cell2mat(substr.runevent(:,2))-time_to_TR1(shift);
    RInd_onset=substr.right_index_onsets-time_to_TR1(shift);
    RMid_onset=substr.right_middle_onsets-time_to_TR1(shift);

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = goal_onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = stim_onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = RInd_onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = RMid_onset;
    %run here to generate SPM.mat
    spm_jobman('run',matlabbatch(1:2));

    %% setting up the contrasts and results/inference jobs
    % load SPM.mat
    spmmat=load(strcat(temp_dir,'SPM.mat'));
    %setup simple main effect contrasts for the four conditions (two
    %onsets, and two button presses), really just the same as the beta
    %because the contrast vectors are all just 1. This should also generate
    %the corresponding t-maps (unthresholded).

    matlabbatch{3}.spm.stats.con.spmmat = {strcat(temp_dir,'SPM.mat')};

    %use spmmat.SPM.xX.name header to find the
    %right columns
    convec=zeros(1,length(spmmat.SPM.xX.name(1,:)));

    [~,goal_main_col]=find(contains(spmmat.SPM.xX.name(1,:),'goalcue*bf(1)'));
    [~,stim_main_col]=find(contains(spmmat.SPM.xX.name(1,:),'stim*bf(1)'));
    [~,rinx_main_col]=find(contains(spmmat.SPM.xX.name(1,:),'right_index*bf(1)'));
    [~,rmid_main_col]=find(contains(spmmat.SPM.xX.name(1,:),'right_middle*bf(1)'));
    [~,constant_col]=find(contains(spmmat.SPM.xX.name(1,:),'constant'));

    %goal contrast
    convec_goal=convec;
    convec_goal(1,goal_main_col)=1;
    convec_goal(1,constant_col)=-1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'goal_vs_const';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = convec_goal;

    %stim contrast
    convec_stim=convec;
    convec_stim(1,stim_main_col)=1;
    convec_stim(1,constant_col)=-1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'stim_vs_const';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = convec_stim;

    %right index contrast
    convec_rind=convec;
    convec_rind(1,rinx_main_col)=1;
    convec_rind(1,constant_col)=-1;
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'rind_vs_const';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = convec_rind;

    %right middle contrast
    convec_rmid=convec;
    convec_rmid(1,rmid_main_col)=1;
    convec_rmid(1,constant_col)=-1;
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'rmid_vs_const';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = convec_rmid;
    
    %run the job, this will update the SPM.mat
    spm_jobman('run',matlabbatch(3));

    %rename ResMS.nii, which contains the mean square of the residual time
    %series after fitting the GLM, so it doesn't get overwritten
    resms=[temp_dir,'/ResMS.nii'];
    newname=['ResMS',num2str(time_to_TR1(shift))];
    movefile(resms,[ResMS_dir,'/',newname,'.nii']);

    % and move each of the t-maps
    t_goal=[temp_dir,'/spmT_0001.nii'];
    newname=['T_goalVSconst',num2str(time_to_TR1(shift))];
    movefile(t_goal,[ResMS_dir,'/',newname,'.nii']);

    t_stim=[temp_dir,'/spmT_0002.nii'];
    newname=['T_stimVSconst',num2str(time_to_TR1(shift))];
    movefile(t_stim,[ResMS_dir,'/',newname,'.nii']);

    t_rind=[temp_dir,'/spmT_0003.nii'];
    newname=['T_rindVSconst',num2str(time_to_TR1(shift))];
    movefile(t_rind,[ResMS_dir,'/',newname,'.nii']);

    t_rmid=[temp_dir,'/spmT_0004.nii'];
    newname=['T_rmidVSconst',num2str(time_to_TR1(shift))];
    movefile(t_rmid,[ResMS_dir,'/',newname,'.nii']);

    %delete the SPM after we got both the residual and t-maps moved
    delete([temp_dir,'/SPM.mat']);
end


end
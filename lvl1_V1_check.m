%% quick constrast to visualize if we see the typical early visual responses
function lvl1_V1_check(fmriprep_dir,derivative_dir,sub,output_dir,TR)

%hard-coded cue duration and stim duration
cue_dr=2;
stim_dr=4;

%add a htag to the sub folder based on sub, run, and rand_v so that
%multiple instances don't try to write and move the same file

%output dir
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% Had to figure out the timing run-by-run because the delay varied across runs
%runkey=fullfile(strcat(fmriprep_dir,'/',sub,'/func/'),'*GoalAttnMemTest*run-01_space-T1w*preproc_bold.nii.gz');
runkey=fullfile(strcat(derivative_dir,'/',sub,'/func/'),strcat('*GoalAttnMemTest*MNI*preproc_bold.nii.gz'));

runfile=dir(runkey);
substr=struct();%put everythin in a struct for easy organize
substr.run=extractfield(runfile,'name');
substr.id=sub;

%unzip the nii.gz files into the temp directory
%gunzip(strcat(fmriprep_dir,'/',sub,'/func/',substr.run),output_dir);

%load the nii files, primarily to get the number
%of time points
substr.runexp=spm_vol(strcat(output_dir,'/',erase(substr.run,'.gz')));
maskfile = strcat(output_dir, '/avg_MNIfunc_BrainMask.nii');
%call smooth function, which is in analyses/pilot/ smooth the unzipped .nii
%files, return smoothed .nii as 1-by-run cells to a field in substr
%substr.runsmooth=crapsmoothspm(strcat(output_dir,'/'),erase(substr.run,'.gz'),[4 4 4]);
substr.runsmooth=spm_vol(strcat(output_dir,'/smoothed_',erase(substr.run,'.gz')));
load('retro_timing_matlabbatch_lvl1.mat');%initialized matlabbatch template MUST HAVE ALL THE NECESSARY FIELDS
sess_temp=matlabbatch{1}.spm.stats.fmri_spec.sess;
%record which condtions each run has, useful for specifying design matrix
%at the end runbycond=cell(length(substr.run),2);%maximam 2 condtions per
%run

%map different naming schemes
% sub_beh=sub;
% sub_beh(4)='_';
% sub_beh=[sub_beh(1:4),age,sub_beh(5:end)];

%load behavioural file
raw=readtable(strcat(fmriprep_dir,'/',sub,'/',sub,'_onsets.csv'));

%loop over runs
for r=1:2%length(runfile)
    matlabbatch{1}.spm.stats.fmri_spec.sess(r)=sess_temp;
    current_run=raw.run==r;
    raw_run=raw(current_run,:);
    %extract the two relevant timing info (goal cue and stim)
    substr.runevent{r}=[raw_run.goal_onset,raw_run.stim_onset];

    %confounds
    conf_name=strcat(fmriprep_dir,'/',sub,'/func/',sub,'_','*run-0',num2str(r),'_desc-confound*.tsv');
    confstruct=dir(conf_name);
    conffile=strcat(confstruct.folder,'/',confstruct.name);
    substr.runconf{r}=tdfread(conffile,'tab');

    %build the cell structure for loading each TR of the image into matlabbatch
    slice=(1:length(substr.runexp{1}));
    slice=cellstr(num2str(slice'));
    slice=cellfun(@strtrim,slice,'UniformOutput',false);%get rid of the white spaces
    comma=repmat(',',(length(substr.runexp{1})-1+1),1);
    comma=cellstr(comma);
    prefix=cell(length(slice),1);
    prefix(:)={substr.runsmooth{r}.fname};%should be the same unique run name repeated for # of TRs
    %prefix=prefix';
    sliceinfo=cellfun(@strcat,prefix,comma,slice,'UniformOutput',false);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = sliceinfo;

    %6 acompcor for WM and CSF, this part assumes that the confounds in the
    %json file are ordered in terms of variance explained
    % Modified to be compatible with fMRIprep 23.0.1, now no longer need to
    % look into the json file
    fn=fieldnames(substr.runconf{r});

    %6 motion regressors since we did not run ICA on these data
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(1).name = 'trans_x';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(1).val = substr.runconf{r}.('trans_x')(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(2).name = 'trans_y';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(2).val = substr.runconf{r}.('trans_y')(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(3).name = 'trans_z';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(3).val = substr.runconf{r}.('trans_z')(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(4).name = 'rot_x';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(4).val = substr.runconf{r}.('rot_x')(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(5).name = 'rot_y';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(5).val = substr.runconf{r}.('rot_y')(1:end);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(6).name = 'rot_z';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(6).val = substr.runconf{r}.('rot_z')(1:end);

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
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(7).name = 'acomp_WM1';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(7).val = substr.runconf{r}.(WM_1)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(8).name = 'acomp_WM2';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(8).val = substr.runconf{r}.(WM_2)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(9).name = 'acomp_WM3';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(9).val = substr.runconf{r}.(WM_3)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(10).name = 'acomp_WM4';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(10).val = substr.runconf{r}.(WM_4)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(11).name = 'acomp_WM5';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(11).val = substr.runconf{r}.(WM_5)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(12).name = 'acomp_WM6';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(12).val = substr.runconf{r}.(WM_6)(1:end);

        w=6+6;%how many regressors we have so far
    else
        for w=1:length(WM_ind)
            WM=fn{WM_ind(w)};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+6).name = strcat('acomp_WM',num2str(w));
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+6).val = substr.runconf{r}.(WM)(1:end);
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

        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+1).name = 'acomp_CSF1';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+1).val = substr.runconf{r}.(CSF_1)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+2).name = 'acomp_CSF2';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+2).val = substr.runconf{r}.(CSF_2)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+3).name = 'acomp_CSF3';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+3).val = substr.runconf{r}.(CSF_3)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+4).name = 'acomp_CSF4';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+4).val = substr.runconf{r}.(CSF_4)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+5).name = 'acomp_CSF5';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+5).val = substr.runconf{r}.(CSF_5)(1:end);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+6).name = 'acomp_CSF6';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+6).val = substr.runconf{r}.(CSF_6)(1:end);
    else
        for c=1:length(CSF_ind)
            CSF=fn{CSF_ind(c)};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+c).name = strcat('acomp_CSF',num2str(c));
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(w+c).val = substr.runconf{r}.(CSF)(1:end);
        end
    end

    %condition names and boxcar duration
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).name = 'goalcue';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).duration = cue_dr;
    %matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).name = 'stim';
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).duration = stim_dr;
    %matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).onset = substr.runevent{r}(:,1);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).onset = substr.runevent{r}(:,2);

    %gotta fill these fields too
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf=128;
    matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {''};
end


%specify run-agnostic fields
matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};%all runs are combined into one
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;%remember to change this according to actual TR in second
matlabbatch{1}.spm.stats.fmri_spec.mask = {maskfile};%specify explicit mask, using subject-specific MNI mask by default

matlabbatch{2}.spm.stats.fmri_est.spmmat = {strcat(output_dir,'/SPM.mat')};


%initil setup for SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');
%run here to generate SPM.mat
spm_jobman('run',matlabbatch(1));
spm_jobman('run',matlabbatch(2));

%% simple contrast of cue vs. baseline
spmmat=load(strcat(output_dir,'/SPM.mat'));
matlabbatch{3}.spm.stats.con.spmmat = {strcat(output_dir,'/SPM.mat')};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'cue';
%use spmmat.SPM.xX.name header to find the
%right columns
convec=zeros(1,length(spmmat.SPM.xX.name(1,:)));
[~,goalcue_col]=find(contains(spmmat.SPM.xX.name(1,:),strcat('goalcue*bf(1)')));
[~,baseline_col]=find(contains(spmmat.SPM.xX.name(1,:),strcat('constant')));
convec(1,goalcue_col)=1;
%convec(1,baseline_col)=-1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = convec;

matlabbatch{3}.spm.stats.con.consess{2}=matlabbatch{3}.spm.stats.con.consess{1};
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'baseline';
%convec(1,goalcue_col)=-1;
convec(1,baseline_col)=1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = convec;

matlabbatch{3}.spm.stats.con.consess{2}=matlabbatch{3}.spm.stats.con.consess{1};
convec=zeros(1,length(spmmat.SPM.xX.name(1,:)));
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'stim';
[~,stim_col]=find(contains(spmmat.SPM.xX.name(1,:),strcat('stim*bf(1)')));
convec(1,stim_col)=1;
%convec(1,baseline_col)=-1;
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = convec;
%run the contrast 
spm_jobman('run',matlabbatch(3));
end
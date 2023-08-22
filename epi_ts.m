%% minimal preprocessing and averaging EPI time series
function epi_ts(fmriprep_dir,behav_dir,ss_list,run_list,mask,output_dir)

%movemeent threshold in mm and highpass threshold in Hz
fd_thr=1;
hpass=0.01;

%ROI mask, need to be in MNI with func resolution
brod_header=spm_vol(mask);
brod=spm_read_vols(brod_header);


%read in ss list and run list
ss_open=fopen(ss_list,'r');
SSID=textscan(ss_open,'%s', 'Delimiter', '\n');
SSID=SSID{1};%this is where all lines are stored
SSID(cellfun('isempty',SSID))=[];%get rid of empty cells/lines

run_open=fopen(run_list,'r');
runs=textscan(run_open,'%s','Delimiter','\n');
runs=runs{1};
runs(cellfun('isempty',runs))=[];


%% Generate the regressor
% Note that there will be some timing offset because of drift over a run,
% for now we are ignoring that and assuming that the each trial begins
% exactly at the start of each epoch. *the actual time recovery algorithm
% does not have this issue, it is simply a limitation of doing this raw
% time-series analyses

TR=2;
mr=20;%micro-time resolution as in SPM
mt_0=10;%micro-time sampling point to get the regressor value
num_trial_per_run=21;
trial_length=20; %trial length in Seconds
s_period=TR/mr; %sampling period (1/s_rate) of the regressor in Seconds
% boxcars with duration equal to stimulus presentation
x1=zeros(num_trial_per_run*trial_length/s_period,1);
x2=zeros(num_trial_per_run*trial_length/s_period,1);
range1=[];
range2=[];
%replace desired entries to 1
for i=1:num_trial_per_run
    trial_range_1=[(8/s_period)+((i-1)*trial_length/s_period):(10/s_period)+((i-1)*trial_length/s_period)];
    trial_range_2=[(16/s_period)+((i-1)*trial_length/s_period):(20/s_period)+((i-1)*trial_length/s_period)];
    range1=[range1,trial_range_1];
    range2=[range2,trial_range_2];
end 
x1(range1)=1;
x2(range2)=1;
% boxcars convolved with canonical HRF
p = [6 16 1 1 6 0 32];%default spm_hrf gamma parameters
hrf=spm_hrf(TR,p,mr);
x1_conv=conv(x1,hrf,"same");
x2_conv=conv(x2,hrf,"same");

%sample the convolved regressors at the sepcified delay
for samp=1:round(trial_length/TR)
    x1_samp(samp)=x1_conv(mt_0+(samp-1)*mr);
    x2_samp(samp)=x2_conv(mt_0+(samp-1)*mr);
end

%shift the regressor earlier by 5 seconds. Since we are dealing with
%periodic data, we can just cut out the first 5 seconds. 
x1_conv_5s=x1_conv(51:end);
x2_conv_5s=x2_conv(51:end);

for samp=1:round(trial_length/TR)
    x1_samp_5s(samp)=x1_conv_5s(mt_0+(samp-1)*mr);
    x2_samp_5s(samp)=x2_conv_5s(mt_0+(samp-1)*mr);
end

%save the regressor plots
figure()
plot(x1_samp,'*','linewidth',2,'Color','b');
hold on
plot(x1_samp_5s,'*','linewidth',2,'Color','r');
saveas(gcf,[output_dir,'/cue_regressor.png']);

figure()
plot(x2_samp,'*','linewidth',2,'Color','b');
hold on
plot(x2_samp_5s,'*','linewidth',2,'Color','r');
saveas(gcf,[output_dir,'/stim_regressor.png']);
%%
for i=1:length(SSID)
    epochs=[];
    for j=1:length(runs)
        %get the run file
        runkey=fullfile(strcat(fmriprep_dir,'/',SSID{i},'/func/'),strcat('*GoalAttnMemTest*run-0',num2str(runs{j}),'_space-MNI*preproc_bold.nii.gz'));
        runfile=dir(runkey);
        if ~isempty(runfile)
            filename=extractfield(runfile,'name');
            gunzip(strcat(fmriprep_dir,'/',SSID{i},'/func/',filename),output_dir);%unzip because spm
            run_header=spm_vol(strcat(output_dir,'/',erase(filename,'.gz')));
            run_tensor=spm_read_vols(run_header{1});%read in the nii file


            %Read in the confound and mark TRs with high movement
            confkey=strcat(fmriprep_dir,'/',SSID{i},'/func/',SSID{i},'_','*run-0',runs{j},'_desc-confound*.tsv');
            confstruct=dir(confkey);
            conffile=strcat(confstruct.folder,'/',confstruct.name);
            runconf=tdfread(conffile,'tab');

            for t=1:size(runconf.framewise_displacement,1)
                fd_num(t)=str2double(runconf.framewise_displacement(t,:));
            end

            fd_over=find(fd_num>=fd_thr);

            %replace TRs with high motion with NAs
            run_tensor(:,:,:,fd_over)=NaN;


            %exclude voxels outside of ROI and those that are all 0 throughout the run to reduce computational cost
            zero_voxels=~any(run_tensor,4);
            v1=logical(zeros(size(brod)));
            v1(brod==17)=1;
            %in V1 and not all zero, in case the V1 mask some regions outside
            %the brain
            v1_nozero=v1&~zero_voxels;
            roi_tensor=[];
            for u=1:size(run_tensor,4)
                vol=run_tensor(:,:,:,u);
                roi_tensor(:,:,:,u)=vol(v1_nozero);
            end
            roi_tensor=squeeze(roi_tensor);

            
            %MATLAB 2017 on Sherlock does not have the highpass() function,
            %so we have to design our own highpass filter. The results are
            %not exactly the same but close enough. 
            highpassFilter=designfilt('highpassfir','StopbandFrequency',0.001,'PassbandFrequency',0.05,'StopbandAttenuation',60,'PassbandRipple',1,'SampleRate',0.5,'DesignMethod','kaiserwin');
            roi_filtered=filtfilt(highpassFilter,roi_tensor');%using zero-phase filter, does require more samples
            %roi_filtered=highpass(roi_tensor',hpass,0.5);            


            %average across voxels, we do not calculate inter-voxel variance
            %because the errorbar needs to be across trials and later across
            %subjects
            avg_ts=mean(roi_filtered,2,'omitnan');

            %Read in event onset file
            raw=readtable(strcat(behav_dir,'/',SSID{i},'_onsets.csv'));
            current_run=raw.run==str2num(runs{j});
            raw_run=raw(current_run,:);

            %calculate the onsets in TR, need to account for drift in time
            rough_tr_pg=round(raw_run.pregoal_onset/2);
            rough_tr_cue=round(raw_run.goal_onset/2);
            %         ISI=diff(rough_tr_pg);
            %         %find where the drift becomes more than 1 sec (from round down to
            %         %round up), should only have 1
            %         drift_trial=find(ISI>10)+1;
            
            % Save mean time series and event onsets for each run
            save([output_dir,'/',SSID{i},'_run-',runs{j},'_avg_ts_goal_onsets.mat'],"avg_ts","rough_tr_cue");

            %Epoch the average time series
            for trial=1:length(rough_tr_pg)
                epochs=[epochs,avg_ts(rough_tr_pg(trial)+1:rough_tr_pg(trial)+10)];
            end
        end
    end

    %save the clean epoch data before averaging
    save([output_dir,'/',SSID{i},'_epochs.mat'],"epochs");

    %calculate mean and std across all trials
    sub_mean(:,i)=mean(epochs,2,'omitnan');
    sub_std=std(epochs,1,2,'omitnan');

    figure()
    plot([1:10],sub_mean(:,i),'+','linewidth',2);
    hold on
    plot([1:10],sub_mean(:,i)+sub_std,'+','linewidth',2);
    plot([1:10],sub_mean(:,i)-sub_std,'+','linewidth',2);
    saveas(gcf,[output_dir,'/submean_substd_',SSID{i},'.png']);

    %calculate std across runs, ignore within-run variability
    actual_run=size(epochs,2)/21;
    run_mean=[];
    for l=1:actual_run
        run_mean(:,l)=mean(epochs(:,(l-1)*21+1:l*21),2,'omitnan');
    end
    std_run=std(run_mean,1,2,'omitnan');
    figure()
    plot([1:10],sub_mean(:,i),'+','linewidth',2);
    hold on
    plot([1:10],sub_mean(:,i)+std_run,'+','linewidth',2);
    plot([1:10],sub_mean(:,i)-std_run,'+','linewidth',2);
    saveas(gcf,[output_dir,'/submean_stdrun_',SSID{i},'.png']);

    %calculate std across trials for each run
    actual_run=size(epochs,2)/21;
    run_mean=[];
    run_std=[];
    for k=1:actual_run
        run_mean=mean(epochs(:,(k-1)*21+1:k*21),2,'omitnan');
        run_std=std(epochs(:,(k-1)*21+1:k*21),1,2,'omitnan');
        figure()
        plot([1:10],run_mean,'+','linewidth',2);
        hold on
        plot([1:10],run_mean+run_std,'+','linewidth',2);
        plot([1:10],run_mean-run_std,'+','linewidth',2);
        saveas(gcf,[output_dir,'/runmean_runstd_',SSID{i},'run-0',num2str(k),'.png']);
    end

end

%grand mean and std across subj
grand_mean=mean(sub_mean,2,'omitnan');
grand_std=std(sub_mean,1,2,'omitnan');
figure()
plot([1:10],grand_mean,'+','linewidth',2);
hold on
plot([1:10],grand_mean+grand_std,'+','linewidth',2);
plot([1:10],grand_mean-grand_std,'+','linewidth',2);
saveas(gcf,[output_dir,'/grandmean_grandstd.png']);
end
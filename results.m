%% compile results and generate figures across subjects, run separately for old and young
function results(ss_list,data_dir,run_list,mask,output_dir)

%hyper param
s_time=-5;
e_time=5;
bin_num=101;
tile=linspace(s_time,e_time,bin_num);
zoom_factor=10;
num_left=floor(((bin_num-1)/zoom_factor)/2);
num_right=num_left;
median_order=8;

% Make output dir
if ~exist(strcat(output_dir,'/figures'),'dir')
    mkdir(strcat(output_dir,'/figures'));
end

%ROI mask
if ~strcmp(mask,'whole')
    brod_header=spm_vol(mask);
    brod=spm_read_vols(brod_header);
end

%read in ss list and run list
ss_open=fopen(ss_list,'r');
SSID=textscan(ss_open,'%s', 'Delimiter', '\n');
SSID=SSID{1};%this is where all lines are stored
SSID(cellfun('isempty',SSID))=[];%get rid of empty cells/lines

run_open=fopen(run_list,'r');
runs=textscan(run_open,'%s','Delimiter','\n');
runs=runs{1};
runs(cellfun('isempty',runs))=[];

% initialize data saving table
headers={'rec_error','run','sub','MetricType'};
%varTypes={'double','double','string','string'};
comp_table=cell2table(cell(0,numel(headers)),'VariableNames',headers);

% Loop over subject, run and fold to generate figures and accuracy measure.
% Need to detect whether a subject has that run through whether the output
% folder has file in it
for i=1:length(SSID)
    %the plots are always based on fold-1 to give a general idea on how the
    %GLM did over the entire search range
    for j=1:length(runs)
        %For now only looking at fold-1,detect whether a folder is empty
        res_dir=[data_dir,'/',SSID{i},'_Run_',runs{j},'_ResMS_fold-1'];
        content=dir(res_dir);
        files=content(~ismember({content.name}, {'.', '..'}));
        
        if ~isempty(files)

            volume_mean_res=[];
            volume_se_res=[];
            if ~strcmp(mask,'whole')
                ROI=find(brod==17);
            else
                ROI='whole';
            end
            for k=1:length(tile)
                tile_str{k}=sprintf('%g', tile(k));
                res_header=spm_vol([res_dir,'/ResMS',tile_str{k},'.nii']);
                res_vol=spm_read_vols(res_header);
                
                tgoal_header=spm_vol([res_dir,'/T_goalVSconst',tile_str{k},'.nii']);
                tgoal_vol=spm_read_vols(tgoal_header);
                tstim_header=spm_vol([res_dir,'/T_stimVSconst',tile_str{k},'.nii']);
                tstim_vol=spm_read_vols(tstim_header);
                trind_header=spm_vol([res_dir,'/T_rindVSconst',tile_str{k},'.nii']);
                trind_vol=spm_read_vols(trind_header);
                trmid_header=spm_vol([res_dir,'/T_rmidVSconst',tile_str{k},'.nii']);
                trmid_vol=spm_read_vols(trmid_header);

                if ~strcmp(ROI,'whole')
                    assert(all(size(res_vol)==size(brod)));%make sure the residule volumne and the ROI file has the same size otherwise the linear index will be wrong
                    ROIres=res_vol(ROI);

                    assert(all(size(tgoal_vol)==size(brod)));
                    ROItgoal=tgoal_vol(ROI);

                    assert(all(size(tstim_vol)==size(brod)));
                    ROItstim=tstim_vol(ROI);

                    assert(all(size(trind_vol)==size(brod)));
                    ROItrind=trind_vol(ROI);

                    assert(all(size(trmid_vol)==size(brod)));
                    ROItrmid=trmid_vol(ROI);
                else
                    ROIres=res_vol;
                    ROItgoal=tgoal_vol;
                    ROItstim=tstim_vol;
                    ROItrind=trind_vol;
                    ROItrmid=trmid_vol;
                end
                volume_mean_res(k)=(sum(sum(sum(ROIres,'omitnan'),'omitnan'),'omitnan'))/sum(sum(sum(~isnan(ROIres))));
                volume_se_res(k)=std(reshape(ROIres,1,[]),1,'omitnan')/sqrt(sum(sum(sum(~isnan(ROIres)))));

                volume_mean_tgoal(k)=sum(sum(sum(ROItgoal,'omitnan'),'omitnan'),'omitnan')/sum(sum(sum(~isnan(ROItgoal))));
                volume_se_tgoal(k)=std(reshape(ROItgoal,1,[]),1,'omitnan')/sqrt(sum(sum(sum(~isnan(ROItgoal)))));

                volume_mean_tstim(k)=sum(sum(sum(ROItstim,'omitnan'),'omitnan'),'omitnan')/sum(sum(sum(~isnan(ROItstim))));
                volume_se_tstim(k)=std(reshape(ROItstim,1,[]),1,'omitnan')/sqrt(sum(sum(sum(~isnan(ROItstim)))));

                volume_mean_trind(k)=sum(sum(sum(ROItrind,'omitnan'),'omitnan'),'omitnan')/sum(sum(sum(~isnan(ROItrind))));
                volume_se_trind(k)=std(reshape(ROItrind,1,[]),1,'omitnan')/sqrt(sum(sum(sum(~isnan(ROItrind)))));

                volume_mean_trmid(k)=sum(sum(sum(ROItrmid,'omitnan'),'omitnan'),'omitnan')/sum(sum(sum(~isnan(ROItrmid))));
                volume_se_trmid(k)=std(reshape(ROItrmid,1,[]),1,'omitnan')/sqrt(sum(sum(sum(~isnan(ROItrmid)))));
                
                disp(k)
            end
            %%
            %plot and save the ROI residual for the current run            
            figure;
            errorbar(tile,volume_mean_res,volume_se_res,'b', 'LineWidth', 1.5);
            title('Mean volumne residual with standard error');
            xlabel('time shift');
            ylabel('Volume-mean residual');
            xticks(linspace(s_time, e_time, 11));
            saveas(gcf,[output_dir,'/figures/',SSID{i},'_',runs{j},'_res_plot.png'])
            close;
            
            figure;
            errorbar(tile,volume_mean_tgoal,volume_se_tgoal,'b', 'LineWidth', 1.5);
            title('Mean volumne T-value for goal>baseline with standard error');
            xlabel('time shift');
            ylabel('Volume-mean T-value');
            xticks(linspace(s_time, e_time, 11));
            saveas(gcf,[output_dir,'/figures/',SSID{i},'_',runs{j},'_tgoal_plot.png'])
            close;

            figure;
            errorbar(tile,volume_mean_tstim,volume_se_tstim,'b', 'LineWidth', 1.5);
            title('Mean volumne T-value for stim>baseline with standard error');
            xlabel('time shift');
            ylabel('Volume-mean T-value');
            xticks(linspace(s_time, e_time, 11));
            saveas(gcf,[output_dir,'/figures/',SSID{i},'_',runs{j},'_tstim_plot.png'])
            close;

            figure;
            errorbar(tile,volume_mean_trind,volume_se_trind,'b', 'LineWidth', 1.5);
            title('Mean volumne T-value for right_index>baseline with standard error');
            xlabel('time shift');
            ylabel('Volume-mean T-value');
            xticks(linspace(s_time, e_time, 11));
            saveas(gcf,[output_dir,'/figures/',SSID{i},'_',runs{j},'_trind_plot.png'])
            close;

            figure;
            errorbar(tile,volume_mean_trmid,volume_se_trmid,'b', 'LineWidth', 1.5);
            title('Mean volumne T-value for right_middle>baseline with standard error');
            xlabel('time shift');
            ylabel('Volume-mean T-value');
            xticks(linspace(s_time, e_time, 11));
            saveas(gcf,[output_dir,'/figures/',SSID{i},'_',runs{j},'_trmid_plot.png'])
            close;

            %global min
            res_min_ind=find(volume_mean_res==min(volume_mean_res));
            res_min_se_ind=find(volume_se_res==min(volume_se_res));
            
            %global max for t value, but still looking for min std
            tgoal_max_ind=find(volume_mean_tgoal==max(volume_mean_tgoal));
            tgoal_min_se_ind=find(volume_se_tgoal==min(volume_se_tgoal));

            tstim_max_ind=find(volume_mean_tstim==max(volume_mean_tstim));
            tstim_min_std_ind=find(volume_se_tstim==min(volume_se_tstim));
            
            trind_max_ind=find(volume_mean_trind==max(volume_mean_trind));
            trind_min_std_ind=find(volume_se_trind==min(volume_se_trind));

            trmid_max_ind=find(volume_mean_trmid==max(volume_mean_trmid));
            trmid_min_std_ind=find(volume_se_trmid==min(volume_se_trmid));

            %check if the min or max is a locally smooth peak
            is_smooth_peak=SmoothPeakDetect(volume_mean_res,'min',res_min_ind,num_left);


            %if not a smooth peak, smooth
            if is_smooth_peak
                volume_sum_res_smooth=volume_mean_res;
                
            else
                volume_sum_res_smooth=medfilt1(volume_mean_res,median_order,'omitnan','truncate');
               
            end
            
            res_min_ind_smooth=find(volume_sum_res_smooth==min(volume_sum_res_smooth));
            
            %peak-based time estimates
            time_est_res=tile(res_min_ind);%minimal mean residual
            time_est_smooth_res=tile(res_min_ind_smooth);%minimal smoothed mean residual
            time_est_res_se=tile(res_min_se_ind);%minimal standard error
            
            time_est_tgoal=tile(tgoal_max_ind);
            time_est_tstim=tile(tstim_max_ind);
            time_est_trind=tile(trind_max_ind);
            time_est_trmid=tile(trmid_max_ind);

            %save smoothed and nonsmoothed global min
            rowdata1={time_est_res,runs{j},SSID{i},'VolResMean'};
            rowdata2={time_est_smooth_res(1),runs{j},SSID{i},'VolResMean_smooth'};%just take the 1st component since the smoothing may give min on multiple time points
            rowdata3={time_est_res_se,runs{j},SSID{i},'VolResSe'};
            rowdata4={time_est_tgoal,runs{j},SSID{i},'VolGoalTvalMean'};
            rowdata5={time_est_tstim,runs{j},SSID{i},'VolStimTvalMean'};
            rowdata6={time_est_trind,runs{j},SSID{i},'VolRIndTvalMean'};
            rowdata7={time_est_trmid,runs{j},SSID{i},'VolRMidTvalMean'};
            comp_table=[comp_table;rowdata1;rowdata2;rowdata3;rowdata4;rowdata5;rowdata6;rowdata7];
        end
    end
    disp(i)
end

%output as csv
save([output_dir,'/comp_table.mat'],'comp_table');
writetable(comp_table,[output_dir,'/comp_table.csv']);%,"WriteMode","overwrite");

end
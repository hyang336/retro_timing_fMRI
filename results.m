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


%ROI mask
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

% initialize data saving table
headers={'rec_error','run','sub','smoothed'};
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

            volume_sum_res=[];

            v1ind=find(brod==17);
            for k=1:length(tile)
                tile_str{k}=sprintf('%g', tile(k));
                res_header=spm_vol([res_dir,'/ResMS',tile_str{k},'.nii']);
                res_vol=spm_read_vols(res_header);
                assert(all(size(res_vol)==size(brod)));%make sure the residule volumne and the ROI file has the same size otherwise the linear index will be wrong
                V1res=res_vol(v1ind);
                volume_sum_res(k)=sum(sum(sum(V1res,'omitnan'),'omitnan'),'omitnan');
            end

            %global min
            min_ind=find(volume_sum_res==min(volume_sum_res));

            %right and left of the min, need to handle the cases where left index
            %may be <0 and right index may be outside the range
            left_df=diff(volume_sum_res(max(1,min_ind-num_left):min_ind));
            right_df=diff(volume_sum_res(min_ind:min(length(volume_sum_res),min_ind+num_right)));
            left_mono_dec=all(left_df<0);
            right_mono_inc=all(right_df>0);


            %if not monotonic, smooth
            if ~left_mono_dec||~right_mono_inc
                volume_sum_res_smooth=medfilt1(volume_sum_res,median_order,'omitnan','truncate');
 
            else
                volume_sum_res_smooth=volume_sum_res;
               
            end
            
            min_ind_smooth=find(volume_sum_res_smooth==min(volume_sum_res_smooth));
            
            time_est=tile(min_ind);
            time_est_smooth=tile(min_ind_smooth);

            %save smoothed and nonsmoothed global min
            rowdata1={time_est,runs{j},SSID{i},'no'};
            rowdata2={time_est_smooth(1),runs{j},SSID{i},'yes'};%just take the 1st component since the smoothing may give min on multiple time points

            comp_table=[comp_table;rowdata1;rowdata2];
        end
    end
    disp(i)
end

%output as csv
writetable(comp_table,[output_dir,'/comp_table.csv']);

end
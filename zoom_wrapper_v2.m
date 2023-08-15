%% same thing as zoom_wrapper but without rand_v

function minoffset=zoom_wrapper_v2(fmriprep_dir,derivative_dir,behav_dir,sub,run,output,TR,mask,varargin)

%% input parser
% Optional input:
% 'bin_num':
% 'start_time':
% 'end_time':
p=inputParser;

%median smoothing order
default.med_order=10;

%microtime resolution and bin number for matlabbatch
default.mt_res=20;
default.mt_0=10;

default.bin_num=101; %number of endpoints between the starting and the end time points (bin number + 1)

default.start_time=-5;
default.end_time=5;

%these isint functions are from the Princeton MVPA toolbox
addParameter(p,'med_order',default.med_order,@isint);
addParameter(p,'mt_res',default.mt_res,@isint);
addParameter(p,'mt_0',default.mt_0,@isint);
addParameter(p,'bin_num',default.bin_num,@isint);
addParameter(p,'start_time',default.start_time,@isnumeric);
addParameter(p,'end_time',default.end_time,@isnumeric);

parse(p,varargin{:});

%shorter variable names
med_order=p.Results.med_order;
mt_res=p.Results.mt_res;
mt_0=p.Results.mt_0;
bin_num=p.Results.bin_num;
s_time=p.Results.start_time;
e_time=p.Results.end_time;

%%
%The GLM is always run using the explicit fmriprep subject-specific
%whole-brain mask in MNI, with implicit mask turned off.

%depending on whether we specify another mask or pass in "whole" as the
%mask parameter, we will calculate the sum of GLM residual within the
%specified region (assuming the specified mask is Brodmann coded).

%Define fmriprep whole-brain subject mask in MNI
sub_wbrain=strcat(fmriprep_dir,'/',sub,'/anat/',sub,'_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz');

%no longer load it in at evaluation stage, but pass it as an explicit mask
%in the model estimation stage
if ~strcmp(mask,'whole')
    brod_header=spm_vol(mask);
    brod=spm_read_vols(brod_header);
end

fold=1;
zoom_in=1;
while zoom_in

    %perform time-shifting GLM
    lvl1_retro_timing_v2(fmriprep_dir,derivative_dir,behav_dir,sub,run,output,TR,fold,sub_wbrain,'bin_num',bin_num,'start_time',s_time,'end_time',e_time,'mt_res',mt_res,'mt_0',mt_0);

    res_dir=strcat(output,'/',sub,'_Run_',num2str(run),'_ResMS_fold-',num2str(fold));
    %bin_num=101;%number of time shift in total for each run for each subject
    tile=linspace(s_time,e_time,bin_num);
    zoom_factor=10;% we are gonna sample n/zoom_factor number of points centered on the min
    num_left=floor(((bin_num-1)/zoom_factor)/2);
    num_right=num_left;

    % get the volume average Res
    volume_sum_res=[];
    for i=1:length(tile)
        tile_str{i}=sprintf('%g', tile(i));
        whole_vol_header=spm_vol([res_dir,'/ResMS',tile_str{i},'.nii']);
        whole_vol=spm_read_vols(whole_vol_header);
        if ~strcmp(mask,'whole')
            volume_sum_res(i)=sum(sum(sum(whole_vol(brod==17),'omitnan'),'omitnan'),'omitnan');%brodmann 17 is V1, assuming the mask is the Brodmann mask
        else
            volume_sum_res(i)=sum(sum(sum(whole_vol,'omitnan'),'omitnan'),'omitnan');%whole-brain avg
        end
    end

    %plot(tile,volume_avg_res);

    %find the global minimum of volume residual
    min_ind=find(volume_sum_res==min(volume_sum_res));
    %min_time=tile(min_ind);

    %right and left of the min, need to handle the cases where left index
    %may be <0 and right index may be outside the range
    left_df=diff(volume_sum_res(max(1,min_ind-num_left):min_ind));
    right_df=diff(volume_sum_res(min_ind:min(length(volume_sum_res),min_ind+num_right)));
    left_mono_dec=all(left_df<0);
    right_mono_inc=all(right_df>0);

    if left_mono_dec&&right_mono_inc
        disp('monotonic on both sides, keep zooming in on global minimum')
        fold=fold+1;
        %new start and end time
        s_time=tile(max(1,min_ind-num_left));
        e_time=tile(min(length(volume_sum_res),min_ind+num_right));
        %Note that by handling the boundaries of the array the temporal
        %resolution of the >1 fold may be higher than 10x the original res


        %change microtime resolution (# of time bins per TR)
        mt_res=mt_res*zoom_factor*2;
        %change reference time bin to half of microtime resolution
        %(aligning regressor to the middle slice)
        mt_0=mt_0*zoom_factor*2;
    else%not monotonic on both sides, smooth the curve to get the minimum estimate
        zoom_in=0;
        disp('no longer monotonic on both sides, stop zooming in, start smoothing')

        %save results
        minind=find(volume_sum_res==min(volume_sum_res));
        minoffset=tile_str{minind};

        %smoothing
        
    end

end
end
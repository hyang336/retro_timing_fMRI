%% plot some figures
% In a perfect world every voxel will have a single min value in terms of
% the ResMS as a function of time shifting. We may need to first filter out
% the informative voxels (e.g. those having only one minimum and the peak
% is tall), then using those voxels to select the correct timing.
function zoom_wrapper(fmriprep_dir,behav_dir,sub,output,TR,age,varargin)

%% input parser
% Optional input:
% 'bin_num':
% 'start_time':
% 'end_time':
p=inputParser;

%microtime resolution and bin number for matlabbatch
default.mt_res=16;
default.mt_0=8;

default.bin_num=101; %number of endpoints between the starting and the end time points (bin number + 1)

default.start_time=0;
default.end_time=10;

%these isint functions are from the Princeton MVPA toolbox
addParameter(p,'mt_res',default.mt_res,@isint);
addParameter(p,'mt_0',default.mt_0,@isint);
addParameter(p,'bin_num',default.bin_num,@isint);
addParameter(p,'start_time',default.start_time,@isnumeric);
addParameter(p,'end_time',default.end_time,@isnumeric);

parse(p,varargin{:});

%shorter variable names
mt_res=p.Results.mt_res;
mt_0=p.Results.mt_0;
bin_num=p.Results.bin_num;
s_time=p.Results.start_time;
e_time=p.Results.end_time;
%%
fold=1;
zoom_in=1;
while zoom_in

    %perform time-shifting GLM
    lvl1_retro_timing(fmriprep_dir,behav_dir,sub,output,TR,age,fold,'bin_num',bin_num,'start_time',s_time,'end_time',e_time,'mt_res',mt_res,'mt_0',mt_0);

    res_dir=strcat(output,'/',sub,'_ResMS_fold-',num2str(fold));
    %bin_num=101;%number of time shift in total for each run for each subject
    tile=linspace(s_time,e_time,bin_num);
    zoom_factor=10;% we are gonna sample n/zoom_factor number of points centered on the min
    num_left=floor(((bin_num-1)/zoom_factor)/2);
    num_right=num_left;

    % get the volume average Res
    volume_avg_res=[];
    for i=1:length(tile)
        tile_str{i}=sprintf('%g', tile(i));
        volume_avg_res(i)=mean(niftiread([res_dir,'/ResMS',tile_str{i},'.nii']),"all","omitnan");
    end

    %plot(tile,volume_avg_res);

    %find the global minimum of volume residual
    min_ind=find(volume_avg_res==min(volume_avg_res));
    %min_time=tile(min_ind);

    %right and left of the min
    left_df=diff(volume_avg_res(min_ind-num_left:min_ind));
    right_df=diff(volume_avg_res(min_ind:min_ind+num_right));
    left_mono_dec=all(left_df<0);
    right_mono_inc=all(right_df>0);

    if left_mono_dec&&right_mono_inc
        disp('monotonic on both sides, keep zooming in on global minimum')
        fold=fold+1;
        %new start and end time
        s_time=tile(min_ind-5);
        e_time=tile(min_ind+5);
        %change microtime resolution (# of time bins per TR)
        mt_res=mt_res*zoom_factor*2;
        %change reference time bin to half of microtime resolution
        %(aligning regressor to the middle slice)
        mt_0=mt_0*zoom_factor*2;
    else%not monotonic on both sides, smooth the curve to get the minimum estimate
        zoom_in=0;
        disp('no longer monotonic on both sides, stop zooming in, start smoothing')
        
        %smoothing
        
    end

end
end
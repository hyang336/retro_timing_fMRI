%% plot some figures
% In a perfect world every voxel will have a single min value in terms of
% the ResMS as a function of time shifting. We may need to first filter out
% the informative voxels (e.g. those having only one minimum and the peak
% is tall), then using those voxels to select the correct timing.
res_dir='C:\Users\haozi\OneDrive\Desktop\Postdoc\Wagner\retro_timing_fMRI_data\lvl1_retro_timing\sub-019_ResMS';
num_shift=101;%number of time shift in total for each run for each subject
tile=linspace(0,10,num_shift);

volume_avg_res=[];
for i=1:length(tile)
    tile_str{i}=sprintf('%g', tile(i));
    volume_avg_res(i)=mean(niftiread([res_dir,'/ResMS',tile_str{i},'.nii']),"all","omitnan");
end

plot(tile,volume_avg_res);
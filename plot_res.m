%% The offset in the res plot could be caused by the fact that they waited for a t-trigger to start a trial, which could accumulated at max 420 ms lag through a run, which will make the regressors inaccurate in the GLM

fold=1;
rand_v=1;
run=2;
%res_dir=strcat(output,'/',sub,'_Rand_',num2str(rand_v),'_Run_',num2str(run),'_ResMS_fold-',num2str(fold));
res_dir=['C:\Users\haozi\Downloads\sub-029_Rand_1_Run_',num2str(run),'_ResMS_fold-1'];
%res_dir='C:\Users\haozi\OneDrive\Desktop\Postdoc\Wagner\retro_timing_fMRI_data\lvl1_retro_timing\sub-029_Rand_1_Run_1_ResMS_fold-1';
mask='C:\Users\haozi\OneDrive\Desktop\Postdoc\Wagner\retro_timing_fMRI_data\rbrodmann.nii';
s_time=-5;
e_time=5;
bin_num=101;
tile=linspace(s_time,e_time,bin_num);
volume_avg_res=[];
brod=niftiread(mask);
v1ind=find(brod==17);
    for i=1:length(tile)
        tile_str{i}=sprintf('%g', tile(i));
        res_vol=niftiread([res_dir,'/ResMS',tile_str{i},'.nii']);
        assert(all(size(res_vol)==size(brod)));%make sure the residule volumne and the ROI file has the same size otherwise the linear index will be wrong
        V1res=res_vol(v1ind);
        volume_avg_res(i)=mean(V1res,"all","omitnan");
    end
%median smooth
volume_avg_res_smooth=medfilt1(volume_avg_res);
figure()
plot(tile,volume_avg_res_smooth);

minind=find(volume_avg_res==min(volume_avg_res));
filename=tile_str{minind}

%%plot the raw time series of the .nii
%brod=niftiread(mask);
epi=niftiread(epi_f);
[v1x,v1y,v1z]=ind2sub(size(brod),find(brod==17));
timeseries=[];
for i=1:length(v1x)
    timeseries=[timeseries;squeeze(epi(v1x(i),v1y(i),v1z(i),:))'];
end
figure(1)%V1
plot([0:2:2*(size(timeseries,2)-1)],mean(timeseries,1));

% A quick fft on the raw time series
V1raw_FFT=fft(mean(timeseries,1));
Fs=0.5;
L=size(timeseries,2);
P2 = abs(V1raw_FFT/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(2)
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

figure(3)%whole-brain
plot([0:2:2*(size(epi,4)-1)],squeeze(mean(epi,[1,2,3]))');



%% Attempts to recover some timing info by leveraging frequency domain

%% What frequencies are contained in the regressor?
TR=2;
mr=20;%micro-time resolution as in SPM
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

%extract phase information from the boxcar regressors
x1_fft=fft(x1);
frequencies=(0:length(x1)-1)*(1/s_period/length(x1));
stem(frequencies, abs(x1_fft));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Frequency Spectrum of Composite Signal');

%extract phase information from the convolved regressors
x1_conv_fft=fft(x1_conv);
frequencies=(0:length(x1_conv)-1)*(1/s_period/length(x1_conv));
stem(frequencies, abs(x1_conv_fft));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Frequency Spectrum of Composite Signal');

%Other than the DC peak at 0, the second and later peaks in the
%frequency spectrum do appear at harmonics of 0.05 Hz

%% what's the frequency and phase content of the post-fMRIprep epi?
mask='C:\Users\haozi\OneDrive\Desktop\Postdoc\Wagner\retro_timing_fMRI_data\rbrodmann.nii';
epi_f='C:\Users\haozi\OneDrive\Desktop\Postdoc\Wagner\random_offsets_for_validation\bids_trimmed_data\sub-029\func\sub-029_task-GoalAttnMemTest_dir-PA_run-08_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz';
%plot the raw time series of the .nii
brod=niftiread(mask);
epi=niftiread(epi_f);
[v1x,v1y,v1z]=ind2sub(size(brod),find(brod==17));%V1
timeseries=[];
for i=1:length(v1x)
    timeseries=[timeseries;squeeze(epi(v1x(i),v1y(i),v1z(i),:))'];
end
figure(1)
plot([0:2:2*(size(timeseries,2)-1)],mean(timeseries,1));

%FFT of the EPI BOLD, the peak position change alot between run 1 and 8 for
%sub-029
epi_fft=fft(mean(timeseries,1)-mean(mean(timeseries,1)));%remove DC before FFT
frequencies=(0:length(epi_fft)-1)*(0.5/length(epi_fft));
stem(frequencies, abs(epi_fft));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Frequency Spectrum of BOLD');




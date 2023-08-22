%% simulate the filter effect

%generate a sythesized signal, with a 0.1 Hz sine component, a Gaussian
%noise component, and a linear drift
t=0:2:450;
puresig=sin(2*pi*0.1*t);
x=sin(2*pi*0.1*t)+normrnd(0,1,1,length(t))-0.01*t;

figure;
plot(x);
title('simulated data');
xlabel('time (s)');
ylabel('Amplitude (a.u.)');

%plot the frequency and phase
fft_x=fft(x);
N = length(x);  % Number of samples
f = (-N/2:N/2-1) * 0.5 / N;

figure;
plot(f, abs(fftshift(fft_x)));%0 frequency at the center
title('Frequency Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

highpassFilter=designfilt('highpassfir','StopbandFrequency',0.001,'PassbandFrequency',0.05,'StopbandAttenuation',60,'PassbandRipple',1,'SampleRate',0.5,'DesignMethod','kaiserwin');
x_filtered=filtfilt(highpassFilter,x');

fft_filt_x=fft(x_filtered);
N = length(x_filtered);  % Number of samples
f = (-N/2:N/2-1) * 0.5 / N;

figure;
plot(f, abs(fftshift(fft_filt_x)));%0 frequency at the center
title('Frequency Magnitude');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


figure;
plot(puresig)
hold on
plot(x_filtered)
title('pure signal and filtered data');
xlabel('time (s)');
ylabel('Amplitude (a.u.)');
legend('pure signal','highpassed data');
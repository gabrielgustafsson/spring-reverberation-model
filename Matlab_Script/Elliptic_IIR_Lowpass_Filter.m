%% Tenth Order Elliptic IIR Low Pass Filter Design

%% Parameters
n = 10; % order of the filter
passband_ripple = 1; % in -dB
stopband_attenuation = 60; %stopband suppresion in -dB
f = 4750; % cut-off frequency in Hz
fs = 44100; % sampling frequency in Hz

format longEng % display the full numbers

% get the filter coefficients
[B,A] = ellip(n,passband_ripple,stopband_attenuation,f/(0.5*fs),'low')


%% Plot the frequency response
figure
freqz(B,A,fs,fs)
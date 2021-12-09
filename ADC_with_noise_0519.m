function [ Time_d,S_D_noise,S_D,length_SD,SNR ] = ADC_with_noise_0519( A_times,A_acc,fmax_c,T )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

	% Initial Setting
	Amp = 2*9.8; % Maximum amplitude of the original "analog" signal,unit:m/s/s
	fmax_c = fmax_c; % Maximum frequency (Hz) of the original "analog" signal (the maximum frequency is technically infinite!)
	fs_c = fmax_c; % A large number to mimic analog signal in discrete-time
	dt_c = 1/fs_c;
	T = T; % Duration of the signal
	N_c = fix(T/dt_c);
	CL = 2.0*9.8; % Clipping level, unit:m/s/s
	Resolution = 16; % The resolution of the ADC (typically 1 or 2 bits greater than technically needed resolution to consider noise)
	% !!!!!note that value of the noise RMS in a specific frequency band, density by the sqrt(frequency band)
    freq_noise = 14; % The frequency corresponding to
    NoiseLevel = 3.9*10^(-7)*9.8*sqrt(fmax_c)/sqrt(freq_noise); % noise level of sensor, unit:m/s/s
	fs_d = fmax_c; % The sampling rate (Hz) of discrete-time signal
	DS = round(fs_c/fs_d); % to make the data size after resample is int
	fmax_d = 0.4*fs_d/2; % The maximum frequency (corner frequency) of the discrete-time signal. It is technically fs_d/2, but signals are usually sampled at a higher rate. 

	% Analog Signal
	Time_c = A_times;
	S_c = A_acc;

	lpFilt = designfilt('lowpassiir','FilterOrder',8,'PassbandFrequency',fmax_d,'PassbandRipple',0.2,'SampleRate',fs_c);
	S_c_Filt = filtfilt(lpFilt,S_c); % "Analog" signal after built-in low-pass filtering (before sampling and quantization)
 
	% Discrete-Time Signal
	dt_d = dt_c*DS;
	N_d = N_c/DS;
	Time_d = [0:dt_d:(N_d-1)*dt_d]';
	if DS<=1
		S_d = resample(S_c_Filt,DS,1);
	else
		S_d = downsample(S_c_Filt,DS);
	end
	consist_len = min(length(Time_d),length(S_d));
	Time_d = Time_d(1:consist_len);
	S_d = S_d(1:consist_len);

	% Check: Reconstructed signal (we must be able to reconstruct the analog signal)
	S_d_re = resample(S_d,DS,1);
	consist_len = min(length(Time_c),length(S_d_re));
	Time_c = Time_c(1:consist_len);
	S_d_re = S_d_re(1:consist_len);
	S_c_Filt = S_c_Filt(1:consist_len);

	% Digital Signal
	AL = linspace(-CL,CL,2^Resolution)';
	for k=1:length(S_d)
		[~,id] = min(abs(AL-S_d(k)));
		S_D(k) = AL(id);
	end
	consist_len = min(length(Time_d),length(S_D));
	Time_d = Time_d(1:consist_len);
	S_D = S_D(1:consist_len)';

	% Digital Signal with noise
	S_D_wgn = NoiseLevel*randn(length(S_D),1); % generate Gaussian white noise, with NoiseLevel 
	Ps=sum(abs(S_D).^2) ;%signal power
	Pn=sum(sum((S_D_wgn).^2));%noise power
	SNR=10*log10(Ps/Pn); % signal to noise ratio
	S_D_noise = S_D + S_D_wgn;
	length_SD = length(S_D_noise);

end


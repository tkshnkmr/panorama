clc
clear
close all

% ====================
% Input data
% ====================

N = 50000; % length of data
L = 500; % length of block data
noverlap = 0.5; % overlap for panorama
K = round(L*noverlap);      % length of overlap
M = floor((N-L)/(L-K))+1;   % Number of segments after
fs = 1; % sampling freq
t = (0:1/fs:(N-1)/fs)';
nfft = 1000;


% ====================
% filter parameters
% ====================
b = 1;
a = [1, -1.8*cos(2*pi*0.2), 0.64];

% ====================
% SNR (dB)
% ====================
parameter_snr = -5;

% % ====================
% % Signal - sum of multiple sinusoids
% % ====================
y1 = cos(2*pi*0.15*t-pi/6);
y2 = 0.25*cos(2*pi*0.25*t+pi/3);
y3 = 0.1*cos(2*pi*0.4*t+pi/8);
str = 'cos(0.15{\omega}t-\pi/6)+1/4*cos(0.25{\omega}t+\pi/3)+1/10*cos(0.4{\omega}t+\pi/8)+noise';
y_sig = y1 + y2 + y3;
% % ====================
% % Signal - sum of multiple sawtooth
% % ====================
% y1 = sawtooth(2*pi*0.15*t-pi/4);
% y2 = 0.25*sawtooth(2*pi*0.25*t+pi/3);
% y3 = 0.1*sawtooth(2*pi*0.4*t+pi/8);
% str = 'sawtooth(0.15{\omega}t-\pi/6)+1/4*sawtooth(0.25{\omega}t+\pi/3)+1/10*sawtooth(0.4{\omega}t+\pi/8)+noise';
% y_sig = y1 + y2 + y3;

% ====================
% Noise - depends on SNR (dB)
% ====================
noise = filter(b,a,randn(N,1));
sigPower = sum(abs(y_sig).^2)/N; % caluculate signal power
noisePower = sum(abs(noise).^2)/N; % calculate noise power
% add noise with certain SNR
noise = noise*sqrt(sigPower/noisePower)*sqrt(10^(-parameter_snr/10));

% ====================
% Input signal
% ====================
y = y_sig + noise;


% ====================
% Panorama
% ====================

[pxx_wel,f_wel] = pwelch(y,L,noverlap,nfft,fs);
pxx_wel = pxx_wel*fs/2/pi;

figure
[pano_func,freq_func_fft,Phase_func_fft] = panorama(y,L,noverlap,nfft,fs,'fft');
subplot(3,1,1)
plot(f_wel*2,10*log10(pxx_wel),':k'); hold on
plot(freq_func_fft*2,10*log10(abs(pano_func)),'Color',[1,0.5,0]); hold on
str = sprintf('N=%d,M=%d,Over=%.2f,nfft=%d,With Hamm win,SNR=%.2fdB',N,M,noverlap,nfft,snr(y_sig,noise));
title(sprintf('fft-based,%s',str));
legend('Pwelch','Panorama')
xlabel('normalised freq');ylabel('dB');

[pano_func,freq_func_conv] = panorama(y,L,noverlap,nfft,fs,'conv');
subplot(3,1,2)
plot(f_wel*2,10*log10(pxx_wel),':k'); hold on
plot(freq_func_conv*2,10*log10(abs(pano_func)),'Color',[0.5,0.5,0.5]); hold on
str = sprintf('N=%d,M=%d,Over=%.2f,nfft=%d,With Hamm win,SNR=%.2fdB',N,M,noverlap,nfft,snr(y_sig,noise));
title(sprintf('conv-based,%s',str));
legend('Pwelch','Panorama')
xlabel('normalised freq');ylabel('dB');

[pano_func,freq_func_sing,Phase_func_sing,Theta_func_sing] = panorama(y,L,noverlap,nfft,fs,'sing');
subplot(3,1,3)
plot(f_wel*2,10*log10(pxx_wel),':k'); hold on
plot(freq_func_sing*2,10*log10(abs(pano_func)),'Color',[0,0.5,1]); hold on
str = sprintf('N=%d,M=%d,Over=%.2f,nfft=%d,With Hamm win,SNR=%.2fdB',N,M,noverlap,nfft,snr(y_sig,noise));
title(sprintf('single realisation,%s',str));
xlabel('normalised freq');ylabel('dB');
legend('Pwelch','Panorama')
set(gcf,'Position',[0 0 1300 700]);



figure
thre = 15;
% ====================
% Phase plot for fft-based
% ====================
subplot(2,1,1)
m_Phase_X = mean(cos(Phase_func_fft),2);
std_Phase_x = std(cos(Phase_func_fft),[],2);
plot(freq_func_fft*2,m_Phase_X(1:nfft/2),'--','Color',[1,0.5,0],'Marker','*','linewidth',1); hold on;
locs = find(std_Phase_x < 2*pi/thre);
for i = 1:length(locs)
    plot(freq_func_fft(locs)*2,m_Phase_X(locs),'*','Color',[1,0.5,0],'linewidth',8);hold on;
end
title(sprintf('cos(Phase)-Thre=2/%d',thre))
ylim([-1 1]);
xlabel('normalised freq');

% ====================
% Theta plot for single-based 
% ====================
subplot(2,1,2)
m_Theta_X = mean(cos(Theta_func_sing),2);
std_Theta_x = std(cos(Theta_func_sing),[],2);
plot(freq_func_sing*2,m_Theta_X(1:nfft/2),'--','Color',[0,0.5,1],'Marker','*','linewidth',1); hold on;
locs = find(std_Theta_x < 2*pi/thre);
for i = 1:length(locs)
    plot(freq_func_sing(locs)*2,m_Theta_X(locs),'*','Color',[0,0.5,1],'linewidth',8);hold on;
end
title(sprintf('cos(Theta)-Thre=2/%d',thre))
ylim([-1 1]);
xlabel('normalised freq');
set(gcf,'Position',[0 0 1300 700]);

%% Apply to EEG data
clc
clear
close all

load('eeg_sample.mat');

y = filtered_eeg;
N = length(filtered_eeg);
L = N/10;
noverlap = 0.5;
nfft = L;
K = round(L*noverlap);      % length of overlap
M = floor((N-L)/(L-K))+1;   % Number of segments after

% ====================
% Panorama
% ====================
[pxx_wel,f_wel] = pwelch(y,L,noverlap,nfft,fs);
pxx_wel = pxx_wel*fs/2/pi;
[pano_fft,freq_func_fft,Phase_func_fft] = panorama(y,L,noverlap,nfft,fs,'fft');
[pano_conv,freq_func_conv] = panorama(y,L,noverlap,nfft,fs,'conv');
[pano_sing,freq_func_sing,Phase_func_sing,Theta_func_sing] = panorama(y,L,noverlap,nfft,fs,'sing');

% ====================
% Plot pwelch and panorama
% ====================
figure
subplot(4,1,1:2)
plot(f_wel,10*log10(pxx_wel),':k'); hold on
plot(freq_func_fft,10*log10(abs(pano_fft)),'Color',[1,0.5,0]); hold on
plot(freq_func_conv,10*log10(abs(pano_conv)),'Color',[0.5,0.5,0.5]); hold on
plot(freq_func_sing,10*log10(abs(pano_sing)),'Color',[0,0.5,1]); hold on
title(sprintf('EEG from Ear, N=%d,M=%d,Over=%.2f,nfft=%d,With Hamm win',N,M,noverlap,nfft));
legend('Pwelch','FFT-based','Conv-based','Single Realisation','Orientation','horizontal','location','southwest');
xlim([25 45])

thre = 15;
% ====================
% Phase plot for fft-based
% ====================
subplot(4,1,3)
m_Phase_X = mean(cos(Phase_func_fft),2);
std_Phase_x = std(cos(Phase_func_fft),[],2);
plot(freq_func_fft*2,m_Phase_X(1:nfft/2),'--','Color',[1,0.5,0],'Marker','*','linewidth',1); hold on;
locs = find(std_Phase_x < 2*pi/thre);
for i = 1:length(locs)
    plot(freq_func_fft(locs)*2,m_Phase_X(locs),'*','Color',[1,0.5,0],'linewidth',8);hold on;
end
title(sprintf('cos(Phase)-Thre=2/%d',thre))
xlim([25 45]);ylim([-1 1]);

% ====================
% Theta plot for single-based 
% ====================
subplot(4,1,4)
m_Theta_X = mean(cos(Theta_func_sing),2);
std_Theta_x = std(cos(Theta_func_sing),[],2);
plot(freq_func_sing*2,m_Theta_X(1:nfft/2),'--','Color',[0,0.5,1],'Marker','*','linewidth',1); hold on;
locs = find(std_Theta_x < 2*pi/thre);
for i = 1:length(locs)
    plot(freq_func_sing(locs)*2,m_Theta_X(locs),'*','Color',[0,0.5,1],'linewidth',8);hold on;
end
title(sprintf('cos(Theta)-Thre=2/%d',thre))
xlim([25 45]);ylim([-1 1]);
xlabel('Freq');
set(gcf,'Position',[0 0 1300 700]);



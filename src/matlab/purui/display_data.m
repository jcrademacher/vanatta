size = [7 7000000];

filename = "fixed_vanatta4x2_nostag_006B_006F_006A_006C_x_001A_004A_004B_004D_txfmr_nicktb_18,5kfc_0,0deg_8bit_pre_16bit_dat_prbs_0,5kbps_usrp_2,5m_depth_005B_purui_tx_60Vrms_81m_81m_1m_2foam_sep_purui_rx_0.dat";

id = fopen(strcat('/home/jradema/Documents/sk/oceans/rx_outputs/River_PAB2_Van_Atta_01-11-2023/',filename),'r');
sig = fread(id,size,'float32').';

ch1 = sig(6e4:end,4);

fs = 192e3;
fc = 18.5e3;

% plot channel 1
%sig = rx_baseband;
t = [0:1/fs:length(ch1)/fs-1/fs];
%sig = awgn(cos(2*pi*fc*t),30);
figure(1);
hold on;
plot(ch1);


window_size = 7e2;
window = chebwin(window_size, 300);
Nfft = 2^10;

[pxx,f] = pwelch(ch1,window,[],Nfft,fs);
% 
% max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
% [maxval,mindex] = max(pxx(max_search)); % max in each row
% carrier_freq = fs/Nfft*max_search(mindex)';


figure(2);

plot(f/1e3,10*log10(pxx)+180);
xlabel("Freq (kHz)");
ylabel("dB re 1uPa^2/Hz");
% xlim([1e3 100e3]);
hold on;
%%
% lpFilt = designfilt('lowpassfir' ...
%                     ,'PassbandFrequency',20e3*2/(fs)...
%                     ,'StopbandFrequency',30e3*2/(fs),'StopbandAttenuation' ...
%                     ,200,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

dfac = 4;
%new_sig = fftfilt(lpFilt,sig);
new_sig=sig;
new_sig = decimate(new_sig,dfac);
fs = fs/dfac;

window_size = floor(length(new_sig)/1);
window = chebwin(window_size);
Nfft = 2^nextpow2(fs*10);

[pxx,f] = pwelch(new_sig,window,[],Nfft,fs);

% max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
% [maxval,mindex] = max(pxx(max_search)); % max in each row
% carrier_freq = fs/Nfft*max_search(mindex)';

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");



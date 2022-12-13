sig_0 = read_complex_binary('../../rx_outputs/River PAB2 Van Atta 12-08-2022/fixed_single_006F_chest_txfmr_nicktb_siggen_18,5kfc_65,0deg_8bit_pre_16bit_dat_prbs_0,5kbps_usrp_2,5m_depth_005B_purui_tx_60Vrms_6m_5m_hphydro_diff_0.dat');
sig = real(sig_0(24:end))-imag(sig_0(24:end));
fs = 2e5;
fc = 18.5e3;

t = [0:1/fs:length(sig)/fs-1/fs];
sig = awgn(cos(2*pi*fc*t),30);
figure(1);
hold on;
plot(sig);


window_size = floor(length(sig)/1);
window = chebwin(window_size);
Nfft = 2^nextpow2(fs*10);

[pxx,f] = pwelch(sig,window,[],Nfft,fs);
% 
% max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
% [maxval,mindex] = max(pxx(max_search)); % max in each row
% carrier_freq = fs/Nfft*max_search(mindex)';

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");
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



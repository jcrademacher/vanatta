sig_0 = read_complex_binary('../../rx_outputs/River PAB2 Van Atta 12-08-2022/fixed_vanatta_006F_006B_chest_txfmr_nicktb_siggen_18,5kfc_-70,0deg_8bit_pre_16bit_dat_prbs_0,5kbps_usrp_2,5m_depth_005B_purui_tx_60Vrms_6m_5m_hphydro_diff_0.dat');
sig = real(sig_0(24:end))-imag(sig_0(24:end));
fs = 2e5;
t = [0:1/fs:length(sig)/fs-1/fs];
figure(1);
hold on;
plot(sig);
fc = 18.5e3;

window_size = floor(length(sig)/1);
window = chebwin(window_size);
Nfft = fs*10;

[pxx,f] = pwelch(sig,window,[],Nfft,fs);
% 
% max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
% [maxval,mindex] = max(pxx(max_search)); % max in each row
% carrier_freq = fs/Nfft*max_search(mindex)';

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");

lpFilt = designfilt('lowpassfir' ...
                    ,'PassbandFrequency',20e3*2/(fs)...
                    ,'StopbandFrequency',25e3*2/(fs),'StopbandAttenuation' ...
                    ,80,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

new_sig = fftfilt(lpFilt,sig);
% new_sig = downsample(new_sig,5);
% fs = fs/5;

window_size = floor(length(new_sig)/1);
window = chebwin(window_size);
Nfft = fs*10;

[pxx,f] = pwelch(new_sig,window,[],Nfft,fs);

% max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
% [maxval,mindex] = max(pxx(max_search)); % max in each row
% carrier_freq = fs/Nfft*max_search(mindex)';

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");



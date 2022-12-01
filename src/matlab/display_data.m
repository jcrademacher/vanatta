sig_0 = read_complex_binary('../../rx_outputs/WHOI Experiments 12-01-2022/fixed_single_chest_006A_txfmr_nicktb_18,5kfc_0,0deg_8bit_pre_16bit_dat_prbs_0,5kbps_usrp_3m_depth_005B_purui_tx_60Vrms_1,9m_1m_hphydro_diff_0.dat');
sig = real(sig_0(24:end))-imag(sig_0(24:end));
fs = 2e5;
t = [0:1/fs:length(sig)/fs-1/fs];
figure(1);
hold on;
plot(sig);
fc = 18.5e3;

window_size = floor(length(sig)/10);
window = chebwin(window_size);
Nfft = fs*10;

[pxx,f] = pwelch(sig,window,[],Nfft,fs,'power');

max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
[maxval,mindex] = max(pxx(max_search)); % max in each row
carrier_freq = fs/Nfft*max_search(mindex)';

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");


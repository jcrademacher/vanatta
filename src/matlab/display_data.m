sig_0 = read_complex_binary('../../rx_outputs/River PAB Van Atta 4 09-15-2022/rx_single_chest_pab_011B_stag9cm_7cm_sp_2,9mtxfmr_+180deg_mosfet_18,5kfc_prbs_0,5kbps_usrp_2,5m_depth_010A_purui_new_tx_6m_5m_hphydro_300mVpp_0.dat');
sig = real(sig_0(20:end))-imag(sig_0(20:end));
fs = 2e5;
t = [0:1/fs:length(sig)/fs-1/fs];
figure(1);
plot(t,sig);

window_size = floor(length(sig)/5);
window = chebwin(window_size);

[pxx,f] = pwelch(sig,window,[],[],fs,'power');

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");


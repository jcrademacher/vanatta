sig = read_complex_binary('../../rx_outputs/River PAB Van Atta 4 08-23-2022/rx_single_chest_pab_011B_nostag_7cm_sp_2,9mtxfmr_+0deg_nx5_18,5kfc_prbs_1kbps_usrp_2,5m_depth_010A_purui_tx_6m_5m_hphydro_0.dat');
sig = real(sig(20:end));
fs = 2e5;
figure(1);
plot(sig);

window_size = floor(length(sig)/5);
window = chebwin(window_size);

[pxx,f] = pwelch(sig,window,[],[],fs);

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");


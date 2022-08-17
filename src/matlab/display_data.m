sig = read_complex_binary('../../rx_outputs/River PAB Van Atta 4 08-15-2022/rx_single_chest_pab_008A_7cm_sp_ind1,8m_+0deg_mosfet_18,5kfc_siggen_data_1kbps_usrp_2,5m_depth_6m_u2b_5m_hphydro_1.dat');
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


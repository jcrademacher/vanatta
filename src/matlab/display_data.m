sig = read_complex_binary('../../rx_outputs/River PAB Van Atta 4 08-18-2022/rx_vanatta4_chest_pab_008A_011B_011A_010B_7cm_sp_ind1,5m_iso_-90deg_nx5_18,5kfc_nonprbs_1kbps_usrp_2,5m_depth_3m_u2b_2m_hphydro_0.dat');
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


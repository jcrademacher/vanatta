sig = read_complex_binary('../../rx_outputs/River PAB Channel Estimate 08-04-2022/rx_array_chest_pab_008A_010B_7cm_sp_ind_+75deg_ts3a_18,5kfc_siggen_data_1kbps_usrp_2,5m_depth_2m_u2b_1m_hphydro_0.dat');
sig = sig(1:end);
plot(real(sig));
figure;
periodogram(real(sig));
max(real(sig))-min(real(sig))


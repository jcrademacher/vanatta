sig = read_complex_binary('../../rx_outputs/River PAB Channel Estimate 07-13-2022/rx_vanatta_chest_pab_005A_007B_ind_0deg_tmux_18,5kfc_siggen_data_usrp_3m_depth_3m_u2b_2m_hphydro_0.dat');
sig = sig(1:end);
plot(real(sig));
figure;
periodogram(real(sig));
max(real(sig))-min(real(sig))


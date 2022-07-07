sig = read_complex_binary('../../rx_outputs/River PAB Van Atta 07-07-2022/rx_vanatta_pab_007B_005A_ind_0deg_tmux_18,5kfc_1kmod_siggen_3m_depth_3m_u2b_2m_hphydro_0.dat');
sig = sig(1:end);
plot(real(sig));
figure;
periodogram(real(sig));
max(real(sig))-min(real(sig))


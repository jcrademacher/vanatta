sig = read_complex_binary('../../rx_outputs/River PAB Van Atta 06-23-2022/rx_backscatter_vanatta_007B_003A_purui_match_0deg_tmux_21,2kfc_1kmod_2m_depth_2m_u2b_1m_midpower_hphydro_diff_0.dat');

plot(real(sig));
figure;
periodogram(real(sig));


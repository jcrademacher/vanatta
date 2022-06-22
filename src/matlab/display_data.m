sig = read_complex_binary('../../rx_outputs/River Switch Van Atta Tests 06-20-2022/rx_backscatter_vanatta_switch_TS5A_pab_003A_006A_ind_0deg_18,5kfc_1kmod_2m_depth_0.dat');

plot(real(sig));
periodogram(real(sig));


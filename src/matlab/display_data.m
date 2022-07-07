sig = read_complex_binary('../../rx_outputs/SMAST Tests 06-24-2022/rx_u2b_tx_prop_loss_18,5kfc_3m_depth_100cm_distance_hydrophone_0.dat');
sig = sig(1:end);
plot(real(sig));
figure;
periodogram(real(sig));
max(real(sig))-min(real(sig))


sig = read_complex_binary('../../rx_outputs/River PAB Channel Estimate 07-20-2022/60Hz_noise_test_004A_conn_008A_conn_tx_conn.dat');
sig = sig(1:end);
plot(real(sig));
figure;
periodogram(real(sig));
max(real(sig))-min(real(sig))


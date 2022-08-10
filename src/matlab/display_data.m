sig = read_complex_binary('../../rx_outputs/River PAB Channel Estimate 07-28-2022/60Hz_test_rx_tx_u2b_0.dat');
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


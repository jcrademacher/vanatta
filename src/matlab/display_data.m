sig_0 = read_complex_binary('../../rx_outputs/River PAB Van Atta 4 09-16-2022/waleed_test_011B_500bps_160mVpp.dat');
sig = real(sig_0(20:end));
fs = 2e5;
t = [0:1/fs:length(sig)/fs-1/fs];
figure(1);
plot(t,sig);

window_size = floor(length(sig)/5);
window = chebwin(window_size);

[pxx,f] = pwelch(sig,window,[],[],fs,'power');

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");


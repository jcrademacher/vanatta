sig_0 = read_complex_binary('../../rx_outputs/test/test_rx.dat');
sig = real(sig_0(24:end));
fs = 2e5;
t = [0:1/fs:length(sig)/fs-1/fs];
figure(1);
plot(sig);
fc = 20e3;

window_size = floor(length(sig));
window = chebwin(window_size);
Nfft = fs*1000;

[pxx,f] = pwelch(sig,window,[],Nfft,fs,'power');

max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
[maxval,mindex] = max(pxx(max_search)); % max in each row
carrier_freq = fs/Nfft*max_search(mindex)';

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");


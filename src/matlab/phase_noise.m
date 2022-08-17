
sig = real(read_complex_binary('../../rx_outputs/River PAB Channel Estimate 07-28-2022/60Hz_test_rx_tx_u2b_0.dat'));
fs = 2e5;

sig = sig(16e4:end);

t = [0:1/fs:length(sig)/fs-1/fs];

window_size = floor(length(sig)/5);
window = chebwin(window_size,120);
Nfft = 2^nextpow2(window_size);

rbw = enbw(window,fs);

[pxx,f] = pwelch(sig,window,[],Nfft,fs);
pxx = pxx/50;

% figure(1);
% plot(sig);

figure(2);
hold on;
plot(f/1e3,10*log10(pxx));
xlabel("Freq (kHz)");

% compute phase noise
offsets = [ceil(3*rbw/10)*10:20:5e3];

[carrier_val,carrier_index] = max(pxx);
carrier_freq = f(carrier_index);

phasenoise = 10*log10(pxx(round(carrier_index+offsets/fs*Nfft))/carrier_val);
figure(3);
hold on;
plot(offsets,phasenoise);
grid on;
grid minor;
xlabel("Freq Offset from Carrier (Hz)");
ylabel("Phase noise (dBc/Hz)");
title(strcat("Phase Noise w/ Carrier Freq = ",num2str(carrier_freq/1e3)," kHz"));





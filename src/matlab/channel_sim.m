fs = 2e5;
fc = 20e3;
fb = 1e3;

fpb1_lp = 7e3;
fsb1_lp = 9e3;

% lowpass after downconversion, before downsampling
lpFilt = designfilt('lowpassfir' ...
                    ,'PassbandFrequency',fpb1_lp*2/fs...
                    ,'StopbandFrequency',fsb1_lp*2/fs,'StopbandAttenuation' ...
                    ,80,'PassbandRipple',0.1);

m = 0.1;
t_len = 1;
t = [0:1/fs:t_len-1/fs];

data = square(2*pi*fb*t,40);
carrier = cos(2*pi*fc*t);

pb_sig = (1+m*data).*carrier;
pb_sig = pb_sig';

load 3m_pb_ch.mat;
channel = ch_pb2;

out = conv(pb_sig,channel);

t_out = [0:1/fs:length(out)/fs-1/fs];
downconv_out = fftfilt(lpFilt, cos(2*pi*fc*t_out)'.*out);

subplot(2,1,1);
plot(pb_sig);
ylim([-2 2]);
subplot(2,1,2);
plot(downconv_out);

t_window = 0.1;
window = chebwin(t_window*fs);
Nfft = 2^nextpow2(length(window));

figure;
hold on;
[pxx,f] = pwelch(data,window,[],Nfft,fs,'power');
plot(f,10*log10(pxx));
grid on;
grid minor;
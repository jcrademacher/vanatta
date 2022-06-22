fs = 2e5;
fc=18.5e3;

fb=1e3;

amp = 0.3;

t_length = 1;
t = [0:1/fs:t_length-1/fs];

r = amp*sin(2*pi*fc*t)+1j*(square(2*pi*fb*t)+1)/2.5;

write_complex_binary(r, strrep("../../tx_outputs/combo_tone_18,5kfc_600mVpp_square_1k_800mVpp_?sec.dat","?",num2str(t_length)));
hold on;
plot(real(r));
plot(imag(r));
fs = 2e5;
fc=50e3;

amp = 0.3;

t_length = 60;
t = [0:1/fs:t_length-1/fs];

y = amp*sin(2*pi*fc*t);

write_complex_binary(y, strrep("../../tx_outputs/tone_80kfc_600mVpp_?sec.dat","?",num2str(t_length)));

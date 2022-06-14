fs = 2e5;
fc=20;

amp = 0.3;

t_length = 1;
t = [0:1/fs:t_length-1/fs];

y = 0*amp*sin(2*pi*fc*t);

write_complex_binary(y, strrep("../../tx_outputs/zero_?sec.dat","?",num2str(t_length)));

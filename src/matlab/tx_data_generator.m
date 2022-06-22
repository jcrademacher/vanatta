%preamble = [0 0 0 1 1 1 0 1 1 0 1];
preamble = [0 0 1 1 1 0 1 0];
%preamble_fm0 = [0 1 0 1 0 0 1 1 0 0 1 0 1 1 0 1];

N_preamble_bits = length(preamble);
N_data_bits = 12;
N_packet_bits = N_data_bits+N_preamble_bits;
prbs_order = 15;

fc = 20e3;
fs = 2e5;
fb = 1e3;
fb_preamble = 100;

fm0_samp = fs/fb;
fm0_samp_preamble = fs/fb_preamble;

N_packets = 1000;
[next_data,seed] = prbs(prbs_order,N_data_bits);

time_to_tx = N_packet_bits*N_packets/fb;
baseband = [];

all_data = [];

[preamble_fm0,state] = generate_fm0_sig2(preamble,fm0_samp_preamble);
[next_data_fm0,state] = generate_fm0_sig2(next_data,fm0_samp,state);
next_packet_fm0 = [preamble_fm0 next_data_fm0];

for pack=1:N_packets
    baseband = [baseband next_packet_fm0];
    all_data = [all_data next_data];
    
    [next_data,seed] = prbs(prbs_order,N_data_bits,seed);
    [preamble_fm0,state] = generate_fm0_sig2(preamble,fm0_samp_preamble,state);
    [next_data_fm0,state] = generate_fm0_sig2(next_data,fm0_samp,state);
    next_packet_fm0 = [preamble_fm0 next_data_fm0];
end

amp = 1/3;
m=0.5;

t = [0:(1/fs):(length(baseband)-1)/fs];

encoded_baseband = amp*(1+m.*baseband);
encoded_baseband_preamble = encoded_baseband(1:N_preamble_bits*fs/fb_preamble);

tx = encoded_baseband.*cos(2*pi*fc*t);
%tx = baseband_long.*cos(2*pi*fc*t_long);
plot(t,encoded_baseband);
%ylim([-1.5 1.5]);

modulated_filename = strcat("../../tx_outputs/tx_modulated_fm0_",strrep(num2str(fc/1e3),'.',','),"kfc_",...
                            strrep(num2str(fb/1e3),'.',','),"kbps_m=",...
                            strrep(num2str(m,2),'.',','),"_a=",...
                            strrep(num2str(amp,2),'.',','),".dat");

data_filename = strcat("../../tx_outputs/data_prbs_order=",strrep(num2str(prbs_order),'.',','),"_len=",...
                            strrep(num2str(N_data_bits),'.',','),"_packets=",...
                            strrep(num2str(N_packets),'.',','),".dat"); 

preamble_filename = strcat("../../tx_outputs/baseband_preamble_fm0_",num2str(N_preamble_bits),'bit_',...
                            strrep(num2str(fb_preamble/1e3),'.',','),"kbps_m=",...
                            strrep(num2str(m,2),'.',','),"_a=",...
                            strrep(num2str(amp,2),'.',','),".dat"); 

write_complex_binary(tx,modulated_filename);
write_complex_binary(all_data,data_filename);
write_complex_binary(encoded_baseband_preamble,preamble_filename);
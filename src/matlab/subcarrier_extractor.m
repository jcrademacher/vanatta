
root = '../../rx_outputs/Sea Grant PAB Van Atta 05-05-2022/rx_sg_+20cm_backscatter_pab2_ind_tmux_22kfc_1kHz_square_0,5m_depth_30cm_dis_20cm_hphydro';

Ntrials = 2;
fmod = 1e3;

if Ntrials > 1
    for fnum=0:Ntrials-1
        sig = read_complex_binary(strcat(root,'_',int2str(fnum),'.dat'));
        
        rx_signals(:,fnum+1) = real(sig);
    end
else
    sig = read_complex_binary(strcat(root,'_0.dat'));
    rx_signals(:,1) = real(sig);
end

rx_len = length(rx_signals(:,1));

fs = 2e5;
%f = linspace(1e3,100e3,1000);
[pxx,w] = periodogram(rx_signals,hamming(rx_len),fs);
f = w/(2*pi)*fs;

[carrier_pow,carrier_index] = max(pxx);
subcarrier_left_index = carrier_index - fmod;
subcarrier_right_index = carrier_index + fmod;

subcarrier_left_window = linspace(subcarrier_left_index(1) - fmod/2,subcarrier_left_index(1) + fmod/2,fmod+1);
subcarrier_right_window = linspace(subcarrier_right_index(1) - fmod/2,subcarrier_right_index(1) + fmod/2,fmod+1);

subcarrier_left_max = max(pxx(subcarrier_left_window,:));
subcarrier_right_max = max(pxx(subcarrier_right_window,:));

avg_subcarrier_pow = mean(cat(2,subcarrier_left_max,subcarrier_right_max));

disp("Average subcarrier power (dB): ");
disp(10*log10(avg_subcarrier_pow));

% carrier_freq = f(carrier_index)';
% subcarrier_left_freq = carrier_freq - fmod;
% subcarrier_right_freq = carrier_freq + fmod;

plot(f,10*log10(pxx(:,1)));


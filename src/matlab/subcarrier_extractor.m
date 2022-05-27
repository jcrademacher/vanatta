% place question mark where degree # should be
folder = '~/Documents/sk/oceans/vanatta/rx_outputs/River PAB Van Atta 05-26-2022/';
file = 'rx_river_backscatter_pab14_ind_array_?deg_tmux_22kfc_1kHz_square_1m_depth_5,8m_dis_4,8m_hphydro';

root = strcat(folder,file);

Ntrials = 3;
fmod = 1e3;

degree_list = linspace(-90,90,13);
Ndeg = length(degree_list);

measured_pattern = zeros(Ndeg,2);
measured_pattern(:,1) = degree_list * pi/180;

for n=1:Ndeg
    filename = strrep(root,"?",num2str(degree_list(n)));

    if Ntrials > 1
        for fnum=0:Ntrials-1
            sig = read_complex_binary(strcat(filename,'_',int2str(fnum),'.dat'));
            
            rx_signals(:,fnum+1) = real(sig);
        end
    else
        sig = read_complex_binary(strcat(filename,'_0.dat'));
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
    
    measured_pattern(n,2) = 10*log10(avg_subcarrier_pow);
%     disp("Average subcarrier power (dB): ");
%     disp(10*log10(avg_subcarrier_pow));
end

polarplot(measured_pattern(:,1),measured_pattern(:,2));
hold on;
rlim([-60 -50]);


%plot(f,10*log10(pxx(:,1)));


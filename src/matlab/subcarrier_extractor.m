clear pxx;

% place question mark where degree # should be
folder = '~/Documents/MIT/sk/oceans/vanatta/rx_outputs/River PAB Van Atta 06-01-2022/';
file = ['rx_river_backscatter_pab_007A_007B_ind_array_?deg_tmux_18,5kfc_1kHz_square_1m_depth_5,8m_dis_4,8m_hphydro_120sec'];

root = strcat(folder,file);

Ntrials = 1;
fmod = 1e3;

degree_list = [-30];
Ndeg = length(degree_list);

measured_pattern = zeros(Ndeg,2);
measured_pattern(:,1) = degree_list * pi/180;

avg_subcar_pow_trial = zeros(Ndeg,Ntrials);

%figure(1);
figure(2);
lg = legend;

for n=1:Ndeg
    filename = strrep(root,"?",num2str(degree_list(n)));
    
    sig = read_complex_binary(strcat(filename,'_0','.dat'));
    rx_signals = real(sig);
    sig_len = length(sig);

    trial_length = floor(sig_len/Ntrials);

    for i=1:Ntrials
       % segment = rx_signals()
    end

%     if Ntrials > 1
%         for fnum=0:Ntrials-1
%             sig = read_complex_binary(strcat(filename,'_',int2str(fnum),'.dat'));
%             
%             rx_signals(:,fnum+1) = real(sig);
%         end
%     else
%         sig = read_complex_binary(strcat(filename,'_0.dat'));
%         rx_signals(:,1) = real(sig);
%     end
%     
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
%     avg_subcar_pow_trial(n,:) = mean(cat(1,subcarrier_left_max,subcarrier_right_max),1);
    
%     measured_pattern(n,2) = 10*log10(avg_subcarrier_pow);
    disp("Average subcarrier power (dB): ");
    disp(10*log10(avg_subcarrier_pow));
    
%     [f,x] = ecdf(10*log10(avg_subcar_pow_trial(n,:)));
%     figure(2);
%     hold on;
%     plot(x,f);
%     lg.String{n} = strcat(num2str(degree_list(n))," deg");
end

% figure(1);
% polarplot(measured_pattern(:,1),measured_pattern(:,2));
% hold on;
% rlim([-80 -50]);
% 
% figure(2);
% xlabel("Subcarrier Power (dB)");
% ylabel("F(x)");



plot(f,10*log10(pxx(:,1)));


clear pxx;

% place question mark where degree # should be
folder = '~/Documents/sk/oceans/vanatta/rx_outputs/River PAB Van Atta 06-01-2022/';
file = ['rx_river_backscatter_pab_007A_007B_ind_vanatta_?deg_tmux_18,5kfc_1kHz_square_1m_depth_5,8m_dis_4,8m_hphydro_120sec'];

root = strcat(folder,file);

fmod = 1e3;
fs = 2e5;
Nsamps = 120*fs;

degree_list = [-30 -15 0 15 30];
Ndeg = length(degree_list);

measured_pattern = zeros(Ndeg,2);
measured_pattern(:,1) = degree_list * pi/180;

trial_length = 100000;
Ntrials = floor(Nsamps/trial_length);

avg_subcarrier_pow = zeros(Ntrials,Ndeg);
rx_signals = zeros(Nsamps,Ndeg);

for n=1:Ndeg
    filename = strrep(root,"?",num2str(degree_list(n)));
    
    sig = read_complex_binary(strcat(filename,'_0','.dat'));
    rx_signals(:,n) = real(sig(1:Nsamps));
end

for i=1:Ntrials
    segment = rx_signals((i-1)*trial_length+1:i*trial_length,:);
%     [pxx,w] = periodogram(segment,hamming(trial_length),fs);
    Nfft = trial_length;
    seg_fft = fft(segment,Nfft);
    pxx = abs(seg_fft(1:Nfft/2,:)).^2;
    f = 0:fs/Nfft:fs/2-fs/Nfft;
    
    [carrier_pow,carrier_index] = max(pxx);
    subcarrier_left_index = carrier_index - fmod*Nfft/fs;
    subcarrier_right_index = carrier_index + fmod*Nfft/fs;
    
    subcarrier_left_window = linspace(subcarrier_left_index(1) - fmod/2*Nfft/fs,subcarrier_left_index(1) + fmod/2*Nfft/fs,fmod*Nfft/fs+1);
    subcarrier_right_window = linspace(subcarrier_right_index(1) - fmod/2*Nfft/fs,subcarrier_right_index(1) + fmod/2*Nfft/fs,fmod*Nfft/fs+1);
    
    subcarrier_left_max = max(pxx(subcarrier_left_window,:),[],1);
    subcarrier_right_max = max(pxx(subcarrier_right_window,:),[],1);
    
    avg_subcarrier_pow(i,:) = mean(cat(1,subcarrier_left_max,subcarrier_right_max));

%     plot(f,pxx);
end

for n=1:Ndeg
    [cdf,x] = ecdf(10*log10(avg_subcarrier_pow(:,n)));
    figure(n);
    plot(x,cdf);
    title(strcat("ECDF of ",num2str(degree_list(n)),' deg'));
    grid on;

    disp(strcat("Avg Pow @ ",num2str(degree_list(n))," deg (dB): "));
    disp(num2str(mean(10*log10(avg_subcarrier_pow(:,n)))));
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
   
%     avg_subcar_pow_trial(n,:) = mean(cat(1,subcarrier_left_max,subcarrier_right_max),1);
    
%     measured_pattern(n,2) = 10*log10(avg_subcarrier_pow);
%     disp("Average subcarrier power (dB): ");
%     disp(10*log10(avg_subcarrier_pow));
    



% figure(1);
% polarplot(measured_pattern(:,1),measured_pattern(:,2));
% hold on;
% rlim([-80 -50]);
% 
% figure(2);
% xlabel("Subcarrier Power (dB)");
% ylabel("F(x)");



% plot(f,10*log10(pxx(:,1)));


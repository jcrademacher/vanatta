fs = 2e5;
fc = 18.5e3;
fb = 1000;
c = 1500;
wc = 2*pi*fc;

init_delay = 50e-3;

dd = 3;
du = 2;
d0 = dd-du;

% A0 is direct path mag, Ad is downlink path mag, Au is uplink path mag
A0 = 0.9;
Ad = 0.1;
Au = 0.2;

tau_0 = d0 / c;
tau_d = dd / c;
tau_u = du / c;

Nel = 1;

preamble1 = [0 0 1 1 1 0 1 0];
expected_preamble1 = generate_fm0_sig2(preamble1,fs/fb);

preamble2 = preamble1;
expected_preamble2 = generate_fm0_sig2(preamble2,fs/fb);

expected_preamble_comb = expected_preamble1;

% m = 0.1;
t_len = 1;
t = [0:1/fs:t_len-1/fs];

%plot(t,data(t-0.1,fb,1));
% n_data_reps1 = 5;
% n_data_reps2 = 5;
% data1 = data(t-tau_u,preamble1,fb,n_data_reps1);
% data2 = data(t-(length(expected_preamble1)*n_data_reps1/fs+1e-3)-tau_u,preamble2,fb,n_data_reps2);
% 
% data_comb_part1 = data(t-(length(expected_preamble1)*n_data_reps1/fs+length(expected_preamble2)*n_data_reps2/fs+2e-3)-tau_u,...
%                     preamble1,fb,n_data_reps1);
% data_comb_part2 = data(t-(length(expected_preamble1)*n_data_reps1/fs+length(expected_preamble2)*n_data_reps2/fs+2e-3)-tau_u,...
%                     preamble2,fb,n_data_reps2);
% 
% data_comb = data_comb_part1+data_comb_part2;
% 
% %yr = A0*cos(wc*(t-tau_0)).*heaviside(t-tau_0) + ... 
%         Ad*Au*data1.*cos(wc*(t-tau_d-tau_u)).*heaviside(t-tau_d-tau_u) + ...
%         Ad*Au*data2.*cos(wc*(t-tau_d-tau_u)).*heaviside(t-tau_d-tau_u) + ...
%         Ad*Au*data_comb.*cos(wc*(t-tau_d-tau_u)).*heaviside(t-tau_d-tau_u);

yr = read_complex_binary('../../rx_outputs/River PAB Channel Estimate 07-15-2022/rx_array_chest_pab_007B_005A_ind_+63,75deg_tmux_18,5kfc_siggen_data_1kbps_usrp_3m_depth_2m_u2b_1m_hphydro_0.dat');        
%plot(t,yr);
% carrier = 0;
% pb_sig = (1+m*data).*carrier;
% pb_sig = pb_sig';

% load 3m_pb_ch.mat;
% channel = ch_pb2;
% figure;
% periodogram(channel);

% out = conv(pb_sig,channel);
% 
% t_out = [0:1/fs:length(out)/fs-1/fs];
% downconv_out = fftfilt(lpFilt, cos(2*pi*fc*t_out)'.*out);

% figure;
% subplot(2,1,1);
% plot(pb_sig);
% ylim([-2 2]);
% subplot(2,1,2);
% plot(downconv_out);

% t_window = 0.1;
% window = chebwin(t_window*fs);
% Nfft = 2^nextpow2(length(window));

% figure;
% hold on;
% [pxx,f] = pwelch(data,window,[],Nfft,fs,'power');
% plot(f,10*log10(pxx));
% grid on;
% grid minor;

% generate low and highpass filters
% lowpass is performed before downsampling as an anti-aliasing filter
% highpass is performed after downsampling as it is a higher order filter
fsb1 = fb/100;
fpb1 = fb/2;
dfac = 1;   % donwsampling factor

fpb1_lp = 7e3;
fsb1_lp = 9e3;

% % highpass for after downsampling
hpFilt = designfilt('highpassfir','PassbandFrequency',fpb1*2/(fs/dfac) ...
                    ,'StopbandFrequency',fsb1*2/(fs/dfac),'StopbandAttenuation',80,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

% lowpass after downconversion, before downsampling
lpFilt = designfilt('lowpassfir' ...
                    ,'PassbandFrequency',fpb1_lp*2/fs...
                    ,'StopbandFrequency',fsb1_lp*2/fs,'StopbandAttenuation' ...
                    ,80,'PassbandRipple',0.1,'DesignMethod','kaiserwin');

%%%% END DESIGN PARAMETERS %%%%

sig = yr(24:end);

rx_len = length(sig);
% Nel x rx_len size matrix of input signals, where each row is time-series on an individual array element
rx_signals = zeros(Nel,rx_len);
rx_signals(1,:) = real(sig);

%%%% CARRIER FREQUENCY AND PHASE EXTRACTION %%%%
% have had some issues with it in the past and since RX and TX USRPs are 
% synchronized in most experiments directly using the known fc works fine

Nfft = 100*fs;
rx_fft = fft(rx_signals',Nfft)';
fft_mag = abs(rx_fft);
max_search = [round(Nfft/fs*(fc-1)):round(Nfft/fs*(fc+1))];
[maxval,mindex] = max(fft_mag(:,max_search),[],2); % max in each row
carrier_phase = angle(rx_fft(max_search(mindex)'));
carrier_freq = fs/Nfft*max_search(mindex)';

carrier_freq = fc;

 % generate the time series and local oscillator
t = [0:1/fs:(rx_len-1)/fs];
lo = exp(1j*(2*pi*carrier_freq*t));
% downconvert
rx_baseband = rx_signals.*lo;
% slice out where data starts
%rx_baseband = rx_baseband(1500:end);

[gdlp,w] = grpdelay(lpFilt);
gdlp = mean(gdlp);

[gdhp,w] = grpdelay(hpFilt);
gdhp = mean(gdhp);

% lowpass filtering both removes the 2fc term and anti-alias filters
% filtfilt used for 0 group delay filtering
rx_baseband = fftfilt(lpFilt,rx_baseband')';
% rx_baseband = rx_baseband(gdlp+1:end);

expected_preamble1 = filtfilt(lpFilt,expected_preamble1')';
expected_preamble2 = filtfilt(lpFilt,expected_preamble2')';
expected_preamble_comb = filtfilt(lpFilt,expected_preamble_comb')';

% remove DC mean
% rx_baseband = rx_baseband - mean(rx_baseband);
rx_baseband = fftfilt(hpFilt,rx_baseband')';
% rx_baseband = rx_baseband(gdhp+1:end);
%expected_preamble1 = fftfilt(hpFilt,expected_preamble1);

rx_baseband = rx_baseband(init_delay*fs+gdlp+gdhp-40:end);

sig_sec = rx_baseband;
Nfft = 2^nextpow2(length(sig_sec));
fft_sig = fft(sig_sec,Nfft);
[subcar_peak,subcar_peak_index] = max(fft_sig);
subcar_peak_f = subcar_peak_index/Nfft*fs;
subcar_peak_phase = angle(subcar_peak);
f = [-fs/2:fs/Nfft:fs/2-fs/Nfft];

t_window = 0.1;
window = chebwin(t_window*fs);
%Nfft = length(window);

figure(1);
hold on;
%[pxx,f] = pwelch(sig_sec,window,[],Nfft,fs,'power');
plot(f,20*log10(abs(fftshift(fft_sig))));
grid on;
grid minor;

% disp("Total Average Power (dB):");
% disp(num2str(10*log10(tot_pow)));

figure(2);
subplot(2,1,1);
plot(real(sig_sec));
subplot(2,1,2);
plot(imag(sig_sec));


%% CORRELATE, ESTIMATE, DECODE, and COMPUTE BER %%

% N_data_bits = 0;
% 
% % resample for even integer fm0_samp value
% % with fs/fb divisble by 2 this code is not run
% if mod(fs/fb,2) ~= 0
%     p = (ceil(fs/fb)+1)*fb;
%     q = fs;
%     
%     rx_baseband = resample(rx_baseband,p,q);
%     
%     fs = p;
% end

fm0_samp = fs/fb;

% DO NOT CHANGE %
% data_len = round(N_data_bits/fb*fs);                % known data length in samples
preamble_len1 = round(length(expected_preamble1));    % known preamble length in samples
preamble_len2 = round(length(expected_preamble2));
packet_len = preamble_len1;                 % known packet length in samples

% number of packets to decode
% usually this needs to be set slightly lower than the total that was
% actually transmitted, I haven't investigated this further and am not sure
% why this is needed
N_trials = 2;
N_packets1 = 5;
N_packets2 = 5;
N_comb_packets = 5;

node_delay = 1e-3;

% % estimates after channel projection for individual elements
% rx_rep_estimates = zeros(Nel,N_packets*(packet_len-preamble_len1));

channel_estimates = zeros(Nel,N_trials*(N_packets1+N_packets2+N_comb_packets));

A_hat = zeros(Nel,N_trials*(N_packets1+N_packets2+N_comb_packets));
ang_hat = zeros(Nel,N_trials*(N_packets1+N_packets2+N_comb_packets));

preamble_starts = zeros(Nel,1);

decode_preamble1 = expected_preamble1-mean(expected_preamble1);
decode_preamble2 = expected_preamble2-mean(expected_preamble2);
decode_preamble_comb = expected_preamble_comb-mean(expected_preamble_comb);

figure;
% CORRELATION AND DECODING %
for pnum=1:N_trials*(N_packets1+N_packets2+N_comb_packets)
    for el=1:Nel
        begdex = (pnum-1)*packet_len+1;
        endex = pnum*packet_len+40*(pnum-1);
        end_preamble_dex = (pnum-1)*packet_len+1+preamble_len1;

%         beg_data_dex = (pnum-1)*data_len+1;
%         end_data_dex = pnum*data_len;
    
%         beg_bit_dex = (pnum-1)*N_data_bits+1;
%         end_bit_dex = pnum*N_data_bits;
        
        if pnum <= N_packets1 || (pnum > N_packets1+N_packets2+N_comb_packets && pnum <= 2*N_packets1+N_packets2+N_comb_packets)
            decode_preamble = decode_preamble1;
        elseif pnum <= N_packets1+N_packets2 || (pnum > 2*N_packets1+N_packets2+N_comb_packets && pnum <= 2*N_packets1+2*N_packets2+N_comb_packets)
            decode_preamble = decode_preamble2;
        else
            decode_preamble = decode_preamble_comb;
        end
        
        % tx norm is length of preamble for binary keying (-1,+1)
        tx_norm = sum(abs(decode_preamble).^2);
        
        if pnum == 1
            % perform cross correlation for packet start
            [rcorr,rlags] = xcorr(real(rx_baseband(el,begdex:end_preamble_dex))',decode_preamble');
            [icorr,ilags] = xcorr(imag(rx_baseband(el,begdex:end_preamble_dex))',decode_preamble');
            % removes tails of correlation 
            corr_tot = rcorr(end_preamble_dex-begdex+1:end)+1j*icorr(end_preamble_dex-begdex+1:end);
            abs_corr = abs(corr_tot);
            % find maximum correlation and begin decoding from there
            [preamble_max,preamble_starts(el)] = max(abs_corr); 
        end
        
        begdex = begdex + floor((pnum-1) / N_packets1)*node_delay*fs;
        endex = endex + floor((pnum-1) / N_packets1)*node_delay*fs;

        clf;
        
        plot(abs_corr/30);
        hold on;
        plot(real(rx_baseband(el,begdex:end_preamble_dex)));
        plot([zeros(1,preamble_starts(el)) decode_preamble/500]);
        %xlim([0 1000]);

        % slice out packet found from correlation
        packet = rx_baseband(el,begdex+preamble_starts(el)-1:endex+preamble_starts(el)-1);
        % remove the mean
        packet = packet - mean(packet(1:preamble_len1));
        % slice out preamble
        packet_preamble = packet(1:preamble_len1);
        % slice out data
        packet_data = packet(preamble_len1+1:end);
        
        %remove edges
        packet_preamble = remove_edges(packet_preamble,preamble1,0.05,fb,fs);
        decode_preamble = remove_edges(decode_preamble,preamble1,0.05,fb,fs);
        % recompute tx_norm since some samples have been removed
        tx_norm = sum(abs(decode_preamble).^2);

        % compute the channel estimate 
        channel_estimates(el,pnum) = sum(packet_preamble.*conj(decode_preamble))/tx_norm;
        A_hat(el,pnum) = 2*abs(channel_estimates(el,pnum));
        ang_hat(el,pnum) = angle(channel_estimates(el,pnum));
        
        % extract estimate of data only
%         rx_rep_estimates(el,beg_data_dex:end_data_dex) = packet_data.*conj(channel_estimates(el,pnum))./abs(channel_estimates(el,pnum)).^2;
%         packet_preamble_proj = packet_preamble.*conj(channel_estimates(el,pnum));
%         preamble_estimate = packet_preamble_proj./abs(channel_estimates(el,pnum)).^2;
%         summed_preamble_estimate(pnum,:) = summed_preamble_estimate(pnum,:) + packet_preamble_proj;
%         % SNR calculation (currently incorrect)
%         noise_est = decode_preamble-preamble_estimate;
%         rx_snr_packet(el,pnum) = mean(abs(preamble_estimate).^2)/mean(abs(noise_est).^2);

%         plot(decode_preamble);
%         hold on;
%         plot(real(preamble_estimate));
        
%         % decode packet
%         if (mod(fm0_samp,2) == 0)
%             % should always enter this branch, camera_decode written by
%             % saad and waleed
%             bits = camera_decode(rx_rep_estimates(el,beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%         elseif mod(floor(fm0_samp),2)==0
%             bits = fm0_decode_new_R12_dec_even(rx_rep_estimates(el,beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%         else
%             bits = fm0_decode_new_R12(rx_rep_estimates(el,beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%         end
%         
%         % place bits into decoded data matrix
%         rx_decoded_data(el,beg_bit_dex:end_bit_dex) = bits;
    end
    
%     % sum array elements 
%     rx_rep_est_summed(beg_data_dex:end_data_dex) = sum(rx_rep_estimates(:,beg_data_dex:end_data_dex).*abs(channel_estimates(:,pnum)).^2,1)./sum(abs(channel_estimates(:,pnum)).^2,1);
%     summed_preamble_estimate(pnum,:) = summed_preamble_estimate(pnum,:) / sum(abs(channel_estimates(:,pnum)).^2,1);
%     summed_noise_est = decode_preamble-summed_preamble_estimate(pnum,:);
%     rx_snr_packet_summed(pnum) = mean(abs(summed_preamble_estimate(pnum,:)).^2)/mean(abs(summed_noise_est).^2);
% 
%     % decode summed packet
%     if (mod(fm0_samp,2) == 0)
%         bits = camera_decode(rx_rep_est_summed(beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%     elseif mod(floor(fm0_samp),2)==0
%         bits = fm0_decode_new_R12_dec_even(rx_rep_est_summed(beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%     else
%         bits = fm0_decode_new_R12(rx_rep_est_summed(beg_data_dex:end_data_dex)',fm0_samp,N_data_bits);
%     end
%     
%     % place summed bits into decoded data matrix
%     rx_decoded_summed_data(beg_bit_dex:end_bit_dex) = bits;
end

mag_percent_error = abs(A_hat - Au*Ad) / (Au*Ad) * 100;
ang_percent_error = abs(ang_hat - angle(exp(1j*2*pi*fc*(tau_u+tau_d)))) / angle(exp(1j*2*pi*fc*(tau_u+tau_d))) * 100;

n1_ch_est = [channel_estimates(1:N_packets1) channel_estimates(N_packets1+N_packets2+N_comb_packets+1:2*N_packets1+N_packets2+N_comb_packets)];
n2_ch_est = [channel_estimates(N_packets1+1:N_packets1+N_packets2) channel_estimates(2*N_packets1+N_packets2+N_comb_packets+1:2*N_packets1+2*N_packets2+N_comb_packets)];
comb_ch_est = [channel_estimates(N_packets1+N_packets2+1:N_packets1+N_packets2+N_comb_packets) channel_estimates(2*N_packets1+2*N_packets2+N_comb_packets+1:end)];

figure;
hold on;
plot(n1_ch_est(1:N_packets1),'o','linewidth',2,'Color','r');
plot(n2_ch_est(1:N_packets2),'o','linewidth',2,'Color','g');
plot(comb_ch_est(1:N_comb_packets),'o','linewidth',2,'Color','b');

plot(n1_ch_est(N_packets1+1:end),'x','linewidth',2,'Color','r');
plot(n2_ch_est(N_packets2+1:end),'x','linewidth',2,'Color','g');
plot(comb_ch_est(N_comb_packets+1:end),'x','linewidth',2,'Color','b');

xlim([-1 1]*1e-3);
ylim([-1 1]*1e-3);

xlabel("Real");
ylabel("Imaginary");

%axis equal;
grid on;
grid minor;
% 
disp("Mean Sum of Single Node Channels Estimates");
h1_h2 = mean(n1_ch_est+n2_ch_est)
disp("Variance Sum of Single Node Channel Estimates");
var_sum_ch_est = var(n1_ch_est+n2_ch_est)

disp("Mean of Combined Channel Estimates");
h_tot = mean(comb_ch_est)
disp("Variance of Combined Channel Estimates");
var_comb_ch_est = var(comb_ch_est)

snr_sum_ch_est1 = abs(mean(n1_ch_est(1:N_packets1)+n2_ch_est(1:N_packets2)))^2/var(n1_ch_est(1:N_packets1)+n2_ch_est(1:N_packets2));
snr_sum_ch_est2 = abs(mean(n1_ch_est(N_packets1+1:end)+n2_ch_est(N_packets2+1:end)))^2/var(n1_ch_est(N_packets1+1:end)+n2_ch_est(N_packets2+1:end));

snr_comb_ch_est1 = abs(mean(comb_ch_est(1:N_comb_packets)))^2/var(comb_ch_est(1:N_comb_packets));
snr_comb_ch_est2 = abs(mean(comb_ch_est(N_comb_packets+1:end)))^2/var(comb_ch_est(N_comb_packets+1:end));

snr_sum_ch_est = 10*log10((snr_sum_ch_est1+snr_sum_ch_est2)/2);
snr_comb_ch_est = 10*log10((snr_comb_ch_est1+snr_comb_ch_est2)/2);

disp("Error")
(h1_h2-h_tot)./abs(h_tot) * 100

disp("Mean Phase Difference")
mean(angle(n1_ch_est./n2_ch_est)/pi*180)
% angle(comb_ch_est)/pi*180

disp("Expected Phase Difference")
angle(exp(1j*(2*2*pi*0.07*sin(-60/180*pi))/(1500/fc)))/pi*180

% (abs(comb_ch_est(2:end))-abs(n1_ch_est+n2_ch_est(2:end)))./abs(comb_ch_est(2:end))

%% EXPORT DATA %%
% t0 = t - init_delay;
% 
% len_packet1 = length(expected_preamble1)*n_data_reps1/fs;
% len_packet2 = length(expected_preamble2)*n_data_reps2/fs;
% 
% data_ch1_1 = data(t0,preamble1,fb,n_data_reps1);
% data_ch2_1 = data(t0-len_packet1-1e-3,preamble2,fb,n_data_reps2);
% data_comb_1 = data(t0-len_packet1-len_packet2-2e-3,preamble2,fb,n_data_reps2);
% 
% data_ch1_2 = data(t0-len_packet1-2*len_packet2-3e-3,preamble1,fb,n_data_reps1);
% data_ch2_2 = data(t0-2*len_packet1-2*len_packet2-4e-3,preamble2,fb,n_data_reps2);
% data_comb_2 = data(t0-2*len_packet1-3*len_packet2-5e-3,preamble2,fb,n_data_reps2);
% 
% data_tot = (data_ch1_1+data_ch1_2+data_comb_1+data_comb_2+1)/2.5+1j*(data_ch2_1+data_ch2_2+data_comb_1+data_comb_2+1)/2.5;
% 
% figure;
% hold on;
% plot(real(data_tot));
% plot(imag(data_tot));
% 
% write_complex_binary(data_tot,"../../tx_outputs/array_channel_estimating_data.dat");

% out = remove_edges(expected_preamble1,preamble1,0.05,fb,fs);
% plot(out);
% hold on;
% plot(expected_preamble1);

function out = data(t,code,fb,nreps)
    fm0_code = generate_fm0_sig2(code,2);
    code_len = length(fm0_code);
    out = zeros(size(t),'like',t);

    for n=1:nreps*code_len
        out(logical((t < n*1/(2*fb)).*(t >= (n-1)*1/(2*fb)))) = fm0_code(mod(n-1,code_len)+1);
    end
end

% removal is a percentage of the bit period
function out = remove_edges(data,expected_code,percent_removal,fb,fs)
    if length(data) ~= length(expected_code)*fs/fb
        error("Data length is not equal to expected length based on fb and fs given")
    end

%     fm0_code = generate_fm0_sig2(expected_code,2);
    
    % array of how many expected edges there are and their indexes (without
    % beginning edge)
    edge_indexes = zeros(length(expected_code)+length(expected_code(expected_code==0)),1);
    
    nedge = 1;
    for n=1:length(expected_code)
        % fill edge center indexes
        bit = expected_code(n);
        
        if bit
            edge_indexes(nedge) = n*fs/fb;
            nedge = nedge + 1;
        else
            edge_indexes(nedge) = n*fs/fb-fs/(2*fb);
            nedge = nedge + 1;
            edge_indexes(nedge) = n*fs/fb;
            nedge = nedge + 1;
        end
    end

    sample_expansion = 2*floor(fs/fb*percent_removal/2)+1;
    samples_to_remove = repelem(edge_indexes,sample_expansion);
    samples_to_remove = samples_to_remove + repmat([-(sample_expansion-1)/2:1:(sample_expansion-1)/2]',length(edge_indexes),1);
    samples_to_remove(end-(sample_expansion-1)/2+1:end) = [];
    samples_to_remove = cat(1,[1:(sample_expansion-1)/2]',samples_to_remove);

    data(samples_to_remove) = [];
    out = data;
end

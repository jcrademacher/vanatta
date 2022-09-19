function [weights,ber_fin_b,ber_fin_a,snr_final] = DFE_500_vanatta(y_rx,complete_bits,num_pkt_ind,num_testing_pkts,const_lvl,W,AAT)
% y_rx = raw data with first sample being first sample of first bit (column
% vector)
% complete_bits = all bits in raw form (0,1) sequence, each row is packet
% including preamble followed by data
% num_pkt_ind = 100 (training packets)
% num_testing_packets = testing packets, test on remaining packets (faster
% if lower)
% const_lvl = 1.126e-4, level of bits in baseband
% W = 0, initially, then you can put weights after return value is given,
% may provide improvement
% AAT = adapt after training, keep as 1

% OUTPUTS
% weights = DFE weights
% ber_fin_b = BER before DFE
% ber_fin_a = BER after DFE
% snr_final = SNR on baseband

% if you want to use the weights upon a second pass, reduce # of training
% packets so weights aren't changed too much

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%num_pkt_ind =1;
%num_testing_pkts = 600;
Fs = 2e5;
%W = 0;
%W = weights;
%load jack_data_vanatta;
ds_ind = 1;%1;%5;
y_rx = downsample(y_rx,ds_ind);
%y_rx = y_rx(1:end-0.1*Fs);
f_factor = 0.99;%[0.98 0.983 0.989 0.999 1];
%num_pkt_ind = 1;%[90:1:105];
ref_tap = 1;ff_tap = 203;%53;%23; 
fb_tap =100;%200;%100;
%const_lvl =1.126e-4;%max(real(y_rx(3.84e6:3.89e6)));%(1.0350e-04); %(3.126e-4);%3.126e-4; 
%num_pkt_ind = 28*ones(1,20);
Fs = 2e5/ds_ind;
data_rate = 500;%5000;%2e3;
fm0samp = Fs/data_rate;
pkt_bits = 24; %50 for waleed 24 for jack
pre_bits = 8; % 10 for waleed 8 for jack
pkt_size = pkt_bits*fm0samp;



row = 1;
for ff = f_factor
    %clearvars -except ff num_pkt row num_pkt_ind f_factor y_rx complete_bits ber_vec_a ber_vec_b cnttt row 
    
    chunk = 1;
    %num_pkt = 5;
    off = 1;
    temp = chunk;
    training_pkt = [];
    training_bits = [];
    
    cnttt = 1;
    
    for num_pkt = num_pkt_ind
        training_pkt = [];
        training_bits = [];
        temp = chunk;
        for ss = 1:num_pkt
%             if(num_testing_pkts == 26)
%                 disp('ola')
%             end
            training_bits = [training_bits complete_bits(temp,:)];
            training_pkt = [training_pkt generate_fm0_sig2([complete_bits(temp,:)],fm0samp).*const_lvl]; 
            temp = temp + 1;
        end
        
        train_symb = downsample(training_pkt,fm0samp/2);

        rx_data = [y_rx((chunk-off)*pkt_size+1:(chunk-off)*pkt_size+num_pkt*pkt_size+num_testing_pkts*pkt_size)];

         acorr = 0.1;
         eqdfe_lms = comm.DecisionFeedbackEqualizer('Algorithm','RLS','InputSamplesPerSymbol',fm0samp/2, 'ReferenceTap',ref_tap, ...
             'ForgettingFactor',ff,'NumForwardTaps',ff_tap,'AdaptAfterTraining',AAT,'NumFeedbackTaps',fb_tap,...
             'Constellation',[min(training_pkt) max(training_pkt)],'InitialInverseCorrelationMatrix',acorr,'InitialWeights',(W),'InitialWeightsSource','Property');
 
         %mxStep = maxstep(eqdfe_lms,rx_data);
         [eq_y err weights] = eqdfe_lms((rx_data),real(train_symb).');
         
        %[eq_y err weights] = dfe_eq_func(ff,fm0samp,(rx_data),real(train_symb),training_pkt);

        %plot(eq_y(length(train_symb):length(train_symb)+100));

        release(eqdfe_lms);
        pkt_len = pkt_bits*2;

        tr_indx = 1:num_pkt*pkt_len;
        te_indx = tr_indx(end)+1:tr_indx(end)+num_testing_pkts*pkt_len;

        tr_ind = (chunk-off)*pkt_size+1:(chunk-off)*pkt_size+num_pkt*pkt_size;
        te_ind = tr_ind(end)+1:tr_ind(end)+pkt_size*num_testing_pkts;
        %figure;
        bit_vec_tx = reshape(complete_bits(chunk+num_pkt:chunk+num_pkt+ ...
            num_testing_pkts-1,:).',1,prod(size(complete_bits(chunk+ ...
            num_pkt:chunk+num_pkt+num_testing_pkts-1,:))));

        act_data = generate_fm0_sig2(bit_vec_tx,fm0samp).*1e-3;

        %% BER/SNR for test equalizer data
        tot_bits_vec = [];
        bits_dat = eq_y(te_indx);
        %rx_data_dec = rx_data(te_ind);
        rx_data_dec = y_rx(te_ind);
        for y_bit = 1:length(eq_y(te_indx))
           if (bits_dat(y_bit)) < 0
             tot_bits_vec = [tot_bits_vec -0.5.*ones(1,fm0samp/2)]; 
          else
             tot_bits_vec = [tot_bits_vec 0.5.*ones(1,fm0samp/2)];
          end
        end


         tot_bits_a = [];
         tot_bits_b = [];
         act_bits_a = [];
         act_bits_b = [];
         snr_vec =[];
         cnt = 1;
         for kk = 1:num_testing_pkts
              %% Before Equalization
              act_data3 = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
              act_data3 = act_data3(pre_bits+1:end-1);
              pkt = rx_data_dec((kk-1)*pkt_size+1:(kk-1)*pkt_size + pkt_size);
             [snr_basis, dec_bits, sig_power , noise_power, h2, n3] = ...
                 fm0_decode_new_R12(pkt(pre_bits*fm0samp+1-fm0samp:end), ...
                 fm0samp,act_data3);
             tot_bits_b = [tot_bits_b dec_bits];
             bits_b = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
             act_bits_b = [act_bits_b bits_b(pre_bits+1:end-1)];
             ber_b(cnt) = sum(abs(bits_b(pre_bits+1:end-1) - dec_bits))/length(dec_bits);
             snr_vec = [snr_vec snr_basis];
             disp(['Packet BER (Before): ' num2str(ber_b(cnt))]);
             %disp(['Packet SNR in dB   : ' num2str(snr_basis)]);

             %% After Equalization
             pkt = tot_bits_vec((kk-1)*pkt_size+1:(kk-1)*pkt_size + pkt_size).';
             [snr_basis, dec_bits, sig_power , noise_power, h2, n3] = ...
                 fm0_decode_new_R12(pkt(pre_bits*fm0samp+1-fm0samp:end), ...
                 fm0samp,act_data3);
             tot_bits_a = [tot_bits_a dec_bits];
             bits_a = bit_vec_tx((kk-1)*pkt_bits+1:(kk-1)*pkt_bits+pkt_bits);
             act_bits_a = [act_bits_a bits_a(pre_bits+1:end-1)];
             ber_a(cnt) = sum(abs(bits_a(pre_bits+1:end-1) - dec_bits))/length(dec_bits);
             disp(['Packet BER (After): ' num2str(ber_a(cnt))]);
             cnt =cnt + 1;
         end

        ber_fin_b = sum(abs(tot_bits_b - act_bits_b))/length(tot_bits_b);
        disp(['<strong>Final BER (Before):</strong> ' num2str(ber_fin_b)]);
        ber_fin_a = sum(abs(tot_bits_a - act_bits_a))/length(tot_bits_a);
        disp(['<strong>Final BER (After):</strong> ' num2str(ber_fin_a)]);
        disp(['<strong>Final SNR        :</strong> ' num2str(mean(snr_vec))]);
        snr_final = mean(snr_vec);
        ber_vec_a(row,cnttt) = ber_fin_a;
        ber_vec_b(row,cnttt) = ber_fin_b;
        cnttt = cnttt + 1;
        
    end
    row = row + 1;
end

end


folder = '~/Documents/sk/oceans/vanatta/rx_outputs/River PAB Phase Tests 06-06-2022/';
file = 'rx_river_single_phase_pab_007B_?_*deg_18,5kfc_1m_depth_hphydro_.dat';
root = strcat(folder,file);

node_list = ["006A","007A","004A"];

fs = 2e5;
Nsamps = 10*fs;

bp = designfilt('bandpassfir', ...       % Response type
       'StopbandFrequency1',17e3, ...    % Frequency constraints
       'PassbandFrequency1',18e3, ...
       'PassbandFrequency2',19e3, ...
       'StopbandFrequency2',20e3, ...
       'DesignMethod','kaiserwin', ...         % Design method
       'StopbandAttenuation1',40, ...         % Design method options
       'StopbandAttenuation2',40, ...
       'PassbandRipple',1, ...
       'SampleRate',fs);               % Sample rate

d = 0.062;
lambda = 1500/18.5e3;

t_window = 10;
window = chebwin(t_window*fs);
Nfft = length(window);
i = 1;

for ang=-45:45:45
    exp_phase_diff = [0];
    act_phase_diff = [0];

    filename = strrep(root,"*",num2str(ang));

    for n=1:length(node_list)/2
        data = read_complex_binary(strcat(filename,"_",num2str(n-1)));
        node1 = real(data(1:Nsamps));
        node2 = imag(data(1:Nsamps));
    
        node1 = filter(bp,node1);
        node2 = filter(bp,node2);
        
        fft1 = fft(node1.*window,Nfft);
        fft2 = fft(node2.*window,Nfft);
        
        fft1 = fft1(1:Nfft/2);
        fft2 = fft2(1:Nfft/2);
        
        max1 = max(fft1);
        max2 = max(fft2);
        
        meas_ang = angle(max2/max1);
%         if n >=2 && ang ~= 0
%             meas_ang = meas_ang + 2*pi*(n-1)*sign(ang);
%         end

        exp_phase_diff = [exp_phase_diff angle(exp(-1j*2*pi*n*d*sin(ang/180*pi)/lambda))];
        act_phase_diff = [act_phase_diff meas_ang];
    end

    subplot(3,1,i);
    plot(exp_phase_diff/pi*180);
    hold on;
    plot(act_phase_diff/pi*180);
    xlabel("n");
    ylabel("Phase (rad)");
    title(strcat("Phase vs. element ",num2str(ang)," deg"));
    ylim([-180 180]);

    i = i+1;
end

disp("Expected phase difference (deg): ");
disp(num2str(exp_phase_diff/pi*180));

disp("Actual phase difference (deg): ");
disp(num2str(act_phase_diff/pi*180));
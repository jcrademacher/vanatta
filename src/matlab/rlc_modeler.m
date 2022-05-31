node = readmatrix("../../impedance/PAB2_RX_1k-60k_IMP_RIVER_ROTATOR_8020.CSV");
node = node(1:201,:);

close all;

fmin = 10e3;
fmax = 30e3;

f = node(:,1);

ifmin = find(f > fmin,1);
ifmax = find(f < fmax,1,'last');

f = f(ifmin:ifmax);
r = node(ifmin:ifmax,2);
x = node(ifmin:ifmax,3);
ydata = [r,x];
zdata = r+1j*x;
% % % 
%rlc(1) = 200;
% rlc(2) = 100e-5;
% rlc(3) = 11e-9;
% rlc(4) = 3e-9;

rlc0 = rlc;

func_min = @(rlc) rlc_model(rlc,f)-ydata;

options = optimoptions('lsqnonlin','StepTolerance',1e-10,...
    'MaxFunctionEvaluations',1e4,'MaxIterations',1e4,'Algorithm','trust-region-reflective');
[rlc,resnorm,residuals,exitflag,output] = lsqnonlin(func_min,rlc0,[0 0 0 0],[],options);

% zopt = F(rlc,f);
% 
imp_opt = rlc_model(rlc,f);
zopt = imp_opt(:,1)+1j*imp_opt(:,2);

figure(1);
subplot(1,2,1);
plot(f/1e3,real(zdata));
hold on;
plot(f/1e3,real(zopt));
title("Measured Real Part vs. LSQ RLC Real Part Fit");
legend("Measured","LSQ Fit");
xlabel("Frequency (kHz)");
ylabel("Real Part (Ohm)");

subplot(1,2,2);
plot(f/1e3,imag(zdata));
hold on;
plot(f/1e3,imag(zopt));
title("Measured Imaginary Part vs. LSQ RLC Imaginary Part Fit");
legend("Measured","LSQ Fit");
xlabel("Frequency (kHz)");
ylabel("Imaginary Part (Ohm)");

disp(strcat("R = ",num2str(rlc(1)), " ohm"));
disp(strcat("L = ",num2str(rlc(2)/1e-6)," uH"));
disp(strcat("C = ",num2str(rlc(3)/1e-9)," nF"));
disp(strcat("C0 = ",num2str(rlc(4)/1e-9)," nF"));

function yout = rlc_model(rlc,fdata)
    
    wdata = 2*pi*fdata;
    yout = zeros(length(fdata),2); % allocate yout

    R1 = rlc(1);
    L1 = rlc(2);
    C1 = rlc(3);
    C0 = rlc(4);
    
%     yout(:,1) = real(R1+1j*wdata*L1+1./(1j*wdata*C1));
%     yout(:,2) = imag(R1+1j*wdata*L1+1./(1j*wdata*C1));
    yout(:,1) = real((R1+1j*wdata*L1+1./(1j*wdata*C1))./(1+1j*wdata*C0.*(R1+1j*wdata*L1+1./(1j*wdata*C1))));
    yout(:,2) = imag((R1+1j*wdata*L1+1./(1j*wdata*C1))./(1+1j*wdata*C0.*(R1+1j*wdata*L1+1./(1j*wdata*C1))));
end

function [c0,rlc1,rlc2] = rlc_modeler(measured_data,fminmax,init_c0,init_rlc1,init_rlc2,ls,display)

    fmin = fminmax(1);
    fmax = fminmax(2);
    
    f = measured_data(:,1);
    
    ifmin = find(f > fmin,1);
    ifmax = find(f < fmax,1,'last');
    
    f = f(ifmin:ifmax);
    r = measured_data(ifmin:ifmax,2);
    x = measured_data(ifmin:ifmax,3);
    ydata = [r';x'];
    zdata = r+1j*x;
    % % % 
%     n = 2;
    if isempty(init_c0) || isempty(init_rlc1) || isempty(init_rlc2)
        rlc0(1) = 300; % ohm
        rlc0(2) = 6; % mH
        rlc0(3) = 15; % nF
        rlc0(4) = 40; % nF

        rlc0(5) = 1;
        rlc0(6) = 1; 
        rlc0(7) = 1; 

        rlc0 = rlc0';
    else
        rlc0 = [init_rlc1';init_c0;init_rlc2'];
    end
    
    func_min = @(rlc) rlc_model(rlc,f)-ydata;
    
    options = optimoptions('lsqnonlin','FunctionTolerance',1e-12,'StepTolerance',1e-12,...
                        'MaxFunctionEvaluations',10e3,'MaxIterations',10e3,'Display','iter',...
                        'Algorithm','trust-region-reflective');
    [rlc,~,~,~,~] = lsqnonlin(func_min,rlc0,[1 1 1 1 0.01 0.01 0.01]',[1e3 1e3 1e3 1e3 10 10 10],options);
    
    %rlc = rlc0;
    imp_opt = rlc_model(rlc,f);
    zopt = imp_opt(1,:)+1j*imp_opt(2,:);
    
    if display
        figure;
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
    
%         disp(strcat("R = ",num2str(rlc(1)), " ohm"));
%         disp(strcat("L = ",num2str(rlc(2))," mH"));
%         disp(strcat("C = ",num2str(rlc(3))," nF"));
%         disp(strcat("C0 = ",num2str(rlc(4))," nF"));
        %disp(strcat("Ls = ",num2str(rlc(5))," mH"));
    end
    
    rlc1(1) = rlc(1);
    rlc1(2) = rlc(2)*1e-3;
    rlc1(3) = rlc(3)*1e-9;
    c0 = rlc(4)*1e-9;

    rlc2(1) = rlc(5);
    rlc2(2) = rlc(6);
    rlc2(3) = rlc(7);

    function yout = rlc_model(rlc,fdata)

        wdata = 2*pi*fdata;
        yout = zeros(2,length(fdata)); % allocate yout
    
        R1 = rlc(1);
        L1 = rlc(2)*1e-3;
        C1 = rlc(3)*1e-9;
        C0 = rlc(4)*1e-9;

        R2 = R1*rlc(5);
        L2 = L1*rlc(6);
        C2 = C1*rlc(7);

        imp_c0 = 1./(1j*wdata*C0);
        imp_res1 = R1+1j*wdata*L1+1./(1j*wdata*C1);
        imp_res2 = R2+1j*wdata*L2+1./(1j*wdata*C2);

        imp = 1j*wdata*ls+par(imp_c0,par(imp_res1,imp_res2));

        yout(1,:) = real(imp);
        yout(2,:) = imag(imp);

        function ztot = par(z1,z2)
            ztot = z1.*z2./(z1+z2);
        end
    end
end

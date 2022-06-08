global params phimin phimax;
%close all;
params.d = 7e-2;
params.c = 1500;
params.f = 18.5e3;
params.N = 2;
params.A = 10^(-40/20);
params.r = 1;
params.phi = 0;
params.do_direction_val = 0;

dmin = 1e-2;
dmax = 20e-2;

phimin = -pi/2;
phimax = pi/2;

fmin = 5e3;
fmax = 60e3;

atot_db = generate_pattern(params);

global Ntheta theta min_r_db pol pol_steer;

min_r_db = -50;
max_r_db = -30;
Ntheta = 1000;
theta = linspace(phimin,phimax,Ntheta); 

f = figure(1);
ax = axes('Parent',f,'position',[0.13 0.39 0.77 0.54]);
%subplot(1,3,3);
pol = polarplot(theta,atot_db);
hold on; 
pol_steer = polarplot(params.phi*ones(1,10),linspace(min_r_db,max_r_db,10),'--');
rlim([min_r_db max_r_db]);

global dslider_label phislider_label fslider_label;

%sliders and labels
dslider = uicontrol('Parent',f,'Style','slider','Position',[100,800,230,23],...
              'value',params.d, 'min',dmin, 'max',dmax,'Tag','d');
dslider_label = uicontrol('Parent',f,'Style','text','Position',[100,770,230,23],...
                'String',strcat('Element separation=',strcat(num2str(params.d/1e-2),'cm')));

phislider = uicontrol('Parent',f,'Style','slider','Position',[100,740,230,23],...
              'value',params.phi, 'min',phimin, 'max',phimax,'Tag','phi');
phislider_label = uicontrol('Parent',f,'Style','text','Position',[100,710,230,23],...
                'String',strcat('Angle=',strcat(num2str(params.phi*180/pi),'deg')));

fslider = uicontrol('Parent',f,'Style','slider','Position',[100,680,230,23],...
              'value',params.f, 'min',fmin, 'max',fmax,'Tag','f');
fslider_label = uicontrol('Parent',f,'Style','text','Position',[100,650,230,23],...
                'String',strcat('Frequency=',strcat(num2str(params.f/1e3),'kHz')));

Ndrop = uicontrol(f,'Style','popupmenu','Tag','N');
Ndrop.Position = [100,620,230,23];
Ndrop.String = {'2','4','6','8','10','12','14','16','18','20'};

Ndrop_label = uicontrol('Parent',f,'Style','text','Position',[100,590,230,23],...
                'String','N elements');

dslider.Callback = @update_plot; 
phislider.Callback = @update_plot; 
fslider.Callback = @update_plot; 
Ndrop.Callback = @update_plot;


% R1 = 400;
% L1 = 100e-6;
% C1 = 30e-9;
% 
% R2 = 65;
% L2 = 150e-6;
% C2 = 25e-9;
% 
% f = 20e3;
% s = 1j*2*pi*f;
% 
% y = 1/(1/(s*C1)+1/(s*C2)+s*(L1+L2));
% Y = [y -y; -y y];
% 
% F = [1/(2*sqrt(R1)) 0; 0 1/(2*sqrt(R2))];
% ZR = [R1 0; 0 R2];
% 
% i = eye(2);
% S = F*(i-ZR*Y)*(i+ZR*Y)^(-1)*F^(-1);
% 
% gamma1 = S(2,2);
% 
% Zr = R2;
% Zl = R1+1/(s*C1)+1/(s*C2)+s*(L1+L2);
% 
% gamma2 = (Zl-conj(Zr))/(Zl+Zr);

function atot = generate_pattern(params)
    global Ntheta theta min_r_db phimin phimax;

    at = zeros(params.N,Ntheta);    % transmitted wave column vector
    an = zeros(params.N,Ntheta);    % transmitted wave at observed point theta (far field)
    sij = get_sparams([55.9 1350e-6 45e-9],[47.5 1310e-6 44.4e-9],params.f);

    lambda = params.c / params.f;
    k = 2*pi/lambda;

    phi = params.phi;

    if params.do_direction_val == 1
         phi = theta;
    end
    
    for n=1:params.N
        at(n,:) = params.A*(sij(n,params.N+1-n)*exp(-1j*k*sin(phi)*(params.N-n)*params.d)+ ...
                    sij(n,n)*exp(-1j*k*sin(phi)*(n-1)*params.d));
        an(n,:) = at(n,:).*exp(1j*k*(params.r-(n-1)*params.d*sin(theta)));
    end
    
    atot = abs(sum(an,1));
    %atot_norm = atot/max(atot);
    atot = 20*log10(atot);

    direction_val = atot(floor((params.phi-phimin)/(phimax-phimin)*(Ntheta-1)+1));

    disp(['Value in direction: ' num2str(direction_val) 'dB']);

    atot(atot < min_r_db) = min_r_db;
end

function a = update_plot(es,ed)
    global dslider_label phislider_label fslider_label pol pol_steer params;

    switch ed.Source.Tag
        case 'd'
            params.d = es.Value;
            dslider_label.String = strcat('Element separation=',strcat(num2str(es.Value/1e-2),'cm'));
        case 'phi'
            params.phi = es.Value;
            phislider_label.String = strcat('Angle=',strcat(num2str(es.Value*180/pi),'deg'));
        case 'f'
            params.f = es.Value;
            fslider_label.String = strcat('Frequency=',strcat(num2str(es.Value/1e3),'kHz'));
        case 'N'
            params.N = 2*es.Value;
        otherwise
    end

    set(pol,'YData',generate_pattern(params));
    set(pol_steer,'XData',params.phi*ones(1,10));
end



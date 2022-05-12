global params;
close all;
params.d = 7.5e-2;
params.c = 1500;
params.f = 20e3;
params.N = 2;
params.A = 1;
params.r = 1;
params.phi = 0;

dmin = 1e-2;
dmax = 20e-2;

phimin = -pi/2;
phimax = pi/2;

fmin = 5e3;
fmax = 60e3;

atot_db = generate_pattern(params);

global Ntheta theta min_r_db pol;

min_r_db = -60;
Ntheta = 1000;
theta = linspace(-pi/2,pi/2,Ntheta); 

f = figure;
ax = axes('Parent',f,'position',[0.13 0.39 0.77 0.54]);
pol = polarplot(theta,atot_db);
rlim([min_r_db inf]);

global dslider_label phislider_label fslider_label;

% sliders and labels
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


function atot = generate_pattern(params)
    global Ntheta theta min_r_db;

    at = zeros(params.N,1);    % transmitted wave column vector
    an = zeros(params.N,Ntheta);    % transmitted wave at observed point theta (far field)
    sij = flip(eye(params.N)); % scattering coefficients anti-diagonal matrix (assumes wave is perfectly received by node)

    sij(1,2) = sij(1,2)*exp(1j*pi/6);

    lambda = params.c / params.f;
    k = 2*pi/lambda;
    
    for n=1:params.N
        at(n) = sij(n,params.N+1-n)*params.A*exp(-1j*k*sin(params.phi)*(params.N-n)*params.d);
        an(n,:) = at(n)*exp(1j*k*(params.r-(params.N-n)*params.d*sin(theta)));
    end
    
    atot = abs(sum(an,1));
    %atot_norm = atot/max(atot);
    atot = 20*log10(atot);

    disp(['Max array val: ' num2str(max(atot)) 'dB']);

    atot(atot < min_r_db) = min_r_db;
end

function a = update_plot(es,ed)
    global dslider_label phislider_label fslider_label pol params;

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
end


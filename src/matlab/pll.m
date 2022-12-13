fs = 2e5;
t_tot = 10e-3;
N = fs*t_tot;

f0 = 3e3;

ph_err = [];

t = zeros(1,N);
vco = zeros(1,N);
ph = zeros(1,N);
ph_est = zeros(1,N);
lp = zeros(1,N);
y = zeros(1,N);
integ = zeros(1,N);

Bn = 0.01*fs;
damp = 1;

k0 = 1;
kd = 1;
kp = 1/(kd*k0)*4*damp/(damp+1/(4*damp))*Bn/fs;
ki = 1/(kd*k0)*4/(damp+1/(4*damp))^2*(Bn/fs)^2;

integ_out = 0;

for n = 1:N-1
    t(n) = t_tot*n/N;
    % input signal
    y(n) = sin(2*pi*f0*t(n)+pi);

    % phase detect
    ph(n) = kd*y(n)*imag(vco(n));

    % loop filter
    integ_out = ki*ph(n)+integ_out;
    lp(n) = kp*ph(n) + integ_out;

    % vco
    ph_est(n+1) = ph_est(n) + k0*lp(n);
    vco(n+1) = exp(-1j*(2*pi*f0*t_tot*(n+1)/N+ph_est(n)));
end


plot(t(1:end-1),y(1:end-1));
hold on;
plot(t(1:end-1),real(vco(1:end-1)));
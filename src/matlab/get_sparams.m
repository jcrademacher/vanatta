function sp = get_sparams(varargin)
% function assumes the basic butterworth-van dyke model of a piezo (static
% C0, parallel RLC network representing mechanical resonance)
if mod(nargin,2) == 0
    error("Odd # of arguments required. Make sure # of nodes is even.")
end

N = nargin-1;
sp = zeros(N,N);
f = varargin{nargin};
w = 2*pi*f;

for n=1:N/2
    rlc_c0_1 = varargin{n};
    rlc_c0_2 = varargin{N-n+1};

    r1 = rlc_c0_1(1);
    l1 = rlc_c0_1(2);
    c1 = rlc_c0_1(3);
    c01 = rlc_c0_1(4);
    
    r2 = rlc_c0_2(1);
    l2 = rlc_c0_2(2);
    c2 = rlc_c0_2(3);
    c02 = rlc_c0_2(4);
    
    Zr = [r1 0; 0 r2];
    F = [1/(2*sqrt(r1)) 0; 0 1/(2*sqrt(r2))];
    
    z12 = 1./(1j*w*(c01+c02));
    z21 = z12;
    z11 = 1./(1j*w*c1)+1j*w*l1+z12;
    z22 = 1./(1j*w*c2)+1j*w*l2+z12;
    
    Z = [z11 z12;z21 z22];
    Y = Z^(-1);
    
    i = eye(2);
    
    sp(n:(N-2*n+1):(N-n+1),n:(N-2*n+1):(N-n+1)) = F*(i-conj(Zr)*Y)*(i+Zr*Y)^(-1)*F^(-1);
end



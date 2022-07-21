function sp = get_sparams(varargin)

if mod(nargin,2) == 0
    error("Odd # of arguments required. Make sure # of nodes is even.")
end

N = nargin-1;
sp = zeros(N,N);
f = varargin{nargin};

for n=1:N/2
    imp1 = varargin{n};
    imp2 = varargin{N-n+1};

    r1 = real(imp1);
    j1 = imag(imp1);
    
    r2 = real(imp2);
    j2 = imag(imp2);
    
    Zr = [r1 0; 0 r2];
    F = [1/(2*sqrt(r1)) 0; 0 1/(2*sqrt(r2))];
    
    y11 = (j1+j2)^(-1);
    
    Y = [y11 -y11;-y11 y11];
    i = eye(2);
    
    sp(n:(N-2*n+1):(N-n+1),n:(N-2*n+1):(N-n+1)) = F*(i-conj(Zr)*Y)*(i+Zr*Y)^(-1)*F^(-1);
end



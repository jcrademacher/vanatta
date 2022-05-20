function sp = get_sparams(varargin)

if mod(nargin,2) == 0
    error("Odd # of arguments required. Make sure # of nodes is even.")
end

N = nargin-1;
sp = zeros(N,N);
f = varargin{nargin};

for n=1:N/2
    node1 = varargin{n};
    node2 = varargin{N-n+1};

    r1 = node1(1);
    l1 = node1(2);
    c1 = node1(3);
    
    r2 = node2(1);
    l2 = node2(2);
    c2 = node2(3);
    
    Zr = [r1 0; 0 r2];
    F = [1/(2*sqrt(r1)) 0; 0 1/(2*sqrt(r2))];
    
    y11 = (1j*2*pi*f*(l1+l2)+(c1+c2)/(1j*2*pi*f*c1*c2))^(-1);
    
    Y = [y11 -y11;-y11 y11];
    i = eye(2);
    
    sp(n:(N-2*n+1):(N-n+1),n:(N-2*n+1):(N-n+1)) = F*(i-conj(Zr)*Y)*(i+Zr*Y)^(-1)*F^(-1);
end



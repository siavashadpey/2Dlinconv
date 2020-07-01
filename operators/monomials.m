function [V, Vx] = monomials(x,p)

% Not to be used for high-orders, since the Vandermonde matrices will be
% ill-conditioned

dim = size(x,2);

switch dim
    case 1
        [V, Vx] = monomials_1D(x,p);
    case 2
        [V, Vx] = monomials_2D(x,p);
end

end

function [V, Vx] = monomials_1D(x,p)

N = size(x,1);

V = zeros(N,p+1);
Vx = V;

for ii=0:p
    kk = ii + 1;
    V(:,kk) = x.^ii;
    Vx(:,kk) = ii*x.^(ii-1);
end

end

function [V, Vx] = monomials_2D(xx,p)

Nc = (p+2)*(p+1)/2;
dim = 2;
x = xx(:,1); y = xx(:,2);
N = size(x,1);

zz = zeros(N,1);
V = zeros(N,Nc);
Vx = zeros(N,Nc,dim);

for jj=0:p
    for ii=0:jj
        kk = jj*(jj+1)/2 + ii + 1;
        V(:,kk)    = x.^ii.*y.^(jj-ii);
        if ii==0
            Vx(:,kk,1) = zz;
        else
            Vx(:,kk,1) = ii*x.^(ii-1).*y.^(jj-ii);
        end
        if (jj-ii)==0
            Vx(:,kk,2) = zz;
        else
            Vx(:,kk,2) = (jj-ii)*x.^ii.*y.^(jj-ii-1);
        end
    end
end

end
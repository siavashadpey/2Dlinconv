function [SBP] = gen_SBP(params)

% Description: Generates reference SBP operator of degree params.p and type
%              params.type

% TODO:
% Use Legendre poly instead of monomials
% Generate Fperm and Vperm more elegantly
% More sig figs for Omega cub rules

% Parameters
p = params.p;
type = params.type;
dim = 2;
nface = 3;
ncar = (p+1)*(p+2)/2;

% Cubature weights and nodes on ref element
[xy, w, Fperm, Vperm] = cub_rule(p,type);

% Quadrature rule and nodes on ref face
[xf, wf] = quad_line(p+1); % p+1 facet nodes
yf = -1*ones(length(xf),1);
xf = [xf yf];

% Vandermonde and derivative matrices 
[V, Vx] = monomials(xy,p); % Change to Legendre poly
[Vf] = monomials(xf,p);    % Change to Legendre poly

% Construct R (for first facet)
nv_f1 = Fperm(:,1); % Indices of volume nodes involved
r = Vf/V(nv_f1,:);

% Normal vector scaled by the Jacobian
N = [0 1 -1; -1 1 0];

% Construct E
nnode = length(w);
E = zeros(nnode,nnode,dim);
for f=1:nface
    nv_f = Fperm(:,f);
    nor = N(:,f);
    E(nv_f,nv_f,1) = E(nv_f,nv_f,1) + r'*diag(nor(1)*wf)*r;
    E(nv_f,nv_f,2) = E(nv_f,nv_f,2) + r'*diag(nor(2)*wf)*r;
end

% EE = zeros(nnode,nnode,dim);
% for f=[2 3 1]
%     nv_f = Fperm(:,f);
%     for i=1:length(nv_f)
%         for j=1:length(nv_f)
%             EE(nv_f(i),nv_f(j),1) = EE(nv_f(i),nv_f(j),1) + r(:,i)'*diag(N(1,f)*wf)*r(:,j);
%             EE(nv_f(i),nv_f(j),2) = EE(nv_f(i),nv_f(j),2) + r(:,i)'*diag(N(2,f)*wf)*r(:,j);
%         end
%     end
% end
% Add case where n > Nc
Nw = nnode - ncar;
if Nw > 0 
    W = null(V'); %??
    W = W(:,1:Nw);
    
    Vtil = [V W];
    
    zz = zeros(Nw,Nw);
    Wx = zeros(nnode, Nw, dim);
    XX = Wx;
    
    XX(:,:,1) = [-Vx(:,:,1)'*diag(w)*W + 0.5*V'*E(:,:,1)*W; zz];
    XX(:,:,2) = [-Vx(:,:,2)'*diag(w)*W + 0.5*V'*E(:,:,2)*W; zz];
    
    Wx(:,:,1) = (Vtil'*diag(w))\(0.5*Vtil'*E(:,:,1)*W + XX(:,:,1));
    Wx(:,:,2) = (Vtil'*diag(w))\(0.5*Vtil'*E(:,:,2)*W + XX(:,:,2));
    
    Vxtil = [Vx Wx];
else
    Vtil = V;
    Vxtil = Vx;
end

% Construct D
D = zeros(nnode,nnode,dim);
D(:,:,1) = Vxtil(:,:,1)/Vtil;
D(:,:,2) = Vxtil(:,:,2)/Vtil;

% Construct Q
Q = zeros(nnode,nnode,dim);
Q(:,:,1) = diag(w)*D(:,:,1);
Q(:,:,2) = diag(w)*D(:,:,2); 

% Outputs
SBP.x = xy;
SBP.xf = xf(:,1);
SBP.h = w;
SBP.Q = Q;
SBP.r = r;
SBP.b = wf;
SBP.N = N(:,1);
SBP.Fperm = Fperm;
SBP.Vperm = Vperm;

end
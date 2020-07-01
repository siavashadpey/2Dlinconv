function [Jac,metr,N] = jactrans(coord,Fperm,ref_SBP)

% Description: computes the Jacobian transformation info and the normal
% components of the facets. These info are used to transform from the ref
% SBP to the physical SBP. 

dim = 2;
nnode = size(coord,1);
nnodef = size(ref_SBP.b,1);
nface = size(Fperm,2);
Jacv = zeros(nnode,dim,dim);
Jacf = zeros(nnodef,dim,nface); 
detJacf = zeros(nnodef,nface);
metr = zeros(nnode,dim,dim);
un = ref_SBP.Fperm;
Vperm = ref_SBP.Vperm;

% H, Q, and D of reference element
Href = diag(ref_SBP.h);
Qref = ref_SBP.Q;
D = zeros(nnode,nnode,dim);
for dd=1:dim
    D(:,:,dd) = Href\Qref(:,:,dd);
end

% Jacobian computed at the volume nodes
for dd=1:dim
    Jacv(:,:,dd) = D(:,:,dd)*coord; 
end
% Determinant of the Jacobian computed at the volume nodes
detJacv = Jacv(:,1,1).*Jacv(:,2,2) - Jacv(:,2,1).*Jacv(:,1,2); 

% Jacobian and determinant of the Jacobian computed at the facet nodes
R = zeros(nnodef,nnode);
R(:,Fperm(:,1)) = ref_SBP.r;
for f=1:nface
    Jacf(:,:,f) = R*D(:,:,1)*coord(Vperm(:,f),:); % Permute the coordinates instead of permuting R and D
    detJacf(:,f) = sqrt(Jacf(:,1,f).^2 + Jacf(:,2,f).^2);
end

% Compute the normal components at each facet node of each facet
N = zeros(nnodef,dim,nface);
for f=1:nface
    N(:,:,f) = [Jacf(:,2,f) -Jacf(:,1,f)]./[detJacf(:,f) detJacf(:,f)]; 
end

% metr is used for transforming the skew-symmetric part of Q, i.e. S.
% (see Crean et al. 2018)
A = Jacv(:,1,1); B = Jacv(:,1,2); C = Jacv(:,2,1); D = Jacv(:,2,2);
metr(:,:,1) = [D -B]; % [dxi/dx  dxi/dy]
metr(:,:,2) = [-C A]; % [deta/dx deta/dy] % includes the Jacobian (it actually cancels out when inverting Jacv)

Jac.jacv = Jacv;
Jac.jacf = Jacf;
Jac.detjacv = detJacv;
Jac.detjacf = detJacf;
end
function ref = gen_ref_tri(x)

% Description: computes linear shape functions and derivative of the shape
% functions evaluated at coordinates x

vtx = [-1  1 -1; 
       -1 -1  1 ]';

[Cinv,~] = monomial_tri(vtx);
[psi,psix] = monomial_tri(x); % In FEM, this is usually the cub nodes

ref.shp = psi/Cinv;  
ref.shpx(:,:,1) = psix(:,:,1)/Cinv;
ref.shpx(:,:,2) = psix(:,:,2)/Cinv; 
ref.x = x;

% vtxf = [-1; 1];
% [Cinv,~] = monomial_line(vtxf);
% [psif,psixf] = monomial_line(xf);
% ref.shpf = psif/Cinv;  
% ref.shpxf = psixf(:,:,1)/Cinv;
% ref.xf = xf;


end
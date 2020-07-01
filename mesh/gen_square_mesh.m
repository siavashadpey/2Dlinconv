function mesh = gen_square_mesh(h,type)
% Discription: generates a square mesh on [1,1]^2 domain using h as nominal
% mesh spacing
% INPUT
%   h: approximate diameter of elements
% OUTPUT
%   mesh: mesh structure

if (nargin < 1)
    error('Specify nominal mesh spacing h.')
end

switch type
    case 'structured'
        ne1d = ceil(1.0/h);
        x1 = linspace(0,1,ne1d+1);
        [xx1,xx2] = ndgrid(x1,x1);
        coord = [xx1(:), xx2(:)];
        con = delaunay(coord);
    case 'unstructured'
        figh = figure;
        fd=@(p) drectangle(p,0,1,0,1);
        [coord,con]=distmesh2d(fd,@huniform,h,[0,0;1,1],[0,0;0,1;1,0;1,1]);
        close(figh);
end


% create boundary edge groups
edge = [con(:,[1,2])   % First  local facet
        con(:,[2,3])   % Second local facet
        con(:,[3,1])]; % Third  local facet
    
tol = 1e-6;

% mid edge coordinate
xe = reshape(mean(reshape(coord(edge(:),:),[size(edge),2]),2),size(edge));

% left
ii = abs(xe(:,1) - 0.0) < tol;
bgrp{1} = sortrows(edge(ii,:),1);

% right
ii = abs(xe(:,1) - 1.0) < tol;
bgrp{2} = sortrows(edge(ii,:),1);

% bottom 
ii = abs(xe(:,2) - 0.0) < tol;
bgrp{3} = sortrows(edge(ii,:),1);

% top
ii = abs(xe(:,2) - 1.0) < tol;
bgrp{4} = sortrows(edge(ii,:),1);

% Element-to-element-to-face info
conf = connect_faces(con);

% Output
mesh.coord = coord;
mesh.con = con;
mesh.conf = conf;
mesh.bgrp = bgrp;
mesh.bgrp = make_bgrp(mesh);
end
function [mesh] = update_mesh_info_old(mesh,ref_tri,nnode)

[nelem, nshp] = size(mesh.con);
nmiss = nnode - nshp;
[ntotal,dim] = size(mesh.coord);
nn = nshp+1:nnode;
ndof = ntotal+nmiss*nelem;
shp = ref_tri.shp(nn,:);

mesh.vtx = mesh.con;
mesh.vtxcoord = mesh.coord;

coord = zeros(ndof,dim); 
coord(1:ntotal,:) = mesh.coord;
con = zeros(nelem,nnode); 
con(:,1:nshp) = mesh.con;

for k=1:nelem
    % EtoN Con
    ignodes = ntotal+1:ntotal+nmiss;
    con(k,nn) = ignodes;
    ntotal = ignodes(end);
    
    % Coord
    ignodes1 = mesh.con(k,:);
    coord1 = mesh.coord(ignodes1,:);
    %coord1 = [-1 -1; 1 -1; -1 1]; % for debugging
    coordn = shp*coord1;
    coord(ignodes,:) = coordn;

end

% Store coordinates and EtoN con table
mesh.con = con;
mesh.coord = coord;

end
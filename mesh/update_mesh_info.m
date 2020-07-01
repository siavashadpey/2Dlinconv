function [mesh] = update_mesh_info(mesh,ref_tri,nnode)

% Description: Updates nodes global indicies and coordinates by simply
% using the elements' vertices coordinates

dim = 2;
[nelem, ~] = size(mesh.con);
ntotal = 0;
nn = 1:nnode;
shp = ref_tri.shp(nn,:);

% Vertices (i.e. old mesh info)
mesh.vtx = mesh.con;
mesh.vtxcoord = mesh.coord;

% New mesh info
con = zeros(nelem,nnode); 
coord = Inf(size(mesh.coord,1),dim);
for k=1:nelem
    % EtoN Con
    ignodes = ntotal+1:ntotal+nnode; % Global indices of current elem
    ig_new = ignodes; 
    n_new = length(ig_new);
    
    % Coordinates of vertices
    igvtx = mesh.vtx(k,:);
    coord1 = mesh.vtxcoord(igvtx,:);

    % Coordinates of element's nodes
    coordn = shp*coord1;
    
    if k>1
        % Identify non-unique nodes
        % Replace their indices with the correpsonding existing node indices
        [tf, ind] = ismember(coordn,coord,'rows'); ind = ind';
        ind_ow = find(tf~=0); % Non-unique nodes
        ignodes(ind_ow) = ind(ind~=0);
        
        % Identify new coordinates
        % Assign them unique node indices (requires updating ig_new)
        [coordn] = setdiff(coordn,coord,'rows','stable');
        n_new = size(coordn,1);
        ig_new = ntotal+1:ntotal+n_new;
        ind_ow = find(tf==0); % Unique nodes
        ignodes(ind_ow) = ig_new;
    end
    
    % Store coordinates of new nodes
    coord(ig_new,:) = coordn;
    
    % Curve coordinates
    %coord(ig_new,1) = coord(ig_new,1) - 1/30*sin(2*pi*coord(ig_new,1)).*cos(2*pi*coord(ig_new,2));
    %coord(ig_new,2) = coord(ig_new,2) - 1/30*cos(2*pi*coord(ig_new,1)).*sin(2*pi*coord(ig_new,2));
    
    % Store element's global node indices
    con(k,nn) = ignodes;
    
    % Update total number of unique nodes
    ntotal = ntotal + n_new;
end

% Store coordinates and EtoN con table
mesh.con = con;
mesh.coord = coord;

end
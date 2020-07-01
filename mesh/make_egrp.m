function mesh = make_egrp(mesh)
% MAKE_EGRP adds edge-group structure to a mesh for DG discretization
% REMARKS: two arrays are added to the mesh
%   t2e: nelem by 3 array of triangle-to-edge connectivity.  t2e(i,j) is
%        the global edge number of the j-th local edge of the i-th element
%   egrp: nedge by 6 array of edge-to-triangle connecivity.
%         egrp(i,[1:2]) is the two vertices that comprise the edge
%         egrp(i,3) is the global element index of the left element
%         egrp(i,4) is the local edge index of the left element
%         egrp(i,5) is the global element index of the right element.  If
%           this a boundary edge, then the number is the negative of the
%           boundary group number.
%         egrp(i,6) is the local edge index of the right element

% Copyright 2018 Masayuki Yano, University of Toronto

tri = mesh.con;
ntri = size(tri,1);
tedge = [tri(:,[1,2])
         tri(:,[2,3])
         tri(:,[3,1])];
[edge,~,ie] = unique(sort(tedge,2),'rows');
t2e = reshape(ie,[ntri,3]); % ie(elem,ledge) is the global edge number
nedge = size(edge,1);

% create edge group structure
egrp = zeros(nedge,6);
for elem = 1:ntri
    for ledge = 1:3
        iedge = t2e(elem,ledge);
        nodes = tedge(elem+ntri*(ledge-1),:);
        if (egrp(iedge,1) == 0)
            egrp(iedge,1:2) = nodes;
            egrp(iedge,3:4) = [elem,ledge];
        else
            egrp(iedge,5:6) = [elem,ledge];
            if any(nodes~=egrp(iedge,[2,1]))
                error('something is wrong');
            end
        end
    end
end

% fill boundary information
for ibgrp = 1:length(mesh.bgrp)
    [~,ie] = intersect(egrp(:,1:2),mesh.bgrp{ibgrp}(:,1:2),'rows');
    egrp(ie,5) = -ibgrp;
end
if any(egrp(:,5)==0) || any(egrp(:,5)<0 & egrp(:,6)~=0)
    error('something is wrong');
end

mesh.t2e = t2e;
mesh.egrp = egrp;

end
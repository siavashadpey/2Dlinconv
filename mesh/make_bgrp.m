function bgrp = make_bgrp(mesh)
% Description: Adds element number and local facet number to each boundary
%              group 

con = mesh.con;
bgrp = mesh.bgrp;
ncon = size(con,1);
edge = [con(:,[1,2])
        con(:,[2,3])
        con(:,[3,1])];
for ibgrp = 1:length(mesh.bgrp)
    bvert = sort(mesh.bgrp{ibgrp}(:,1:2),2);
    edge = sort(edge,2);
    [~,ib] = intersect(edge,bvert,'rows');
    [belem,bledge] = ind2sub([ncon,3],ib);
    bgrp{ibgrp} = [belem, bledge];
end
end


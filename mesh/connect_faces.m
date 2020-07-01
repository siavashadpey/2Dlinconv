function [conf]= connect_faces(con)
% Description: : Local facet nb i of element k coincides with local facet
%                nb j of element v where 
%                v = conf(k,i,1) and j = conf(k,i,2)
% Adapted from Nodal Discontinuous Galerkin Methods: Algorithms, Analysis, 
%              and Applications by Hesthaven and Warburton (www.nudg.org/)
nface=3;
K = size(con,1);
nnode = max(max(con));

% Edges
edge = [con(:,[1,2])
        con(:,[2,3])
        con(:,[3,1])];
edge = sort(edge,2)-1;

% Initialize element to element and Element to faces connectivity
EToE = (1:K)'*ones(1,nface); EToF = ones(K,1)*(1:nface);

% Uniquely number each set of three faces by their node numbers 
id = edge(:,1)*nnode + edge(:,2)+1;
spNodeToNode=[id, (1:nface*K)', EToE(:), EToF(:)];

% Sort by global face number
sorted=sortrows(spNodeToNode,1);

% Find matches in the sorted face list
[indices,~]=find( sorted(1:(end-1),1)==sorted(2:end,1) );

% Make links reflexive 
matchL = [sorted(indices,:)   ;sorted(indices+1,:)];
matchR = [sorted(indices+1,:) ;sorted(indices,:)];

% Insert matches
EToE(matchL(:,2)) = matchR(:,3); EToF(matchL(:,2)) = matchR(:,4);

conf(:,:,1) = EToE; conf(:,:,2) = EToF;

end
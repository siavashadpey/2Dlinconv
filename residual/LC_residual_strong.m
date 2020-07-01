function Res = LC_residual_strong(mesh,ref_SBP,U,eqn,t)

% Useful parameters
Fperm = ref_SBP.Fperm;           % Needed for SBP R matrix
nelem = size(mesh.con,1);        % Nb of elements
nnode = size(ref_SBP.h,1);       % Nb of volume nodes/element
nnodef = size(ref_SBP.b,1);      % Nb of facet nodes/face
nface = size(ref_SBP.Fperm,2);   % Nb of faces/element
Res_k = zeros(nnode,nelem);      % Residual matrix (nnode x nelem)

UK = reshape(U,[nnode,nelem]); 

% Loop through each element
% Construct physical SBP operators
% Add volume terms and surface terms
for k=1:nelem
    ignodes = mesh.con(k,:)';      % Global indices of volume nodes
    coord = mesh.coord(ignodes,:); % Coord of volume nodes

    uk = UK(:,k);                  % Solution at volume nodes
    spd_a = eqn.a(coord);          % Wave speed
    
    % Construct physical SBP operators
    [jac,~,N] = jactrans(coord,Fperm,ref_SBP);
    [H,~,~,Dx,Dy,~,~] = gen_physSBP(ref_SBP,jac);
    
    % Add volume terms to residual
    fx = spd_a(:,1).*uk; fy = spd_a(:,2).*uk;
    Res_k(:,k) = - (Dx*fx + Dy*fy); % Strong form
    
    % Compute surface terms
    SAT = zeros(nnode,1);
    % Loop through each face
    for f_k=1:nface
        detJacf = jac.detjacf(:,f_k);
        Rk = zeros(nnodef,nnode); 
        FPk = Fperm(:,f_k)'; 
        Rk(:,FPk) = ref_SBP.r;    % Projection operator of facet f_k
        uk_f = Rk*uk;             % Solution at facet nodes
        xxyyq = Rk*coord;         % Coordinates of facet nodes
        spd_a = eqn.a(xxyyq);     % Wave speed at facet nodes
        v = mesh.conf(k,f_k,1);   % Neighbouring element index
        if v==k % Non-periodic boundary interface
            if (eqn.isPBC)
                continue
            end
            uv_f = eqn.BC(xxyyq,t);
        else % Inter-element coupling
            ig_v = mesh.con(v,:)';        % Neighbouring element node indices
            coord_v = mesh.coord(ig_v,:); % Neighbouring element's node coord
            f_v = mesh.conf(k,f_k,2);     % Neighbouring element local facet index
            FPv = Fperm(:,f_v)'; 
            Rv = zeros(nnodef,nnode); Rv(:,FPv) = ref_SBP.r;
            uv = UK(:,v);
            uv_f = Rv*uv; uv_f = flipud(uv_f); % Need to flip vector since CCW is used for both elements
%             xxyyq_v = Rv*coord_v; xxyyq_v = flipud(xxyyq_v);
%             diff = norm(xxyyq_v - xxyyq);
%             if diff > 1E-12
%                 error('Facet coordinates do not coincide on element %d with local face %d. \n',k,f_k)
%             end
        end
        % Numerical flux (Roe scheme)
        a_n = spd_a(:,1).*N(:,1,f_k) + spd_a(:,2).*N(:,2,f_k);
        fs = 0.5*(a_n.*(uk_f+uv_f) + abs(a_n).*(uk_f-uv_f));
        fn = a_n.*uk_f;
        % Integrate flux over surface
        cont = Rk'*diag(detJacf.*ref_SBP.b)*(fs - fn);  
        SAT = SAT + cont;
    end
    % Add SAT to residual
    Res_k(:,k) = Res_k(:,k) - H\SAT;
end

% Periodic BCs (assuming structured square grid)
if (eqn.isPBC)
    for ibgrp=[1,3]
        bgrpK = mesh.bgrp{ibgrp};
        bgrpV = mesh.bgrp{ibgrp+1};
        % Loop through each facet
        for iedge=1:size(bgrpK,1)
            k = bgrpK(iedge,1); f_k = bgrpK(iedge,2); % Facet's element nb and its local facet index
            v = bgrpV(iedge,1); f_v = bgrpV(iedge,2); % Facet's neighbouring element nb and its local facet index

            ignodes = mesh.con(k,:)';      
            coord = mesh.coord(ignodes,:); 

            uk = UK(:,k);

            % Calculate H and the determinant of the Jacobians at the facets
            [jac_k,~,N] = jactrans(coord,Fperm,ref_SBP);
            [Hk] = gen_physSBP(ref_SBP,jac_k);

            detJacf = jac_k.detjacf(:,f_k);
            Rk = zeros(nnodef,nnode); 
            FPk = Fperm(:,f_k); 
            Rk(:,FPk) = ref_SBP.r;
            uk_f = Rk*uk;
            xxyyq = Rk*coord;
            spd_a = eqn.a(xxyyq);
 
            FPv = Fperm(:,f_v);
            Rv = zeros(nnodef,nnode); Rv(:,FPv) = ref_SBP.r;
            uv = UK(:,v);
            uv_f = Rv*uv; uv_f = flipud(uv_f);
            
            % Calculate Hv
            [jac_v,~,~] = jactrans(coord_v,Fperm,ref_SBP); % Fperm is not used for jac_v
            [Hv] = gen_physSBP(ref_SBP,jac_v);

            a_n = spd_a(:,1).*N(:,1,f_k) + spd_a(:,2).*N(:,2,f_k);
            fs = 0.5*(a_n.*(uk_f+uv_f) + abs(a_n).*(uk_f-uv_f));
            fn_k = a_n.*uk_f;
            fn_v = a_n.*uv_f;
            
            % Integrate flux over surface;
            contK = Rk'*diag(detJacf.*ref_SBP.b)*(fs - fn_k);
            contV = -Rv'*diag(detJacf.*ref_SBP.b)*(flipud(fs - fn_v));
            
            Res_k(:,k) = Res_k(:,k) - Hk\contK;
            Res_k(:,v) = Res_k(:,v) - Hv\contV;
        end
    end
end

% Output residual as a vector
Res = Res_k(:);

end
function err_nor = err_H(mesh,ref_SBP,eqn,t,U)

% Description: Outputs the H-norm of the error in the numerical solution.
% The error is normalized with the H-norm of the exact solution.

% Useful parameters
nelem = size(mesh.con,1);
nnode = size(ref_SBP.h,1);
Fperm = ref_SBP.Fperm;
UK = reshape(U,[nnode,nelem]);

err = 0; 
u0_H = 0; 

% Compute the square of the H-norm of the error and the exact solution
for k=1:nelem
    ignodes = mesh.con(k,:)';
    coord = mesh.coord(ignodes,:);
    Uk = UK(:,k);
    Uex = eqn.u(coord,t);
    err_k = Uk - Uex;
    % Calculate physical H
    [jac,~,~] = jactrans(coord,Fperm,ref_SBP);
    [H] = gen_physSBP(ref_SBP,jac);
    
    err = err + err_k'*H*err_k;
    %u0_H = u0_H + Uex'*H*Uex;
    
end

% Normalize error
err_nor = sqrt(err); %/sqrt(u0_H);
    

end
function [H,Qx,Qy,Dx,Dy,Ex,Ey] = gen_physSBP(ref_SBP,jac)

% Description: Generates a physical SBP operator from a reference SBP
% operator by using the Jacobian transformation info. 

href = ref_SBP.h;
Qref = ref_SBP.Q;
jacv = jac.jacv;
detJacv = jac.detjacv;

nnode = size(href,1);
dim = 2;

H = diag(href.*detJacv);

Dref = zeros(nnode,nnode,dim); Eref = Dref;
for d=1:dim
    Dref(:,:,d) = diag(1./href)*Qref(:,:,d);
    Eref(:,:,d) = Qref(:,:,d) + Qref(:,:,d)';
end


% J11 = repmat(jacv(:,1,1),[1 nnode]); 
% J12 = repmat(jacv(:,1,2),[1 nnode]); 
% J21 = repmat(jacv(:,2,1),[1 nnode]); 
% J22 = repmat(jacv(:,2,2),[1 nnode]); 
% JJ  = repmat(detJacv, [1 nnode]);
% Dx = 1./JJ.*( J22.*Dref(:,:,1) - J21.*Dref(:,:,2));
% Dy = 1./JJ.*(-J12.*Dref(:,:,1) + J11.*Dref(:,:,2));
% Ex =  J22.*Eref(:,:,1) - J21.*Eref(:,:,2);
% Ey = -J12.*Eref(:,:,1) + J11.*Eref(:,:,2);

J11 = diag(jacv(:,1,1)); 
J12 = diag(jacv(:,1,2)); 
J21 = diag(jacv(:,2,1)); 
J22 = diag(jacv(:,2,2)); 
JJinv  = diag(1./detJacv);
Dx = JJinv*( J22*Dref(:,:,1) - J21*Dref(:,:,2));
Dy = JJinv*(-J12*Dref(:,:,1) + J11*Dref(:,:,2));
Ex =  J22*Eref(:,:,1) - J21*Eref(:,:,2);
Ey = -J12*Eref(:,:,1) + J11*Eref(:,:,2);

Qx = H*Dx; Qy = H*Dy;

end
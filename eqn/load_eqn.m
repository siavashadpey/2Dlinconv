function [eqn] = load_eqn(eqntyp,params)

switch eqntyp
    case 'const_lin_conv'
        eqn.a = @(x) [params(1), params(2)]; 
        a = [params(1); params(2)];
        eqn.u = @(x,t) sin(2*pi.*(x(:,1)-a(1).*t)).*sin(2*pi.*(x(:,2)-a(2).*t));
        eqn.BC = eqn.u;
        eqn.isPBC = false;
        eqn.isexact = true;
    case 'const_per_lin_conv'
        eqn.a = @(x) [params(1), params(2)]; 
        a = [params(1); params(2)];
        eqn.u = @(x,t) sin(2*pi/2.*(x(:,1)-a(1).*t)).*sin(2*pi/2.*(x(:,2)-a(2).*t));
        eqn.isPBC = true;
        eqn.isexact = true;
    case 'circ_lin_conv'
        eqn.a = @(x) [-x(:,2), x(:,1)];
        eqn.u = @(x,t) exp(-sum(bsxfun(@minus,x,[0.3,0.0]).^2,2)/0.2^2);
        eqn.BC = @(x,t) 0;
        eqn.isPBC = false;
        eqn.isexact = false;
    case 'bump'
        eqn.a = @(x) [params(1), params(2)]; 
        eqn.u = @(x,t) bump_fn(x);
        eqn.isPBC = true;
        eqn.isexact = true;
end

end

function u = bump_fn(x)
% From first multiD SBP paper. Hicken et al.
rr = sqrt((x(:,1)-.5).^2 + (x(:,2)-.5).^2);
u(rr<=0.5) = 1-(4.*rr(rr<=0.5).^2-1).^5;
u(rr>0.5) = 1;
u = u';
end


function [err] = lin_conv_2D(SBP,dx,eqntyp,spd_a,form,tf)
close all

% Input arguments
if nargin<1
    SBP.p = 3; SBP.type = 'omega';   % SBP operator info
end
if nargin<2
    dx = 1/12;                      % Mesh spacing
end
if nargin<3
    eqntyp = 'const_lin_conv';                 % Problem to solve
end
if nargin<4
    spd_a.x = .5; spd_a.y = spd_a.x;        % Wave speed
end
if nargin<5
    form = 'weak';                 % Strong or weak forms
end
if nargin<6
    tf = .5;                          % Final time
end

is_debug = false; % To plot mesh

% SBP operators on reference element (H, Q, R, B, Nxy, xy)
ref_SBP = gen_SBP(SBP); 

% Generate linear shape functions
% Only used for calculating physical location of cubature nodes as a fn of
% x_vtx (linear mapping)
ref_tri = gen_ref_tri(ref_SBP.x);

% Generate mesh 
mesh_typ = 'structured';
mesh = gen_square_mesh(dx,mesh_typ);

% Load PDE eqn
switch eqntyp
    case {'const_lin_conv','const_per_lin_conv','bump'}
        param = [spd_a.x spd_a.y];
        eqn = load_eqn(eqntyp,param);
    case 'circ_lin_conv'
        eqn = load_eqn(eqntyp);
        mesh.coord = bsxfun(@plus,2*mesh.coord,[-1,-1]);
end

% Get element-to-node connectivity and coordinates of all nodes
nnode = size(ref_SBP.h,1);
mesh = update_mesh_info(mesh,ref_tri,nnode); % Doesn't change elem indices nor local facet indices

% Plot mesh
if (is_debug)
    opt.number_nodes = true; opt.number_elems = false; 
    plot_mesh(mesh,opt); 
end

% Choose weak or strong form (algebraically the same)
% Note: strong form is more efficient, since there is no need to loop
% through the elements and invert the H matrix on the LHS of the eq
switch form
    case 'weak'
        fres = @(u,t) LC_residual_weak(mesh,ref_SBP,u,eqn,t);
    case 'strong'
        fres = @(u,t) LC_residual_strong(mesh,ref_SBP,u,eqn,t);
end

% Time step
% CFL = .05;   % Conservative CFL
% h = CFL*dx;  % This should technically include the wave speed
h = 0.005; 

% Solve PDE
nn = mesh.con'; nn = nn(:); xxyy = mesh.coord(nn,:); t0 = 0.00;
U0 = eqn.u(xxyy,t0);          % Initial condition
U = RK4(fres,t0,tf,h,U0,mesh);

% Plot solution
nn = mesh.con'; nn = nn(:); xxyy = mesh.coord(nn,:);
figure(2),plot_field(mesh,U) % Numerical solution
if (eqn.isexact)
    figure(3),plot_field(mesh,eqn.u(xxyy,tf)) % Exact solution (if applicable)
end

% Outputs:
% H-norm of the error (if applicable)
if (eqn.isexact)
    err = err_H(mesh,ref_SBP,eqn,tf,U)
end
end
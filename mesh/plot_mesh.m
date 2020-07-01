function h = plot_mesh(mesh,opt)
% PLOT_MESH plots linear mesh (quadratic info is ignored)
% INPUT
%   mesh: mesh structure
%   opt: optional structure; all fields below are optional
%     opt.number_nodes: if true, show node numbers
%     opt.number_elems: if true, show element numbers
%     opt.n_ref: number of refinement for curved meshes
% OUTPUT
%   h: handle to the patch object

% Copyright 2018 Masayuki Yano, University of Toronto

if (nargin < 2)
    opt = [];
end

colm = 'none'; % mesh color
cole = 'k'; % edge color

h = patch('vertices',mesh.vtxcoord,'faces',mesh.vtx(:,1:3),...
        'facecolor',colm,'edgecolor',cole);
hold on;


% plot nodes
if isfield(opt,'number_nodes') && (opt.number_nodes ~= false)
    hn = plot(mesh.coord(:,1),mesh.coord(:,2),'ob');
    %set(hn,'markersize',20,'markerfacecolor','w');
    set(hn,'markersize',10,'markerfacecolor','b');
    %htn = zeros(size(mesh.coord,1),1);
    %for i = 1:size(mesh.coord,1)
    %    htn(i) = text(mesh.coord(i,1),mesh.coord(i,2),num2str(i));
    %end
    %set(htn,'fontsize',16,'Color','blue','horizontalalignment','center','verticalalignment','middle');
end

% plot elements
if isfield(opt,'number_elems') && (opt.number_elems ~= false)
    con = mesh.con(:,1:3);
    ncon = size(con,1);
    xcon = reshape(mean(reshape(mesh.coord(con(:),:),[size(con),2]),2),[ncon,2]);
    hte = zeros(ncon,1);
    for i = 1:ncon
        hte(i) = text(xcon(i,1),xcon(i,2),num2str(i));
    end
    set(hte,'fontsize',16,'Color','r','edgecolor','r','backgroundcolor','w','horizontalalignment','center','verticalalignment','middle');
end

% mark origin
if isfield(opt,'mark_origin') && (opt.mark_origin ~= false)
    a = 0.3;
    con = mesh.con(:,1:3);
    ncon = size(con,1);
    xcon = reshape(mean(reshape(mesh.coord(con(:),:),[size(con),2]),2),[ncon,2]);
    xrv = a*xcon + (1-a)*mesh.coord(con(:,1),:);
    hrv = plot(xrv(:,1),xrv(:,2),'ok');
    set(hrv,'markersize',8,'markerfacecolor','k');
end
end


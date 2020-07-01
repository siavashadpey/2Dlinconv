function plot_field(mesh,U)

% Description: plots a 2D field 

N = 100;
nn = mesh.con'; nn = nn(:); xxyy = mesh.coord(nn,:);
zz = scatteredInterpolant(xxyy(:,1),xxyy(:,2),U);
xmin = min(mesh.coord(:,1)); xmax = max(mesh.coord(:,1));
ymin = min(mesh.coord(:,2)); ymax = max(mesh.coord(:,2));
xx = linspace(xmin,xmax,N); yy = linspace(ymin,ymax,N); 
[xgrid,ygrid] = meshgrid(xx,yy);
zgrid = zz(xgrid,ygrid);
surf(xgrid,ygrid,zgrid), colorbar
xlabel('x')
ylabel('y')
zlabel('solution')
end
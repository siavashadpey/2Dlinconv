function [shp,shpx] = monomial_tri(x)

N = size(x,1);

zero = zeros(N,1);
one = ones(N,1);

shp = [one x(:,1) x(:,2)];
shpx(:,:,1) = [zero one zero];
shpx(:,:,2) = [zero zero one];
end
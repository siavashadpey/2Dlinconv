function [shp,shpx] = monomial_line(x)

N = size(x,1);

zero = zeros(N,1);
one = ones(N,1);

shp = [one x];
shpx = [zero one];
end
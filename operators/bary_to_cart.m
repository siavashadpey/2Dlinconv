function ccoord = bary_to_cart(bcoord)

% Description: converts barycentric coordinates to cartesian coordinates. 

% Cartesian coords of triangle vertices
L = [-1  1 -1; 
     -1 -1  1];

n = size(bcoord,1);

ccoord = NaN(2,n);
for i=1:n
    ccoord(:,i) = bcoord(i,1)*L(:,1) + bcoord(i,2)*L(:,2) + bcoord(i,3)*L(:,3);
end

ccoord = ccoord';

end

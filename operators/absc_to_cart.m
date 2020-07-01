function [xy, w] = absc_to_cart(S)

% Description: converts cubature rule given in barycentric coords to
% cubature rules in cartesian coords. 

xy = [];
w  = [];

fields = fieldnames(S);
for i = 1:length(fields)
    field_nam = num2str(fields{i});
    field = S.(fields{i});
    if contains(field_nam,'s111')
        a = field.absc;
        w_s = field.w;
        b = [a(1), a(2), 1 - a(1) - a(2)];
        bcoord = perms(b); % 6 points
    elseif contains(field_nam,'s21')
        a = field.absc;
        w_s = field.w;
        b = [a, a, 1 - 2*a];
        bcoord = perms(b); % 3 points
    elseif contains(field_nam,'s3')
        a = field.absc;
        w_s = field.w;
        b = [a a a]; 
        bcoord = b; % 1 point
    else
        continue
    end
    bcoord = unique(bcoord,'rows','stable');
    xy_s = bary_to_cart(bcoord);
    n = size(xy_s,1);
    w_s = w_s*ones(n,1);
    
    xy = [xy; xy_s];
    w = [w; w_s];
end

vtx = [-1 1 -1; -1 -1 1]';
plot(vtx(:,1),vtx(:,2),'-',xy(:,1),xy(:,2),'o')
end
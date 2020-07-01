% Description: run convergence study tests

eqntyp = 'const_lin_conv';
spd_a.x = 0.5; spd_a.y = spd_a.x;
form = 'weak';
tf = 0.5; 

type = {'omega','omega','omega'}; M = length(type);
p = {1,2,3};
fol = 'data/';
n = (2:4)'; N = length(n);
h = 1./(4*n);
for j=1:M
    err = zeros(N,1);
    SBP.p = p{j}; SBP.type = type{j};
    for i=1:N
        err(i) = lin_conv_2D(SBP,h(i),eqntyp,spd_a,form,tf);
        err(i)
    end
    figure
    loglog(h,err,'o--')
    c = polyfit(log(h),log(err),1); 
    c(1) 
    fnam = strcat(fol,'conv_stud_',num2str(SBP.type),'_',num2str(SBP.p),'_fe_affine_grid.dat');
    header = strcat('# h , err. Conv rate: ', num2str(c(1)));
    print_data(fnam, header, [h err]);
end

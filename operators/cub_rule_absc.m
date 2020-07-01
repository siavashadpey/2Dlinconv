function [xy, w] = cub_rule_absc(p,type)

% Description: generates cubature rules in cartesian coord based off of
% cubature rules in barycentric coordinates 

type = lower(type);

switch type
    case 'diage'
        [S] = cub_rule_absc_diage(p);
end

[xy, w] = absc_to_cart(S);
end

function [S] = cub_rule_absc_diage(p)

switch p
    case 1
        S.s111.absc = [0, 1/2 - sqrt(3)/6];
        S.s111.w    = 1/3;
    case 2
        S.s111.absc = [0, 1/2 - sqrt(15)/10];
        S.s111.w    = 1/12;
        
        S.s21.absc  = 1/2;
        S.s21.w     = 1/5;
        
        S.s3.absc   = 1/3;
        S.s3.w      = 9/10;
    case 3
        S.s111_1.absc = [0, 0.330009478207572];
        S.s111_1.w    = 0.080913081365980;
        
        S.s111_2.absc = [0, 0.0694318442029737];
        S.s111_2.w    = 0.030198029745131;
        
        S.s111_3.absc = [0.1870738791912771, 0.5841571139756569];
        S.s111_3.w    = 0.222222222222222;
    case 4
        S.s111_1.absc = [0, 0.230765344947159];
        S.s111_1.w    = 0.041060919360858;
        
        S.s111_2.absc = [0, 0.046910077030668];
        S.s111_2.w    = 0.013202630162003;
        
        S.s21_1.absc  = 1/2;
        S.s21_1.w     = 0.037074169667900;
        
        S.s21_2.absc  = 0.4384239524408185;
        S.s21_2.w     = 0.249473464579547;
        
        S.s21_3.absc  = 0.1394337314154536;
        S.s21_3.w     = 0.210858659241689;
        
        S.s3.absc     = 1/3;
        S.s3.w        = 0.182199822395427;
        
end

end
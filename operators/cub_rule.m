function [xy, w, Fperm, Vperm] = cub_rule(p,type)

% Description: Holds cubature rules 

% TODO: add missing Fperm and Vperm matrices

type = lower(type);

switch type
    case 'gamma'
        [xy, w, Fperm, Vperm] = cub_rule_gamma(p);
    case 'omega'
        [xy, w, Fperm, Vperm] = cub_rule_omega(p);
    case 'diage'
        [xy, w, Fperm, Vperm] = cub_rule_diage(p);
    otherwise
        error('Unsupported SBP type. \n')
end

end

function [xy, w, Fperm, Vperm] = cub_rule_gamma(p)
% TODO: add Fperm and Vperm for p=3 and 4
switch p 
    case 1
        xy  = [-1.0               ,  1.0               , -1.0                 ;
               -1.0               , -1.0               ,  1.0               ]'; 
        w   = [ 0.6666666666666666,  0.6666666666666666,  0.6666666666666666]';
        Fperm = [1 2; 2 3; 3 1]';
        Vperm = [1 2 3; 2 3 1; 3 1 2]';
    case 2
        xy  = [-1.0                ,  1.0                , -1.0                ,  0.0                , 0.0                , -1.0                , -0.3333333333333333  ;
               -1.0                , -1.0                ,  1.0                , -1.0                , 0.0                ,  0.0                , -0.3333333333333333]';
        w   = [ 0.09999999999999999,  0.09999999999999999,  0.09999999999999999,  0.26666666666666666, 0.26666666666666666,  0.26666666666666666,  0.9000000000000002]';
        Fperm = [1 2 4; 2 3 5; 3 1 6]';
        Vperm = [1 2 3 4 5 6 7; 2 3 1 5 6 4 7; 3 1 2 6 4 5 7]';
    case 3 
        xy = [-1.0                ,  1.0                , -1.0                , -0.4130608881819197 ,  0.4130608881819197 ,  0.4130608881819197 , -0.4130608881819197 , -1.0                , -1.0               ,  -0.5853096486728182,  0.17061929734563638, -0.5853096486728182   ;
              -1.0                , -1.0                ,  1.0                , -1.0                , -1.0                , -0.4130608881819197 ,  0.4130608881819197 ,  0.4130608881819197 , -0.4130608881819197,  -0.5853096486728182, -0.5853096486728182 ,  0.17061929734563636]';
        w  = [ 0.02974582604964118,  0.02974582604964118,  0.02974582604964118,  0.09768336246810204,  0.09768336246810204,  0.09768336246810204,  0.09768336246810204,  0.09768336246810204, 0.09768336246810204,   0.4415541156808217,  0.4415541156808217 ,  0.4415541156808217 ]';
    case 4
        xy = [-1.0                 ,  1.0                 , -1.0                 , -0.5773502691896257,   0.0                ,  0.5773502691896257 ,  0.5773502691896257 , 0.0                , -0.5773502691896257 , -1.0                , -1.0                , -1.0                , -0.7384168123405102, -0.1504720765483788,  0.4768336246810203, -0.1504720765483788 , -0.7384168123405102 , -0.6990558469032424   ;
              -1.0                 , -1.0                 ,  1.0                 , -1.0               ,  -1.0                , -1.0                , -0.5773502691896257 , 0.0                ,  0.5773502691896257 ,  0.5773502691896257 ,  0.0                , -0.5773502691896257 , -0.7384168123405102, -0.6990558469032424, -0.7384168123405102, -0.15047207654837885,  0.47683362468102025, -0.15047207654837885]';
        w  = [ 0.012698412698412695,  0.012698412698412695,  0.012698412698412695,  0.04285714285714284,  0.05079365079365077,  0.04285714285714284,  0.04285714285714284, 0.05079365079365077,  0.04285714285714284,  0.04285714285714284,  0.05079365079365077,  0.04285714285714284,  0.2023354595827503,  0.3151248578775673,  0.2023354595827503,  0.3151248578775673 ,  0.2023354595827503 ,  0.3151248578775673 ]';
    otherwise
        error('Unsupported p degree for Gamma SBP cubature rule. \n')
end

end

function [xy, w, Fperm, Vperm] = cub_rule_omega(p)

switch p 
    case 1
        xy = [-0.6666666666666667,  0.3333333333333335, -0.6666666666666667  ;
              -0.6666666666666667, -0.6666666666666667,  0.3333333333333334]';
        w  = [ 0.6666666666666666,  0.6666666666666666,  0.6666666666666666]';
        Fperm = [1 2 3; 2 3 1; 3 1 2]';
    case 2
        xy = [-0.8168475729804585,  0.633695145960917 , -0.8168475729804585, -0.10810301816807022, -0.10810301816807022, -0.7837939636638596   ;
              -0.8168475729804585, -0.8168475729804585,  0.633695145960917 , -0.7837939636638596 , -0.10810301816807022, -0.10810301816807022]';
        w  = [ 0.2199034873106437,  0.2199034873106437,  0.2199034873106437, 0.44676317935602283 ,  0.44676317935602283,  0.44676317935602283]';
        Fperm = [1 2 3 4 5 6; 2 3 1 5 6 4; 3 1 2 6 4 5]'; 
    case 3
        xy = [-0.8613766937233731 ,  0.7227533874467464 , -0.8613766937233732 , -0.3773557414367168 ,  0.23519630088171353,  0.2351963008817135 , -0.3773557414367168 , -0.8578405594449967 , -0.8578405594449967 , -0.3333333333333333   ;
              -0.8613766937233733 , -0.8613766937233733 ,  0.7227533874467464 , -0.8578405594449967 , -0.8578405594449967 , -0.3773557414367168 ,  0.2351963008817135 ,  0.2351963008817135 , -0.3773557414367168 , -0.3333333333333333 ]';
        w =  [ 0.11550472674301035,  0.11550472674301035,  0.11550472674301035,  0.20924480696331949,  0.20924480696331949,  0.20924480696331949,  0.20924480696331949,  0.20924480696331949,  0.20924480696331949,  0.39801697799105223]';
        Fperm = [1 2 3 4 5 6 7 8 9 10; 2 3 1 6 7 8 9 4 5 10; 3 1 2 8 9 4 5 6 7 10]'; 
    case 4
        xy = [-0.9156687711811358  ,  0.8313375423622715  , -0.9156687711811358  , -0.5768758827238152 , -0.05141062176497946,  0.48091319998088583,  0.4809131999808858 , -0.05141062176497946, -0.5768758827238151 , -0.9040373172570707 , -0.8971787564700411 , -0.9040373172570706 , -0.5158280524810427,  0.031656104962085374, -0.5158280524810427    ;
              -0.9156687711811358  , -0.9156687711811358  ,  0.8313375423622715  , -0.9040373172570706 , -0.8971787564700411 , -0.9040373172570706 , -0.5768758827238152 , -0.05141062176497946,  0.4809131999808858 ,  0.4809131999808858 , -0.05141062176497946, -0.5768758827238152 , -0.5158280524810426, -0.5158280524810426  ,  0.031656104962085374]';
        w  = [ 0.045386157905236924,  0.045386157905236924,  0.045386157905236924,  0.11055758694624952,  0.145828414950907  ,  0.11055758694624952,  0.11055758694624952,  0.145828414950907  ,  0.11055758694624952,  0.11055758694624952,  0.145828414950907  ,  0.11055758694624952,  0.2543369199180238,  0.2543369199180238  ,  0.2543369199180238  ]';  
        Fperm = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; 2 3 1 7 8 9 10 11 12 4 5 6 14 15 13; 3 1 2 10 11 12 4 5 6 7 8 9 15 13 14]';
    otherwise
        error('Unsupported p degree for Omega SBP cubature rule. \n')
end
Vperm = Fperm;
end

function [xy, w, Fperm, Vperm] = cub_rule_diage(p)

switch p 
    case 1
        xy = [ -0.577350269189626   0.577350269189626   0.577350269189626  -0.577350269189626  -1.000000000000000   -1.000000000000000;
               -1.000000000000000  -1.000000000000000  -0.577350269189626   0.577350269189626   0.577350269189626   -0.577350269189626]';
        w  = [0.333333333333333    0.333333333333333   0.333333333333333   0.333333333333333   0.333333333333333   0.333333333333333]'; 
        Fperm = [1 2; 3 4; 5 6]';
        Vperm = [ 1 2 3 4 5 6; 3 4 5 6 1 2; 5 6 1 2 3 4]';
    case 2
        xy = [-0.774596669241483                   0  0.774596669241483  0.774596669241483                  0 -0.774596669241483  -1.000000000000000 -1.000000000000000 -1.000000000000000 -0.333333333333333
              -1.000000000000000  -1.000000000000000 -1.000000000000000 -0.774596669241483                  0  0.774596669241483   0.774596669241483                  0 -0.774596669241483 -0.333333333333333]';
        w  = [ 0.083333333333333   0.200000000000000  0.083333333333333  0.083333333333333  0.200000000000000  0.083333333333333   0.083333333333333  0.200000000000000  0.083333333333333  0.900000000000000]';
        Fperm = [1 2 3; 4 5 6; 7 8 9]';
        Vperm = [1 2 3 4 5 6 7 8 9 10; 4 5 6 7 8 9 1 2 3 10; 7 8 9 1 2 3 4 5 6 10]';
    case 3
        xy = [-0.861136311594053 ,-0.339981043584856 , 0.339981043584856 , 0.861136311594053 , 0.861136311594053 , 0.339981043584856 ,-0.339981043584856 ,-0.861136311594053 ,-1                 ,-1                 ,-1                 ,-1                 ,-0.542461986333868, 0.168314227951314, 0.168314227951314,-0.542461986333868,-0.625852241617446,-0.625852241617446  ;
              -1                 ,-1                 ,-1                 ,-1                 ,-0.861136311594053 ,-0.339981043584856 , 0.339981043584856 , 0.861136311594053 , 0.861136311594053 , 0.339981043584856 ,-0.339981043584856 ,-0.861136311594053 ,-0.625852241617446,-0.625852241617446,-0.542461986333868, 0.168314227951314, 0.168314227951314,-0.542461986333868]';
        w  = [ 0.0301980297451310, 0.0809130813659800, 0.0809130813659800, 0.0301980297451310, 0.0301980297451310, 0.0809130813659800, 0.0809130813659800, 0.0301980297451310, 0.0301980297451310, 0.0809130813659800, 0.0809130813659800, 0.0301980297451310, 0.222222222222222, 0.222222222222222, 0.222222222222222, 0.222222222222222, 0.222222222222222, 0.222222222222222]';
        Fperm = [1 2 3 4; 5 6 7 8; 9 10 11 12]';
        Vperm = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18; 5 6 7 8 9 10 11 12 1 2 3 4 15 16 17 18 13 14; 9 10 11 12 1 2 3 4 5 6 7 8 17 18 13 14 15 16]'; 
    case 4
        xy = [-0.906179845938664 ,-0.538469310105682 , 0                 , 0.538469310105682 , 0.906179845938664 , 0.906179845938664 , 0.538469310105682 ,0                 ,-0.538469310105682 ,-0.906179845938664 ,-1                 ,-1                 ,-1                 ,-1                 ,-1                 ,-0.721132537169093,-0.123152095118363, 0.442265074338186,-0.123152095118363,-0.721132537169093,-0.753695809763274,-0.333333333333333  ;
              -1                 ,-1                 ,-1                 ,-1                 ,-1                 ,-0.906179845938664 ,-0.538469310105682 ,0                 , 0.538469310105682 , 0.906179845938664 , 0.906179845938664 , 0.538469310105682 , 0                 ,-0.538469310105682 ,-0.906179845938664 ,-0.721132537169093,-0.753695809763274,-0.721132537169093,-0.123152095118363, 0.442265074338186,-0.123152095118363,-0.333333333333333]';
        w  = [ 0.0132026301620030, 0.0410609193608580, 0.0370741696679000, 0.0410609193608580, 0.0132026301620030, 0.0132026301620030, 0.0410609193608580,0.0370741696679000, 0.0410609193608580, 0.0132026301620030, 0.0132026301620030, 0.0410609193608580, 0.0370741696679000, 0.0410609193608580, 0.0132026301620030, 0.210858659241689, 0.249473464579547, 0.210858659241689, 0.249473464579547, 0.210858659241689, 0.249473464579547, 0.182199822395427]';
        Fperm = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15]';
        Vperm = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; 6 7 8 9 10 11 12 13 14 15 1 2 3 4 5 18 19 20 21 16 17 22; 11 12 13 14 15 1 2 3 4 5 6 7 8 9 10 20 21 16 17 18 19 22]';
    otherwise 
        error('Unsupported p degree for diagonal E SBP cubature rule. \n')
end

end


function setup()
dirs = {'distmesh','mesh','ref','main','residual','solver',...
        'operators','eqn','pprocess','run'};
d0=fileparts([pwd,filesep]);
d0=[d0,filesep];
for k = 1:length(dirs)
addpath([d0,dirs{k}]);
end
fprintf('Folders added to path.\n')



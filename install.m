function install
cocodir = fileparts(mfilename('fullpath'));
addpath(fullfile(cocodir, 'toolbox'));
addpath(fullfile(cocodir, 'toolbox', 'MyPHI'));
addpath(fullfile(cocodir, 'toolbox', 'Initial_solution'));
addpath(fullfile(cocodir, 'toolbox', 'MyMethods_forQP'));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The COCO setup function must also be run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function startup
% cocodir = fileparts(fileparts(mfilename('fullpath')));
% addpath(fullfile(cocodir, 'toolbox'));
% cp = coco_path;
% addpath(cp{:});
% end
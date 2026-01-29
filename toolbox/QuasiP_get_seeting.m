function data = QuasiP_get_seeting(prob,tbid,data)

% QUASIP_GET_SETTINGS   Read 'QuasiP' toolbox instance settings.
%
% Merge user-supplied toolbox settings with default values.
%
% DATA = QUASIP_GET_SETTINGS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data strcture.

defaults.order   = 1;

par_opts.parallel = false;    % integration in parallel computation
par_opts.ncores   = 1;        % number of cores in parallel computation
defaults.par_opts = par_opts;

copts = coco_get(prob, tbid);
copts = coco_merge(defaults, copts);

data.order     = copts.order;

assert(isfield(copts, 'ode_opts'), ...
  '%s: The d and e (number of base frequencies) is not defined');
data.ode_opts     = copts.ode_opts;

assert(isfield(copts, 'DOF'), ...
  '%s: The number of DOF is not defined');
data.DOF     = copts.DOF;

if ~isfield(copts, 'PHI_INITIAL') % defaults.method: VCHB
    
    methods = repmat({'HB'}, 1, data.ode_opts.d);
    k = repmat({1:1:5}, 1, data.ode_opts.d);
    U = repmat({2*5+1}, 1, data.ode_opts.d);
    S = repmat({2^5}, 1, data.ode_opts.d);
    S_w = repmat({2^8}, 1, data.ode_opts.d);
    
    copts.PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ... 
    'd',       data.ode_opts.d,...
    'e',       data.ode_opts.e,...
    'method',  { methods   }, ...   
    'isCO1',   0,...                     
    'k',       { k    }, ...        
    'U',       { U    }, ...    
    'S',       { S    }, ...
    'S_wave',  { S_w  });
end
PHI = COM_PHI(copts.PHI_INITIAL);
data.PHI     = PHI;

if copts.par_opts.parallel
    % activate parallel computing
    activate_parallel(copts.par_opts.ncores)
end

end

function activate_parallel(varargin)
%ACTIVATE_PARALLEL This function starts parallel pool for the first time
%and detects the number of cores avaliable

try    
    pp1 = gcp('nocreate');
    if isempty(pp1)
        h = helpdlg('Starting parallel pool for the first time and detecting number of available cores.', 'Info');
        disp('----------------------------------------------------------------')
        disp('Starting parallel pool for the first time and detecting number')
        disp('of available cores.')
        disp('----------------------------------------------------------------')
        defaultProfile = parallel.defaultClusterProfile;
        myCluster = parcluster(defaultProfile);
        if numel(varargin)>=1 && isnumeric(varargin{1})
            parpool(myCluster,varargin{1});
        else
            parpool(myCluster);
        end
        pp2 = gcp('nocreate');
        cpuNum =pp2.NumWorkers;
        save('cluster_info.mat','cpuNum')

        if isvalid(h)
            close(h)
        end
    end
catch ME
    rethrow(ME)
end

end



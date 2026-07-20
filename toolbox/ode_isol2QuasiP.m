function prob = ode_isol2QuasiP(prob, oid, varargin)
% ODE_ISOL2QUASIP   Append 'QuasiP' instance constructed from initial data.
%
% Parse input sequence to construct toolbox data and initial solution guess
% and use this to construct an instance of 'QuasiP'.
%
% PROB     = ODE_ISOL2QUASIP(PROB, OID, @F, @DFDX, X_MF, OM)
% VARARGIN = { @DFDP|[], @L|[], PNAMES|[], P0|[] }
%
% PROB   - Continuation problem structure.
% OID    - Object instance identifier (string).

% ODEs:  - L(p) = F(x,[dx],p,theta)
%          =A(p)ddx + B(p)dx + C(p)x

% @F     - Function handle to Nonlinear part and Excitation.
% @DFDX  - Function handle to Jacobian w.r.t. x and [dx].
% @DFDP  - Optional function handle to Jacobian w.r.t. p.
%          if F(x,[dx],theta).
% @L     - Optional function handle to Linear part and Jacobian w.r.t. p.
%          if B = eye.
% X_MF   - Initial amplitude of x in multi-frequency domain
% OM     - Initial base frequeincies of x. (also be defined as new continuation parameters p)
% PNAMES - Optional string label or cell array of string labels for
%          continuation parameters tracking p.
% P0     - Optional initial p.
tbid = coco_get_id(oid, 'QuasiP');
data = struct();
data = QuasiP_get_seeting(prob, tbid, data);
%----------------------------------------------------------
grammar   = 'F DFDX [DFDP [L [DFDTH]]] X_MF OM [PNAMES [P0]]';
args_spec = {
    'F',          '',        '@',      'fhan',   [],  'read',  {}
    'DFDX',       '',        '@',   'dfdxhan',   [],  'read',  {}
    'DFDP',       '',     '@|[]',   'dfdphan',   [],  'read',  {}
    'L',          '',     '@|[]',      'lhan',   [],  'read',  {}
    'DFDTH',      '',     '@|[]',  'dfdthhan',   [],  'read',  {}    
    'X_MF',       '',    '[num]',      'X_MF',   [],  'read',  {}
    'OM',         '',    '[num]',        'OM',   [],  'read',  {}
    'PNAMES', 'cell',    '{str}',    'pnames',   {},  'read',  {}
    'P0',         '',    '[num]',        'p0',   [],  'read',  {}
    };
[args, ~] = coco_parse(grammar, args_spec, [], varargin{:});
%----------------------------------------------------------
data.f   = args.fhan;
data.fx  = args.dfdxhan;
data.fp  = args.dfdphan;
data.l   = args.lhan;
data.fth = args.dfdthhan;
X_MF     = args.X_MF;
OM       = args.OM;
p0       = args.p0;
%----------------------------------------------------------
data = quasip_post_process(tbid, data, args, p0);
quasip_print_equation_table(data);
%----------------------------------------------------------
sol = QuasiP_init_sol(data,X_MF,OM,p0);
QuasiP_arg_check(tbid, data, X_MF, OM, p0);
data = QuasiP_init_data(data, OM, p0);
prob = ode_construct_QuasiP(prob, tbid, data, sol);
end
%==========================================================================
%==========================================================================
%==========================================================================
function data = quasip_post_process(tbid, data, args, p0)
if strcmp(tbid, 'QuasiP')
    prefix = '';
else
    parts = strsplit(tbid, '.');
    prefix = [parts{1}, '.'];
end
data.my_tbid = tbid;
data.prefix  = prefix;
%----------------------------------------------------------
if ~isempty(args.pnames) && iscell(args.pnames)
    for i = 1:length(args.pnames)
        data.pnames{i} = [prefix, args.pnames{i}];
    end
else
    data.pnames = args.pnames;
end
%----------------------------------------------------------
data.kOMnames = cell(1, data.ode_opts.e);
for ii = 1 : data.ode_opts.e
    data.kOMnames{ii} = sprintf([prefix 'k_OM%d'], ii);
end
%----------------------------------------------------------
nUOM = data.ode_opts.d - data.ode_opts.e;
data.uOMnames = cell(1, nUOM);
for ii = 1 : nUOM
    data.uOMnames{ii} = sprintf([prefix 'u_om%d'], data.ode_opts.e + ii);
end
%----------------------------------------------------------
if isempty(data.l)
    data.ODE = 1;
    data.isL = false;
else
    [isABC,isABCwithp] = data.l(p0);
    if isABC(1) == 1
        data.ODE = 2;
    else
        data.ODE = 1;
    end
    data.isL        = true;
    data.isABC      = isABC;
    data.isABCwithp = isABCwithp;
end
%----------------------------------------------------------
data.isp      = ~isempty(p0);
data.isfwithp = ~isempty(data.fp);
data.isfwithu = data.ode_opts.isfwithu;
%----------------------------------------------------------
if data.stable_opts.issta
    assert(~(data.ODE == 1 && data.isfwithu), ...
        'Cannot determine the stability of dx = f(x,dx)');
end
switch data.PHI.parametrization_tori
    case 'MT'
        switch data.PHI.type
            case 'SHT'
                data.PHI = MT_SHT_COM_PHI_Supplement(data);
        end
end
end
%==========================================================================
%==========================================================================
%==========================================================================
function quasip_print_equation_table(data)
if data.isL
    terms = [];names = {'A', 'B', 'C'};exps = {'*ddx', '*dx', '*x'};
    for i = 1:3
        if data.isABC(i)
            term = names{i};
            if data.isABCwithp(i)
                term = [term '(p)'];
            end
            term = [term exps{i}];
            if isempty(terms)
                terms = term;
            else
                terms = [terms ' + ' term];
            end
        end
    end
    terms = [terms ' = f(x'];
else
    terms = 'dx = f(x';
end
%----------------------------------------------------------
if data.isfwithu
    terms = [terms ',dx'];
end
if data.isfwithp
    terms = [terms ',p'];
end
for ii = 1:data.ode_opts.e
    terms = [terms, sprintf(',theta%d', ii)];
end
terms = [terms, ')'];
%----------------------------------------------------------
headers = {'N','d','e','pnames','kOMs','uoms','para_tori','type_method','period','L','pS'};
switch data.PHI.parametrization_tori
    case{'MP'}
        values  = {
            num2str(data.ode_opts.n), ...
            num2str(data.ode_opts.d), ...
            num2str(data.ode_opts.e), ...
            strjoin(data.pnames, '+'), ...
            strjoin(data.kOMnames, '+'), ...
            strjoin(data.uOMnames, '+'), ...
            data.PHI.parametrization_tori, ...
            strcat(string(data.PHI.type), '_', strjoin(data.PHI.method,'+'), '_(', strjoin(string(data.PHI.OmOd(:)),'>-'), ')'), ...
            strjoin(string(data.PHI.P), '+'), ...
            num2str(data.PHI.L), ...
            strjoin(string(data.PHI.pS), '+')
            };
    case{'MT'}
        switch data.PHI.type
            case{'DM'}

            case{'SHT'}
                pS = cell(1,data.ode_opts.d);
                if data.PHI.IsFixedStep
                    pS{data.PHI.Omt}   = data.PHI.P{data.PHI.Omt}*data.PHI.Nsteps;
                    for ii = 1:data.ode_opts.d-1
                        pS{data.PHI.Omnot(ii)} = data.PHI.pSnot{ii};
                    end
                else
                    pS{data.PHI.Omt}   = 'UN';
                    for ii = 1:data.ode_opts.d-1
                        pS{data.PHI.Omnot(ii)} = data.PHI.pSnot{ii};
                    end
                end
                values  = {
                    num2str(data.PHI.N), ...
                    num2str(data.ode_opts.d), ...
                    num2str(data.ode_opts.e), ...
                    strjoin(data.pnames, '+'), ...
                    strjoin(data.kOMnames, '+'), ...
                    strjoin(data.uOMnames, '+'), ...
                    data.PHI.parametrization_tori, ...
                    [data.PHI.type '_' data.PHI.method ' |' func2str(data.PHI.ODEsolver) ' |Omt=' num2str(data.PHI.Omt) ' |M=' num2str(data.PHI.M_SHT)], ...
                    strjoin(string(data.PHI.P), '+'), ...
                    num2str(data.PHI.L), ...
                    strjoin(string(pS), '+')
                    };
        end
end
%----------------------------------------------------------
values_char = cellfun(@(v) char(v), values, 'UniformOutput', false);
col_widths = cellfun(@(h,v) max(strlength(h), strlength(v))+4, headers, values_char);
col_widths(end) = col_widths(end)-4;
total_width = sum(col_widths);
line = repmat('=', 1, total_width);
term_str = [' ' terms ' '];
fill_len = total_width  - length(term_str);
left_fill = repmat('>', 1, ceil(fill_len/2));
right_fill = repmat('<', 1, floor(fill_len/2));
title_line = [left_fill term_str right_fill];
%----------------------------------------------------------
disp(line);
disp(title_line);
disp(line);
fmt = strjoin(compose('%%-%ds', col_widths), '');
fprintf([fmt '\n'], headers{:});
fprintf([fmt '\n'], values{:});
disp(line);
end
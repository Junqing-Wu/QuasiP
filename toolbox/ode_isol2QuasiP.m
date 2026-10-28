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
% @DFDTH - Optional function handle to Jacobian w.r.t. theta
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
pnames   = args.pnames;
%----------------------------------------------------------
data = QuasiP_pnames_update(tbid, data, pnames);
data = QuasiP_system_identify(data,p0);
%----------------------------------------------------------
if data.stable_opts.issta
    assert(~(data.ODE == 1 && data.isfwithu), ...
        'Cannot determine the stability of dx = f(x,dx)');
end
switch data.PHI.param
    case 'MT'
        switch data.PHI.method
            case 'SHT'
                if data.ode_opts.e > 0
                assert(~isempty(data.fth), ...
                    '%s: e>0, data.fth is empty!', ...
                    tbid);
                end
                data.PHI = MT_SHT_COM_PHI_Supplement(data);
        end
end
QuasiP_print_equation_table(data);
%----------------------------------------------------------
sol = QuasiP_init_sol(data,X_MF,OM,p0);
QuasiP_arg_check(tbid, data, X_MF, OM, p0);
data = QuasiP_init_data(data, OM, p0);
prob = ode_construct_QuasiP(prob, tbid, data, sol);
end
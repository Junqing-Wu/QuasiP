function prob = ode_isol2QuasiP(prob, oid, varargin)
% ODE_ISOL2QUASIP   Append 'QuasiP' instance constructed from initial data.
%
% Parse input sequence to construct toolbox data and initial solution guess
% and use this to construct an instance of 'QuasiP'.
%
% PROB     = ODE_ISOL2FORWARD(PROB, OID, @F @DFDX @DFDP T0 T X0 X1 PNAMES P0)
% VARARGIN = { @DFDT|[], @DFDXDX|[], @DFDXDP|[], @DFDPDP|[] }
%
% PROB   - Continuation problem structure.
% OID    - Object instance identifier (string).
% @F     - Function handle to vector field.
% @DFDX  - Optional function handle to Jacobian w.r.t. problem variables.
% @DFDP  - Optional function handle to Jacobian w.r.t. problem parameters.
% T0     - Initial time.
% T      - Time period.
% X0     - Initial state.
% X1     - Final state.
% PNAMES - Optional string label or cell array of string labels for
%          continuation parameters tracking problem parameters.
% P0     - Initial solution guess for problem parameters.

tbid = coco_get_id(oid, 'QuasiP'); % Create toolbox instance identifier

data = struct();
data = QuasiP_get_seeting(prob, tbid, data);       % Get toolbox settings
switch data.order
    case 1
        grammar   = 'B DBDP F DFDX DFDP [DFDT [DFDXDX [DFDXDP [DFDPDP]]]] X_MF OM [PNAMES] P0';
        args_spec = {
               'B',     '',     '@',      'bhan', [], 'read', {}
            'DBDP',     '',     '@',   'dbdphan', [], 'read', {}            
               'F',     '',     '@',      'fhan', [], 'read', {}
            'DFDX',     '',     '@',   'dfdxhan', [], 'read', {}
            'DFDP',     '',     '@',   'dfdphan', [], 'read', {}
            'DFDT',     '',  '@|[]',   'dfdthan', [], 'read', {}
          'DFDXDX',     '',  '@|[]', 'dfdxdxhan', [], 'read', {}
          'DFDXDP',     '',  '@|[]', 'dfdxdphan', [], 'read', {}
          'DFDPDP',     '',  '@|[]', 'dfdpdphan', [], 'read', {}
            'X_MF',     '', '[num]',      'X_MF', [], 'read', {}
              'OM',     '', '[num]',        'OM', [], 'read', {}
          'PNAMES', 'cell', '{str}',    'pnames', {}, 'read', {}
              'P0',     '', '[num]',        'p0', [], 'read', {}
          };
        [args, ~] = coco_parse(grammar, args_spec, [], varargin{:});
        data.b   = args.bhan;
        data.bp  = args.dbdphan;
        data.f   = args.fhan;
        data.fx  = args.dfdxhan;
        data.fp  = args.dfdphan;
        data.ft  = args.dfdthan;
        data.fxx = args.dfdxdxhan;
        data.fxp = args.dfdxdphan;
        data.fpp = args.dfdpdphan;
        
    case 2
        
            
    otherwise
        error('The order of ODE should be either 1 or 2.');
end

X_MF  = args.X_MF;
OM    = args.OM;
p0    = args.p0;
data.pnames = args.pnames;
data.nL = data.DOF*data.PHI.L;

QuasiP_arg_check(tbid, data, X_MF, OM, p0);          % Validate input
data = QuasiP_init_data(data, OM, p0);               % Build toolbox data
sol  = QuasiP_init_sol(data,X_MF,OM,p0);             % Build initial solution guess
prob = ode_construct_QuasiP(prob, tbid, data, sol);  % Append continuation problem

end

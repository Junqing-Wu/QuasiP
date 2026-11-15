function prob = ode_isol2QuasiP(prob, oid, varargin)
%ODE_ISOL2QUASIP  Start continuation of a d-torus from an initial guess.
%
% This is the public Starter/Generator entry point of QuasiP. It reads the
% model callbacks and Xk_MF, identifies the ODE representation, constructs
% PHI for MP/m-VCF or MT/FSE-SHT, maps Xk_MF to C, and registers the NAE
% residual and OPC rows with COCO.
%
% The model is l(ddx,dx,x,p)=f(x,dx,theta_1^e,p), with
% l=A(p)*ddx+B(p)*dx+C(p)*x. If A is absent, the first-order form is used.
% The torus is 2*P(i)*pi-periodic in theta_i; P(i) is the period-multiplication
% number that permits torus-doubling and higher torus-period branches.
% The COCO unknown vector is u=[C(:); Om(:); p(:)], where C is the NL-vector
% of undetermined coefficients, Om is omega_1^d, and p is the model vector.
%
% PROB = ODE_ISOL2QUASIP(PROB, OID, @F, @DFDX, [@DFDP], [@L], [@DFDTH], Xk_MF, OM, [PNAMES], [P0])
%
% On input:
%   PROB   : COCO problem structure containing QuasiP options.
%   OID    : Object identifier; the instance is coco_get_id(OID,'QuasiP').
%   @F     : Callback for f(x,dx,theta,p), or f(x,theta,p) when isfwithu=false.
%   @DFDX  : Jacobian of f with respect to x and, when present, dx.
%   @DFDP  : Optional Jacobian of f with respect to p; [] when unavailable.
%   @L     : Optional callback returning A, B, C and their p-derivatives.
%   @DFDTH : Optional derivative with respect to excitation coordinates.
%   Xk_MF  : Struct with k_wave=k_tilde and Xc=X_c(k_tilde). Rows of k_wave
%            are harmonic combinations compatible with P(i), k(i), and PHI.
%   OM     : Initial omega_1^d; the first e entries are known excitation
%            frequencies and the remaining entries are unknown response frequencies.
%   PNAMES : Names of the model parameters p exposed to COCO.
%   P0     : Initial model parameters; use [] when no model parameter is present.
%
% On output:
%   PROB   : COCO problem with R(C,Om,p)=0, OPC rows, monitors, and slots.
%
% See also: QUASIP_INIT_SOL, ODE_CONSTRUCT_QUASIP, QUASIP_GET_SEETING.

tbid = coco_get_id(oid, 'QuasiP');
data = struct();
data = QuasiP_get_seeting(prob, tbid, data);
%----------------------------------------------------------
grammar   = 'F DFDX [DFDP [L [DFDTH]]] Xk_MF OM [PNAMES [P0]]';
args_spec = {
             'F',          '',        '@',      'fhan',   [],  'read',  {}
             'DFDX',       '',        '@',   'dfdxhan',   [],  'read',  {}
             'DFDP',       '',     '@|[]',   'dfdphan',   [],  'read',  {}
             'L',          '',     '@|[]',      'lhan',   [],  'read',  {}
             'DFDTH',      '',     '@|[]',  'dfdthhan',   [],  'read',  {}    
             'Xk_MF',      '',     'stct',     'Xk_MF',   [],  'read',  {}
             'OM',         '',    '[num]',        'OM',   [],  'read',  {}
             'PNAMES', 'cell',    '{str}',    'pnames',   {},  'read',  {}
             'P0',         '',    '[num]',        'p0',   [],  'read',  {}
            };
% Parse the callback handles, initial coefficients, frequencies, and parameters
% according to the calling sequence documented above.
[args, ~] = coco_parse(grammar, args_spec, [], varargin{:});
%----------------------------------------------------------
% Copy parsed callbacks and initial data into the internal DATA structure.
% Subsequent private residual functions access only DATA and the solution vector u.
data.f   = args.fhan;
data.fx  = args.dfdxhan;
data.fp  = args.dfdphan;
data.l   = args.lhan;
data.fth = args.dfdthhan;
Xk_MF    = args.Xk_MF;
OM       = args.OM;
p0       = args.p0;
pnames   = args.pnames;
%----------------------------------------------------------
% Build COCO-visible names for model parameters and base frequencies.
data = QuasiP_pnames_update(tbid, data, pnames);
% Identify first- or second-order form and record which Jacobians are available.
data = QuasiP_system_identify(data,p0);
%----------------------------------------------------------
% Stability post-processing requires compatible model callbacks.
if data.sta_opts.issta
    assert(~(data.ODE == 1 && data.isfwithu), ...
        'Cannot determine the stability of dx = f(x,dx)');
end
% Add formulation-specific data. MT-SHT needs supplementary trajectory and
% phase-shift information before its residual can be evaluated.
switch data.PHI.param
    case 'MT'
        switch data.PHI.method
            case 'SHT'
                if data.ode_opts.e > 0
                assert(~isempty(data.fth),'%s: e>0, data.fth is empty!', tbid);
                end
                data.PHI = MT_SHT_COM_PHI_Supplement(data);
        end
end
% Print the assembled governing-equation summary for reproducibility.
QuasiP_print_equation_table(data);
%----------------------------------------------------------
% Convert the user initial guess into the internal unknown vector
% sol.u = [C(:); Om(:); p0(:)].
sol = QuasiP_init_sol(data,Xk_MF,OM,p0);
% Check dimensions, frequency ordering, and consistency of the initial guess.
QuasiP_arg_check(tbid, data, Xk_MF, OM, p0);
% Store index ranges for C, Om, p, excitation frequencies, and unknown frequencies.
data = QuasiP_init_data(data, OM, p0);
% Register the residual, phase condition, monitors, and post-processing with COCO.
prob = ode_construct_QuasiP(prob, tbid, data, sol);
end
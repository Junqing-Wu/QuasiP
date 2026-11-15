function prob = ode_construct_QuasiP(prob, tbid, data, sol)
%ODE_CONSTRUCT_QUASIP  Register the QuasiP Generator and Processor in COCO.
%
% Add the NAE residual R(C,Om,p)=0, add d-e OPC equations when response
% frequencies are unknown, expose p, known frequencies kOMs, and unknown
% frequencies uOMs as COCO parameters, and install the save/post-processing
% slots. This is the QuasiP counterpart of the COCO po construction routine.
%
% On input:
%   prob : COCO problem structure.
%   tbid : QuasiP instance identifier.
%   data : Initialized Generator/Processor data.
%   sol  : Starter solution with sol.u=[C(:); Om(:); p(:)] and optional sol.t0.
%
% On output:
%   prob : COCO problem containing the registered QuasiP equations.
%
% A nonempty sol.t0 is passed to COCO as the branch-switching tangent.
has_t0 = isfield(sol, 't0') && ~isempty(sol.t0);
if data.PHI.isvar
    if ~has_t0
        prob = coco_add_func(prob, tbid, @QuasiP, data, 'zero', ...
          'u0', sol.u, 'F+dF');
    else
        prob = coco_add_func(prob, tbid, @QuasiP, data, 'zero', ...
          'u0', sol.u, 't0', sol.t0, 'F+dF');
    end
else
    if ~has_t0
        prob = coco_add_func(prob, tbid, @QuasiP, data, 'zero', ...
          'u0', sol.u);
    else
        prob = coco_add_func(prob, tbid, @QuasiP, data, 'zero', ...
          'u0', sol.u, 't0', sol.t0);
    end
end
%----------------------------------------------------------
% Add one orthogonality phase condition for every unknown response frequency.
% Expose the unknown response frequencies.
if data.ode_opts.d > data.ode_opts.e
    pcfid = coco_get_id(tbid, 'OPC');
    uidx   = coco_get_func_data(prob, tbid, 'uidx');
    prob  = coco_add_func(prob, pcfid, @QuasiP_PC, data, 'zero', 'uidx', uidx, 'F+dF');
end
%----------------------------------------------------------
% Expose model parameters as COCO continuation parameters.
if ~isempty(data.pnames)
    pfid = coco_get_id(tbid, 'pars');
    uidx = coco_get_func_data(prob, tbid, 'uidx');
    prob = coco_add_pars(prob, pfid, uidx(data.p_idx), data.pnames);
end
%----------------------------------------------------------
% Expose the prescribed excitation frequencies.
if data.ode_opts.e > 0
    komfid = coco_get_id(tbid, 'kOMs');
    uidx   = coco_get_func_data(prob, tbid, 'uidx');
    prob   = coco_add_pars(prob, komfid, uidx(data.kOM_idx), data.kOMnames);
end
%----------------------------------------------------------
% Add one orthogonality phase condition for every unknown response frequency.
% Expose the unknown response frequencies.
if data.ode_opts.d > data.ode_opts.e
    uomfid = coco_get_id(tbid, 'uOMs');
    uidx   = coco_get_func_data(prob, tbid, 'uidx');
    prob   = coco_add_pars(prob, uomfid, uidx(data.uOM_idx), data.uOMnames);
end
%----------------------------------------------------------
if data.sta_opts.issta && data.ode_opts.d == 1
    prob = QuasiP_d1_test(prob, tbid, data);
end
%----------------------------------------------------------
prob     = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');
PostPtid = coco_get_id(tbid, 'PostP');
prob     = coco_add_slot(prob, PostPtid, @PostP_Amp_Sta_bddat, data, 'bddat');
%----------------------------------------------------------
end
function prob = ode_construct_QuasiP(prob, tbid, data, sol)
% ODE_CONSTRUCT_QUASIP  Append an instance of 'QuasiP' to problem.
%
% Add R(C,OM,P)=0 and PHASE conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters.
%
% PROB = ODE_CONSTRUCT_QUASIP(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.   
%----------------------------------------------------------
switch data.PHI.parametrization_tori
    case{'MP'}
        switch data.PHI.type
            case{'mVCF'}
                prob = coco_add_func(prob, tbid, @MP_mVCF_QuasiP, data, 'zero', 'u0', sol.u, 'F+dF');
            case{'sCCF'}
        end
    case{'MT'}
        switch data.PHI.type
            case{'DM'}

            case{'SHT'}
                if data.PHI.isvar
                   prob = coco_add_func(prob, tbid, @MT_SHT_QuasiP, data, 'zero', 'u0', sol.u, 'F+dF');  
                else
                   prob = coco_add_func(prob, tbid, @MT_SHT_QuasiP, data, 'zero', 'u0', sol.u);  
                end
        end
end
if data.ode_opts.d > data.ode_opts.e
    pcfid = coco_get_id(tbid, 'OPC');
    uidx   = coco_get_func_data(prob, tbid, 'uidx');
    prob  = coco_add_func(prob, pcfid, @QuasiP_PC, data, 'zero', 'uidx', uidx, 'F+dF');
end
%----------------------------------------------------------
if ~isempty(data.pnames)
    pfid = coco_get_id(tbid, 'pars');
    uidx = coco_get_func_data(prob, tbid, 'uidx');
    prob = coco_add_pars(prob, pfid, uidx(data.p_idx), data.pnames);
end
%----------------------------------------------------------
if data.ode_opts.e > 0
    komfid = coco_get_id(tbid, 'kOMs');
    uidx   = coco_get_func_data(prob, tbid, 'uidx');
    prob   = coco_add_pars(prob, komfid, uidx(data.kOM_idx), data.kOMnames);
end
%----------------------------------------------------------
if data.ode_opts.d > data.ode_opts.e
    uomfid = coco_get_id(tbid, 'uOMs');
    uidx   = coco_get_func_data(prob, tbid, 'uidx');
    prob   = coco_add_pars(prob, uomfid, uidx(data.uOM_idx), data.uOMnames);
end
%----------------------------------------------------------
prob     = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');
PostPtid = coco_get_id(tbid, 'PostP');
prob     = coco_add_slot(prob, PostPtid, @PostP_forQP, data, 'bddat') ;
%----------------------------------------------------------
end




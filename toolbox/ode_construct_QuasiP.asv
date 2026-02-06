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

switch data.PHI.parametrization_tori
    case{'MP'}
        switch data.PHI.type
            case{'mVCF'}
                prob = coco_add_func(prob, tbid, @MP_mVCF_QuasiP, data, 'zero', 'u0', sol.u);
            case{'sCCF'}

        end
    case{'MT'}
        switch data.PHI.type
            case{'DM'}

            case{'SHOOT'}

        end
end

if ~isempty(data.pnames)
    pfid = coco_get_id(tbid, 'pars');
    uidx = coco_get_func_data(prob, tbid, 'uidx');
    prob = coco_add_pars(prob, pfid, uidx(data.p_idx), data.pnames);
end

if data.amp_opts.amp
    amid = coco_get_id(tbid, 'amp');
    uidx = coco_get_func_data(prob, tbid, 'uidx');
    outdof = data.amp_opts.outdof;
    numoutdof    = numel(outdof);
    ampNames     = cell(1, numoutdof);
    for k = 1:numoutdof
        ampNames{k} = strcat('amp',num2str(outdof(k)));
    end
    switch data.PHI.parametrization_tori
        case{'MP'}
            switch data.PHI.type
                case{'mVCF'}
                    prob = coco_add_func(prob, amid, @MP_mVCF_amp, data, 'regular', ampNames, 'uidx', uidx);
                case{'sCCF'}

            end
        case{'MT'}
            switch data.PHI.type
                case{'DM'}

                case{'SHOOT'}

            end
    end
end


prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end




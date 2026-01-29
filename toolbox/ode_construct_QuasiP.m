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

if data.ode_opts.d > data.ode_opts.e
    kk = data.ode_opts.e+1;
    for ii = 1:data.ode_opts.d-data.ode_opts.e
        PC_id = coco_get_id(tbid, sprintf('PC_%d', kk));
        uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
        data.PCnumber = kk;
        prob = coco_add_func(prob, PC_id, @QuasiP_PC, data, 'zero','uidx',uidx(data.C_idx));
        kk = kk+1;
    end
end

if ~isempty(data.pnames)
  pfid = coco_get_id(tbid, 'pars');
  uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
  prob = coco_add_pars(prob, pfid, uidx(data.p_idx), data.pnames);
end

prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data, y, J] = QuasiP_PC(prob, data, u)
   y = 0; J =1;
end



function [data, y, J] = QuasiP_PC(prob, data, u)
% Generator phase equations:
%   y_PC = [OPC_(e+1),...,OPC_d]^T = 0.
% Each row removes one phase-shift null direction.
%QUASIP_PC  Dispatch phase-condition evaluation to the selected formulation.
%
% Select the MP m-VCF or MT SHT phase-condition routine. The phase conditions
% remove the phase-shift nullspace associated with unknown response frequencies.
%
% On input:
%   prob : COCO problem structure.
%   data : QuasiP model and PHI data.
%   u    : Unknown vector [C(:);Om(:);p(:)].
%
% On output:
%   data : Updated formulation data.
%   y    : Phase-condition residual.
%   J    : Phase-condition Jacobian.
switch data.PHI.param
    case{'MP'}
        switch data.PHI.method
            case{'mVCF'}
                [data, y, J] = MP_mVCF_QuasiP_PC(prob, data, u);
            case{'sCCF'}
        end
    case{'MT'}
        switch data.PHI.method
            case{'DM'}

            case{'SHT'}
            switch data.PHI.type
                case 'HT'
                    [data, y, J] = MT_SHT_HT_QuasiP_PC(prob, data, u);
                case 'MF'
                    [data, y, J] = MT_SHT_MF_QuasiP_PC(prob, data, u);
            end                
        end
end
end






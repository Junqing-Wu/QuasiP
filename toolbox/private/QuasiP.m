function [data, y, J] = QuasiP(prob, data, u)
switch data.PHI.parametrization_tori
    case{'MP'}
        switch data.PHI.type
            case{'mVCF'}
                [data, y, J] = MP_mVCF_QuasiP(prob, data, u);
            case{'sCCF'}
        end
    case{'MT'}
        switch data.PHI.type
            case{'DM'}

            case{'SHT'}
                [data, y, J] = MT_SHT_QuasiP(prob, data, u);
        end
end
end
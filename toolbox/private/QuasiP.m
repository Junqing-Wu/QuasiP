function [data, y, J] = QuasiP(prob, data, u)
switch data.PHI.param
    case{'MP'}
        switch data.PHI.method
            case{'mVCF'}
                [data, y, J] = MP_mVCF_QuasiP(prob, data, u);
            case{'sCCF'}
        end
    case{'MT'}
        switch data.PHI.method
            case{'DM'}

            case{'SHT'}
                [data, y, J] = MT_SHT_QuasiP(prob, data, u);
        end
end
end
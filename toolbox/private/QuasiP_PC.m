function [data, y, J] = QuasiP_PC(prob, data, u)
switch data.PHI.parametrization_tori
    case{'MP'}
        switch data.PHI.type
            case{'mVCF'}
                [data, y, J] = MP_mVCF_QuasiP_PC(prob, data, u);
            case{'sCCF'}
        end
    case{'MT'}
        switch data.PHI.type
            case{'DM'}

            case{'SHT'}
            switch data.PHI.method
                case 'HT'
                    [data, y, J] = MT_SHT_HT_QuasiP_PC(prob, data, u);
                case 'MF'
                    [data, y, J] = MT_SHT_MF_QuasiP_PC(prob, data, u);
            end                
        end
end
end







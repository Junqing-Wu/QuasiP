function Amp = PostP_amp(data,Cii,Omii,outdof,isuamp)
switch data.PHI.parametrization_tori
    case{'MP'}
        switch data.PHI.type
            case{'mVCF'}
                Amp = MP_mVCF_PostP_amp(data,Cii,Omii,outdof,isuamp);
            case{'sCCF'}
        end
    case{'MT'}
        switch data.PHI.type
            case{'DM'}
            case{'SHOOT'}
        end
end
end
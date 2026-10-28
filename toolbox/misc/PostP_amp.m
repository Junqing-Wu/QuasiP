function Amp = PostP_amp(data,Cii,Omii,outdof,isuamp)
switch data.PHI.param
    case{'MP'}
        switch data.PHI.method
            case{'mVCF'}
                Amp = MP_mVCF_PostP_amp(data,Cii,Omii,outdof,isuamp);
            case{'sCCF'}
        end
    case{'MT'}
        switch data.PHI.method
            case{'DM'}
            case{'SHT'}
                Amp =  MT_SHT_PostP_amp(data,Omii,outdof,isuamp);
        end
end
end
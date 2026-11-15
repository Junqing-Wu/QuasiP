function Amp = PostP_amp(data,Cii,Omii,outdof,isuamp)
%POSTP_AMP  Dispatch amplitude processing to the selected Generator method.
%
% Select MP_mVCF_PostP_amp for the multi-period m-VCF representation or
% MT_SHT_PostP_amp for the multi-trajectory FSE-SHT representation.
%
% On input:
%   data    : QuasiP data containing PHI and amplitude settings.
%   Cii     : Current coefficient block.
%   Omii    : Current base-frequency vector.
%   outdof  : Requested state-component indices.
%   isuamp  : Displacement/velocity amplitude flag.
%
% On output:
%   Amp     : Formulation-specific amplitude data.

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
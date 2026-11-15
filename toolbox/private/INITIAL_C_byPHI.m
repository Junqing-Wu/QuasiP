function C = INITIAL_C_byPHI(data,Xk_MF,OM)
%INITIAL_C_BYPHI  Dispatch initial-coefficient construction by parameterization.
%
% The same user-level Xk_MF representation can initialize either the MP m-VCF
% formulation or the MT SHT formulation. This dispatcher selects the conversion
% routine associated with data.PHI.param and data.PHI.method.
%
% On input:
%   data  : QuasiP data containing the selected PHI parameterization.
%   Xk_MF : Initial multi-frequency response representation.
%   OM    : Initial base-frequency vector, used by shooting initialization.
%
% On output:
%   C     : Initial coefficient or shooting-data vector.
PHI = data.PHI;
switch PHI.param
  case{'MP'} 
      switch PHI.method
          case{'mVCF'}
               C = MP_mVCF_INITIAL_C_byPHI(data,Xk_MF);
          case{'sCCF'}
               C = MP_sCCF_INITIAL_C_byPHI(data,Xk_MF,OM);     
      end
  case{'MT'} 
      switch PHI.method
          case{'DM'}
               C = MT_DM_INITIAL_C_byPHI(data,Xk_MF,OM);
          case{'SHT'}
               C = MT_SHT_INITIAL_C_byPHI(data,Xk_MF,OM);
      end
end
end
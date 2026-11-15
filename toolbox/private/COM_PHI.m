function PHI = COM_PHI(PHI_INITIAL)
% Parameterization dispatch equation:
%   PHI = COM_PHI(PHI_INITIAL).
% The selected PHI.param/PHI.method determines the Generator representation.
%COM_PHI  Select the discretization constructor for the QuasiP parameterization.
%
% Dispatch to MP_mVCF_COM_PHI or MT_SHT_COM_PHI according to PHI_INITIAL.param
% and PHI_INITIAL.method. This keeps the public settings structure independent
% of the private data layout used by each numerical formulation.
%
% On input:
%   PHI_INITIAL : User-specified parameterization and discretization settings.
%
% On output:
%   PHI : Complete private discretization structure.
switch PHI_INITIAL.param
  case{'MP'} 
      switch PHI_INITIAL.method
          case{'mVCF'}
               PHI = MP_mVCF_COM_PHI(PHI_INITIAL);
          case{'sCCF'}
               PHI = MP_sCCF_COM_PHI(PHI_INITIAL);     
      end
  case{'MT'} 
      switch PHI_INITIAL.method
          case{'DM'}
               PHI = MT_DM_COM_PHI(PHI_INITIAL);
          case{'SHT'}
               PHI = MT_SHT_COM_PHI(PHI_INITIAL);     
      end
end
end
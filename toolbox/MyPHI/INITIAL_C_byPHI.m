function C = INITIAL_C_byPHI(data,X_MF,OM)
PHI = data.PHI;
switch PHI.parametrization_tori
  case{'MP'} 
      switch PHI.type
          case{'mVCF'}
               C = MP_mVCF_INITIAL_C_byPHI(data,X_MF,OM);
          case{'sCCF'}
               C = MP_sCCF_INITIAL_C_byPHI(data,X_MF,OM);     
      end
  case{'MT'} 
      switch PHI_INITIAL.type
          case{'DM'}
               C = MT_DM_INITIAL_C_byPHI(data,X_MF,OM);
          case{'SHOOT'}
               C = MT_SHOOT_INITIAL_C_byPHI(data,X_MF,OM);     
      end
end
end
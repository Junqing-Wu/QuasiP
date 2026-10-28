function C = INITIAL_C_byPHI(data,X_MF,OM)
PHI = data.PHI;
switch PHI.param
  case{'MP'} 
      switch PHI.method
          case{'mVCF'}
               C = MP_mVCF_INITIAL_C_byPHI(data,X_MF);
          case{'sCCF'}
               C = MP_sCCF_INITIAL_C_byPHI(data,X_MF,OM);     
      end
  case{'MT'} 
      switch PHI.method
          case{'DM'}
               C = MT_DM_INITIAL_C_byPHI(data,X_MF,OM);
          case{'SHT'}
               C = MT_SHT_INITIAL_C_byPHI(data,X_MF,OM);
      end
end
end
function PHI = COM_PHI(PHI_INITIAL)
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

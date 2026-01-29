function PHI = COM_PHI(PHI_INITIAL)

switch PHI_INITIAL.parametrization_tori
  case{'MP'} 
      switch PHI_INITIAL.type
          case{'mVCF'}
               PHI = MP_mVCF_COM_PHI(PHI_INITIAL);
          case{'sCCF'}
               PHI = MP_sCCF_COM_PHI(PHI_INITIAL);     
      end
  case{'MT'} 
      switch PHI_INITIAL.type
          case{'DM'}
               PHI = MT_DM_COM_PHI(PHI_INITIAL);
          case{'SHOOT'}
               PHI = MT_SHOOT_COM_PHI(PHI_INITIAL);     
      end
end
end
MODULE Saha
   USE GlobalDeclarations
   USE PhysicsDeclarations
   USE DataDeclarations
   IMPLICIT NONE
CONTAINS


  SUBROUTINE SahaSetTemp(q, form, T)  !Cannot use q(iE)
    REAL(KIND=qPREC), INTENT(INOUT), DIMENSION(:) :: q
    INTEGER :: form
    REAL(KIND=qPREC) :: T
    SELECT CASE (form)
    CASE (CONSERVATIVE_FORM)
       q(iE)=3d0/2d0*q(1)*T/TempScale
       q(iE)=q(iE)+half*SUM(q(m_low:m_high)**2)/q(1)
       IF (lMHD) q(iE)=q(iE)+half*SUM(q(iBx:iBz)**2)
    CASE (SOURCE_FORM)
       q(iE)=3d0/2d0*q(1)*T/TempScale
    CASE (PRIMITIVE_FORM)
       q(iE)=q(1)*T/TempScale
    CASE (ROE_FORM) 
       q(iE)=5d0/2d0*q(1)*T/TempScale
       IF (lMHD) q(iE)=q(iE) + sum(q(iBx:iBz)**2)
       q(iE)=q(iE)+half*q(1)*sum(q(m_low:m_high)**2)
    END SELECT  
  END SUBROUTINE SahaSetTemp


  SUBROUTINE SahaSetPress(q, form, P)  
    REAL(KIND=qPREC), INTENT(INOUT), DIMENSION(:) :: q
    INTEGER :: form
    REAL(KIND=qPREC) :: P
    SELECT CASE (form)
    CASE (CONSERVATIVE_FORM)
       q(iE)=3d0/2d0*P
       q(iE)=q(iE)+half*SUM(q(m_low:m_high)**2)/q(1)
       IF (lMHD) q(iE)=q(iE)+half*SUM(q(iBx:iBz)**2)
    CASE (SOURCE_FORM)
       q(iE)=3d0/2d0*P
    CASE (PRIMITIVE_FORM)
       q(iE)=P
    CASE (ROE_FORM) 
       q(iE)=5d0/2d0*P
       IF (lMHD) q(iE)=q(iE) + sum(q(iBx:iBz)**2)
       q(iE)=q(iE)+half*q(1)*sum(q(m_low:m_high)**2)
    END SELECT  
  END SUBROUTINE SahaSetPress


  PURE FUNCTION SahadEdT(q, form)
    REAL(KIND=qPREC), INTENT(IN), DIMENSION(:) :: q
    INTEGER, INTENT(IN) :: form
    REAL(KIND=qPREC) :: SahadEdT
    SELECT CASE (form)
    CASE (CONSERVATIVE_FORM)
       SahadEdT=3d0/2d0*q(1)/TempScale
    CASE (SOURCE_FORM)
       SahadEdT=3d0/2d0*q(1)/TempScale
    CASE (PRIMITIVE_FORM)
       SahadEdT=3d0/2d0*q(1)/TempScale
    CASE (ROE_FORM) 
       SahadEdT=3d0/2d0*q(1)/TempScale
    END SELECT
    
  END FUNCTION SahadEdT


  PURE FUNCTION SahaPress(q, form)
    REAL(KIND=qPREC), INTENT(IN), DIMENSION(:) :: q
    INTEGER, INTENT(IN) :: form
    REAL(KIND=qPREC) :: SahaPress
    SELECT CASE (form)
    CASE (CONSERVATIVE_FORM)
       SahaPress=q(iE)-half*SUM(q(m_low:m_high)**2)/q(1)
       IF (lMHD) SahaPress=SahaPress-half*SUM(q(iBx:iBz)**2)
       SahaPress=2d0/3d0*SahaPress
    CASE (SOURCE_FORM)
       SahaPress=2d0/3d0*q(iE)
    CASE (PRIMITIVE_FORM)
       SahaPress=q(iE)
    CASE (ROE_FORM) 
       SahaPress=q(1)*q(iE)-half*q(1)*sum(q(m_low:m_high)**2)
       IF (lMHD) SahaPress=SahaPress-sum(q(iBx:iBz)**2)
       SahaPress=2d0/5d0*SahaPress
    END SELECT
  END FUNCTION SahaPress


  PURE FUNCTION SahaInternalEnergy(q, form)
    REAL(KIND=qPREC), INTENT(IN), DIMENSION(:) :: q
    INTEGER, INTENT(IN) :: form
    REAL(KIND=qPREC) :: SahaInternalEnergy
    SELECT CASE (form)
    CASE (CONSERVATIVE_FORM)
       SahaInternalEnergy=q(iE)-half*SUM(q(m_low:m_high)**2)/q(1)
       IF (lMHD) SahaInternalEnergy=SahaInternalEnergy-half*SUM(q(iBx:iBz)**2)
    CASE (SOURCE_FORM)
       SahaInternalEnergy=q(iE)
    CASE (PRIMITIVE_FORM)
       SahaInternalEnergy=3d0/2d0*q(iE)
    CASE (ROE_FORM) 
       SahaInternalEnergy=q(1)*q(iE)-half*q(1)*sum(q(m_low:m_high)**2)
       IF (lMHD) SahaInternalEnergy=SahaInternalEnergy-sum(q(iBx:iBz)**2)
       SahaInternalEnergy=3d0/5d0*SahaInternalEnergy
    END SELECT    
  END FUNCTION SahaInternalEnergy
  

  PURE FUNCTION SahaTemperature(q, form)
    REAL(KIND=qPREC), INTENT(IN), DIMENSION(:) :: q
    INTEGER, INTENT(IN) :: form
    REAL(KIND=qPREC) :: SahaTemperature
    SahaTemperature=SahaPress(q, form)/q(1)*TempScale
  END FUNCTION SahaTemperature

  
  PURE FUNCTION SahaSoundSpeed(q, form)
    REAL(KIND=qPREC), INTENT(IN), DIMENSION(:) :: q
    INTEGER, INTENT(IN) :: form
    REAL(KIND=qPREC) :: SahaSoundSpeed
    SahaSoundSpeed=sqrt(5d0/3d0*SahaPress(q, form)/q(1))
  END FUNCTION SahaSoundSpeed

END MODULE Saha

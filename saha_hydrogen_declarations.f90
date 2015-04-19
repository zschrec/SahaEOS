!
!Declarations Needed For Hydrogen Saha
!

Module SahaHydrogenDeclarations
    IMPLICIT NONE
    SAVE
    PUBLIC
    
    !digits of accuracy
    INTEGER, PARAMETER :: ikind=selected_REAL_KIND(p=20)
    !number of elements to consider
    INTEGER, PARAMETER :: elements=1
    !temperature to iterate over
    INTEGER  :: T
    !*************************************
    !Physical Constants - Begin
    !*************************************
    !kg - mass of electron
    REAL(KIND=ikind), PARAMETER :: m_e = 9.10938291e-31
    ! m^2 kg s^-1 - planck's constant
    REAL (KIND=ikind), PARAMETER :: h = 6.62606957e-34
    ! cm^2 g s^-1 - planck's constant
    REAL (KIND=ikind), PARAMETER :: h_cgs = 6.62606957e-27
    !kg/m^3 - density
    REAL(KIND=ikind), PARAMETER :: rho = 1.0e-3
    !g/cm^3 - density
    REAL(KIND=ikind), PARAMETER :: rho_cgs = rho*1.0e-3
    !m^2 kg s^-2 K^-1 - boltzmann's constant
    REAL(KIND=ikind), PARAMETER :: k_kg = 1.3806488e-23
    !cm^2 g s^-2 K^-1 - boltzmann's constant
    REAL(KIND=ikind), PARAMETER :: k_cgs = 1.3806488e-16
    !eV K^-1 - boltzmann's constant
    REAL(KIND=ikind), PARAMETER :: k_eV = 8.6173324e-5 
    !Avogadro's number
    REAL(KIND=ikind), PARAMETER :: N_o = 6.0221413e23
    !pi
    REAL(KIND=ikind), PARAMETER :: pi = 4.0*atan(1.0)
    !constant in ratio calculation
    REAL (KIND=ikind), PARAMETER :: const = ((2*pi*m_e*k_kg)/(h**2))**(3.0/2.0)
    !*************************************
    !Physical Constants - End 
    !*************************************
 
    !arrays of constants for elements/ionization levels
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: x = (/ 1 /)
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: A = (/ 1.0080 /)
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: A_cgs = (/ 1.673e-24 /)
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: A_kg = (/ 1.673e-27 /)
    REAL (KIND=ikind),DIMENSION(elements) :: v_e
    !hydrogen ionization constants
    REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h = (/ 13.5984 /)
    REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h_cgs = (/ 2.1787e-11 /)
    REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h_kg = (/ 2.1787e-18 /)
    REAL (KIND=ikind), DIMENSION(2), PARAMETER :: B_h = (/ 10.0**0.3, 10.0**0.0 /)
    !iteration variables
    REAL (KIND=ikind) :: N_e, N_e1
    REAL (KIND=ikind), PARAMETER :: tol  = 1.0e-6
    INTEGER (KIND=ikind), PARAMETER :: lower_T  = 100
    INTEGER (KIND=ikind), PARAMETER :: upper_T  = 20000
    INTEGER (KIND=ikind), PARAMETER :: increment_T  = 100
    !record variables
    REAL (KIND=ikind) :: N_e_initial, N_e_last
    !ionization degree - hydrogen
    REAL (KIND=ikind), DIMENSION(1) :: y_h
    !ionization ratio - hydrogen
    REAL (KIND=ikind), DIMENSION(1) :: ratio_h

CONTAINS

      FUNCTION Pressure(ion_frac,T)
        REAL(KIND=ikind) :: Pressure, ion_frac 
        INTEGER  :: T
        Pressure=(1+ion_frac)*(k_cgs/A_cgs(1))*rho_cgs*T
      END FUNCTION Pressure

      FUNCTION InternalEnergy(ion_frac,T)
        REAL(KIND=ikind) :: InternalEnergy, ion_frac 
        INTEGER  :: T
        InternalEnergy=1.5*Pressure(ion_frac,T)+ion_frac*(X_h_cgs(1)/A_cgs(1))*rho_cgs
      END FUNCTION InternalEnergy

      FUNCTION SpecificHeatConstantVolume(ion_frac,T)
        REAL(KIND=ikind) :: SpecificHeatConstantVolume, ion_frac
        INTEGER  :: T
        SpecificHeatConstantVolume=(k_cgs/A_cgs(1))*(1.5*(1.0+ion_frac) + &
              ((1.5+ (X_h_cgs(1)/(k_cgs*T)))**2.0)*((ion_frac*(1.0-ion_frac))/(2.0-ion_frac)))
      END FUNCTION SpecificHeatConstantVolume


END Module SahaHydrogenDeclarations



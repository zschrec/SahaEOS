
!
! Declarations file for stand alone solver modules
!

MODULE SahaHydrogenInternalEnergyEOSDeclarations
  IMPLICIT NONE
  SAVE
  PUBLIC
  
  !digits of accuracy
  INTEGER, PARAMETER :: ikind=SELECTED_REAL_KIND(p=20)
  !number of elements to consider
  INTEGER, PARAMETER :: elements=1
  ! K - temperature
  REAL (KIND=ikind) :: T
  
  INTEGER :: here_count = 1
  INTEGER :: there_count = 1
  INTEGER :: ham_count = 0 
  INTEGER :: banana_count = 0
  INTEGER :: iter_count = 0

  !*************************************
  !*************************************
  !Physical Constants - Begin
  !*************************************
  !g - mass of electron
  REAL(KIND=ikind), PARAMETER :: m_e = 9.10938291e-28
  ! cm^2 g s^-1 - planck's constant
  REAL (KIND=ikind), PARAMETER :: h = 6.62606957e-27
  !g/cm^3 - density
  REAL(KIND=ikind), PARAMETER :: rho = 1.0e-9
  !cm^2 g s^-2 K^-1 - boltzmann's constant
  REAL(KIND=ikind), PARAMETER :: k = 1.3806488e-16
  !eV K^-1 - boltzmann's constant
  REAL(KIND=ikind), PARAMETER :: k_eV = 8.6173324e-5 
  !Avogadro's number
  REAL(KIND=ikind), PARAMETER :: N_o = 6.0221413e23
  !pi
  REAL(KIND=ikind), PARAMETER :: pi = 4.0*atan(1.0)
  !constant in ratio calculation
  REAL (KIND=ikind), PARAMETER :: const = ((2*pi*m_e*k)/(h**2))**(3.0/2.0)
  !*************************************
  !Physical Constants - End 
  !*************************************

  !arrays of constants for elements/ionization levels
  ! g of element per g of mixture
  REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: x = (/ 1.0 /)
  ! g/mol
  REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: A = (/ 1.0080 /)
  ! g - mass of atoms
  REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: m = (/ 1.673e-24 /)
  ! number of atom of an element
  REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: n = N_o*rho*(/ x(1)/A(1) /)
  ! average electron contribution
  REAL (KIND=ikind),DIMENSION(elements) :: v_e
  
  !hydrogen ionization constants
  ! eV - ionization energies
  REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h_eV = (/ 13.5984 /)
  REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h = (/ 2.1787e-11 /)
  ! partition functions
  REAL (KIND=ikind), DIMENSION(2), PARAMETER :: B_h = (/ 10.0**0.3, 10.0**0.0 /)
  !electron density
  REAL (KIND=ikind) :: N_e = 1.0e3 
!  REAL (KIND=ikind) :: N_e = 1.0 
  REAL (KIND=ikind) :: N_e1 = 0
  REAL (KIND=ikind) :: g_N_e
  !iteration variables
  INTEGER :: U
 REAL, PARAMETER :: tol  = 1.0e-9
!  REAL, PARAMETER :: tol  = 1.0e-6
  ! CGS 
  INTEGER, PARAMETER :: lower_U  = 100
  INTEGER, PARAMETER :: upper_U  = 18000
  INTEGER, PARAMETER :: increment_U = 100


  !record variables
  REAL (KIND=ikind) :: N_e_initial, N_e_last, N_e1_last
  REAL (KIND=ikind) :: T_last = 0.0
  REAL (KIND=ikind) :: T_last_last = 0.0
  REAL (KIND=ikind) :: N_e_last_temp, N_e1_last_temp
  !ionization degree array - hydrogen
  REAL (KIND=ikind), DIMENSION(1) :: y_h
  !ionization ratio - hydrogen
  REAL (KIND=ikind), DIMENSION(1) :: ratio_h

CONTAINS

      FUNCTION Pressure(ion_frac,T)
        REAL (KIND=ikind) :: Pressure, ion_frac 
        REAL (KIND=ikind) :: T
        Pressure=(1+ion_frac)*(k/m(1))*rho*T
      END FUNCTION Pressure

      FUNCTION InternalEnergy(ion_frac,T)
        REAL(KIND=ikind) :: InternalEnergy, ion_frac
        REAL (KIND=ikind) :: T
        InternalEnergy=1.5*Pressure(ion_frac,T)+ion_frac*(X_h(1)/m(1))*rho
      END FUNCTION InternalEnergy

      FUNCTION SpecificHeatConstantVolume(ion_frac,T)
        REAL(KIND=ikind) :: SpecificHeatConstantVolume, ion_frac
        REAL (KIND=ikind) :: T
        SpecificHeatConstantVolume=(k/m(1))*(1.5*(1.0+ion_frac) + &
              ((1.5+ (X_h(1)/(k*T)))**2.0)*((ion_frac*(1.0-ion_frac))/(2.0-ion_frac)))
      END FUNCTION SpecificHeatConstantVolume


END MODULE SahaHydrogenInternalEnergyEOSDeclarations


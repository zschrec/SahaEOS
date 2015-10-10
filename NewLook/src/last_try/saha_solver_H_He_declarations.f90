
!
! Declarations file for stand alone solver modules
!

MODULE SahaSolverHHeDeclarations
  IMPLICIT NONE
  SAVE
  PUBLIC
  
  !digits of accuracy
  INTEGER, PARAMETER :: ikind=SELECTED_REAL_KIND(p=20)
  !number of elements to consider
  INTEGER, PARAMETER :: elements=2
  !temperature to iterate over
  INTEGER :: T

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
  ! just hydrogen
  !REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: x = (/ 1.0, 0.0 /)
  ! only hydrogen and helium
  REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: x = (/ 0.7, 0.3 /)
  ! g/mol
  REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: A = (/ 1.0080, 4.002602 /)
  ! average electron contribution
  REAL (KIND=ikind),DIMENSION(elements) :: v_e
  
  !hydrogen ionization constants
  ! eV - ionization energies
  REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h_eV = (/ 13.5984 /)
  REAL (KIND=ikind), DIMENSION(2), PARAMETER :: X_he_eV = (/ 24.587, 54.416 /)
  ! partition functions
  REAL (KIND=ikind), DIMENSION(2), PARAMETER :: B_h = (/ 10.0**0.3, 10.0**0.0 /)
  REAL (KIND=ikind), DIMENSION(3), PARAMETER :: B_he = (/ 10.0**0.0, 10.0**0.3, 1.0/)
  !iteration variables
  REAL (KIND=ikind) :: N_e = 1e20
  REAL (KIND=ikind) :: N_e1 = 0
  REAL (KIND=ikind), PARAMETER :: tol  = 1.0e-6
  ! Kelvin
  INTEGER (KIND=ikind), PARAMETER :: lower_T  = 100
  INTEGER (KIND=ikind), PARAMETER :: upper_T  = 60000
  !INTEGER (KIND=ikind), PARAMETER :: upper_T  = 20000
  INTEGER (KIND=ikind), PARAMETER :: increment_T  = 100

  !record variables
  REAL (KIND=ikind) :: N_e_initial, N_e_last
  !ionization degree array - hydrogen
  REAL (KIND=ikind), DIMENSION(1) :: y_h
  !ionization ratio - hydrogen
  REAL (KIND=ikind), DIMENSION(1) :: ratio_h
  !ionization degree array - helium
  REAL (KIND=ikind), DIMENSION(2) :: y_he
  !ionization ratio - helium
  REAL (KIND=ikind), DIMENSION(2) :: ratio_he

END MODULE SahaSolverHHeDeclarations


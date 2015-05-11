!
!Declarations Needed For Hydrogen Saha
!

Module SahaHydrogenEOSAstroBEARDeclarations
  IMPLICIT NONE
  SAVE
  PUBLIC
    
  !digits of accuracy
  INTEGER, PARAMETER :: ikind=selected_REAL_KIND(p=20)
  !number of elements to consider
  INTEGER, PARAMETER :: elements=1
  !temperature to iterate over
  INTEGER (KIND=selected_INT_KIND(20))  :: T
  ! temperature to be used in iterations
  REAL (KIND=ikind) :: T_real
  ! temperature index
  INTEGER  :: T_index
  ! density to be used in iterations
  REAL :: rho_iter 
  ! density index
  INTEGER :: rho_index
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
  REAL(KIND=ikind) :: rho_enter
  REAL(KIND=ikind) :: internal_energy_enter
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
  ! average electron contribution
  REAL (KIND=ikind),DIMENSION(elements) :: v_e

  !hydrogen ionization constants
  ! eV - ionization energies
  REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h_eV = (/ 13.5984 /)
  REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h = (/ 2.1787e-11 /)
  ! partition functions
  REAL (KIND=ikind), DIMENSION(2), PARAMETER :: B_h = (/ 10.0**0.3, 10.0**0.0 /)
  
  !iteration variables
  REAL (KIND=ikind) :: N_e = 1e20
  REAL (KIND=ikind) :: N_e1 = 0
  REAL (KIND=ikind), PARAMETER :: tol  = 1.0e-6
  ! Kelvin
  INTEGER (KIND=ikind), PARAMETER :: lower_T  = 100
  INTEGER (KIND=ikind), PARAMETER :: upper_T  = 20000
  INTEGER (KIND=ikind), PARAMETER :: increment_T  = 100
  !number of temperatures
  INTEGER, PARAMETER :: num_temps = (upper_T - lower_T)/increment_T + 1
  ! g/cm^3  
  INTEGER, PARAMETER :: lower_rho = -10 ! 1.0e-10 
  INTEGER, PARAMETER :: upper_rho = -8 ! 1.0e-8
  REAL, PARAMETER :: increment_rho = 0.05 ! 10^0.5
  ! number of densities
  INTEGER, PARAMETER :: num_rho = (upper_rho - lower_rho)/increment_rho + 1
  !array of internal energies
  REAL(KIND=ikind), DIMENSION(num_temps, num_rho) :: internal_energies 
  
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
    Pressure=(1+ion_frac)*(k/m(1))*rho*T
  END FUNCTION Pressure

  FUNCTION InternalEnergy(ion_frac,T)
    REAL(KIND=ikind) :: InternalEnergy, ion_frac
    INTEGER  :: T
    InternalEnergy=1.5*Pressure(ion_frac,T)+ion_frac*(X_h(1)/m(1))*rho
  END FUNCTION InternalEnergy

  FUNCTION Pressure_RHO(ion_frac,T,current_rho)
    REAL(KIND=ikind) :: Pressure_RHO, ion_frac,current_rho, T
    Pressure_RHO=(1+ion_frac)*(k/m(1))*current_rho*T
  END FUNCTION Pressure_RHO

  FUNCTION InternalEnergy_RHO(ion_frac,T,current_rho)
    REAL(KIND=ikind) :: InternalEnergy_RHO, ion_frac, current_rho, T
    InternalEnergy_RHO=1.5*Pressure_RHO(ion_frac,T,current_rho)+ &
          ion_frac*(X_h(1)/m(1))*current_rho
  END FUNCTION InternalEnergy_RHO

  FUNCTION SpecificHeatConstantVolume(ion_frac,T)
    REAL(KIND=ikind) :: SpecificHeatConstantVolume, ion_frac, T
    SpecificHeatConstantVolume=(k/m(1))*(1.5*(1.0+ion_frac) + &
          ((1.5+ (X_h(1)/(k*T)))**2.0)*((ion_frac*(1.0-ion_frac))/(2.0-ion_frac)))
  END FUNCTION SpecificHeatConstantVolume

  FUNCTION GetTemperatureInternalEnergy(internal_energies, rho, internal_energy)
    REAL(KIND=ikind) :: GetTemperatureInternalEnergy, internal_energy, rho, m
    REAL(KIND=ikind), DIMENSION(num_temps, num_rho) :: internal_energies 
    REAL(KIND=ikind) :: approx_rho_1, approx_T_1, approx_Ug_1, offset_1
    REAL(KIND=ikind) :: approx_rho_2, approx_T_2, approx_Ug_2, offset_2
    
    rho_index = 1

    rho_iter = 1.0*lower_rho
    DO WHILE (rho_iter <= upper_rho)
                
      IF ( 10.0**rho_iter > rho) THEN
                
        T_index=1
        DO T = lower_T, upper_T, increment_T
              
          IF (internal_energies(T_index,rho_index) > internal_energy) THEN
            
            approx_rho_1 = (-1.0*increment_T)*&
                (LOG10(internal_energies(T_index,rho_index))- &
                  LOG10(internal_energies(T_index,rho_index-1)))
                
            approx_rho_2 = (-1.0*increment_T)*&
                (LOG10(internal_energies(T_index-1,rho_index))- &
                  LOG10(internal_energies(T_index-1,rho_index-1)))

            approx_T_1 = (1.0*increment_rho)*&
                (LOG10(internal_energies(T_index,rho_index-1))- &
                  LOG10(internal_energies(T_index,rho_index)))

            approx_T_2 = (1.0*increment_rho)*&
                (LOG10(internal_energies(T_index-1,rho_index-1))- &
                  LOG10(internal_energies(T_index,rho_index-1)))

            approx_Ug_1 = (increment_T*1.0)*(1.0*increment_rho)
            approx_Ug_2 = (increment_T*1.0)*(1.0*increment_rho)

            offset_1 = -1.0* &
              (approx_rho_1*(rho_iter-increment_rho)+approx_T_1* &
                (T-increment_T)+approx_Ug_1* &
                LOG10(internal_energies(T_index-1,rho_index-1)))

            offset_2 = -1.0* &
              (approx_rho_2*rho_iter+approx_T_2*T+approx_Ug_2* &
                LOG10(internal_energies(T_index,rho_index)))

            GetTemperatureInternalEnergy=((-1.0*(offset_1+approx_rho_1*LOG10(rho)+approx_Ug_1*&
              LOG10(internal_energy))/approx_T_1) + &
              (-1.0*(offset_2+approx_rho_2*LOG10(rho)+approx_Ug_2*&
              LOG10(internal_energy))/approx_T_2))/2.0

            RETURN 
              
          END IF 

          T_index = T_index+1
        END DO
      END IF
          
      rho_index = rho_index+1
      rho_iter = rho_iter+increment_rho
    END DO

  END FUNCTION GetTemperatureInternalEnergy


END Module SahaHydrogenEOSAstroBEARDeclarations



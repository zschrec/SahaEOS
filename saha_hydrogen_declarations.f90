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
    INTEGER (KIND=selected_INT_KIND(20))  :: T
    REAL (KIND=ikind) :: T_real
    INTEGER  :: T_iter
    !*************************************
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
    REAL(KIND=ikind) :: rho_enter
    REAL(KIND=ikind) :: internal_energy_enter
    INTEGER :: rho_iter 
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
    REAL (KIND=ikind) :: N_e = 1e20
    REAL (KIND=ikind) :: N_e1 = 0
    REAL (KIND=ikind), PARAMETER :: tol  = 1.0e-6
    INTEGER (KIND=ikind), PARAMETER :: lower_T  = 100
    INTEGER (KIND=ikind), PARAMETER :: upper_T  = 20000
    INTEGER (KIND=ikind), PARAMETER :: increment_T  = 100
    !number of densities
    INTEGER, PARAMETER :: num_rho = 5
    !array of densities
    REAL(KIND=ikind),DIMENSION(num_rho),PARAMETER :: rhos=(/ 1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1 /) 
    !number of temperatures
    INTEGER, PARAMETER :: num_temps = (upper_T - lower_T)/increment_T + 1
    !array of temperatures
    INTEGER, DIMENSION(num_temps) :: temps 
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
        Pressure=(1+ion_frac)*(k_cgs/A_cgs(1))*rho_cgs*T
      END FUNCTION Pressure

      FUNCTION InternalEnergy(ion_frac,T)
        REAL(KIND=ikind) :: InternalEnergy, ion_frac
        INTEGER  :: T
        InternalEnergy=1.5*Pressure(ion_frac,T)+ion_frac*(X_h_cgs(1)/A_cgs(1))*rho_cgs
      END FUNCTION InternalEnergy

      FUNCTION Pressure_RHO(ion_frac,T,current_rho_cgs)
        REAL(KIND=ikind) :: Pressure_RHO, ion_frac,current_rho_cgs, T
        !INTEGER  :: T
        Pressure_RHO=(1+ion_frac)*(k_cgs/A_cgs(1))*current_rho_cgs*T
      END FUNCTION Pressure_RHO

      FUNCTION InternalEnergy_RHO(ion_frac,T,current_rho_cgs)
        REAL(KIND=ikind) :: InternalEnergy_RHO, ion_frac, current_rho_cgs, T
        !INTEGER  :: T
        InternalEnergy_RHO=1.5*Pressure_RHO(ion_frac,T,current_rho_cgs)+ &
              ion_frac*(X_h_cgs(1)/A_cgs(1))*current_rho_cgs
      END FUNCTION InternalEnergy_RHO

      FUNCTION SpecificHeatConstantVolume(ion_frac,T)
        REAL(KIND=ikind) :: SpecificHeatConstantVolume, ion_frac, T
        !INTEGER  :: T
        SpecificHeatConstantVolume=(k_cgs/A_cgs(1))*(1.5*(1.0+ion_frac) + &
              ((1.5+ (X_h_cgs(1)/(k_cgs*T)))**2.0)*((ion_frac*(1.0-ion_frac))/(2.0-ion_frac)))
      END FUNCTION SpecificHeatConstantVolume

      FUNCTION GetTemperatureInternalEnergy(internal_energies, rho, internal_energy)
        REAL(KIND=ikind) :: GetTemperatureInternalEnergy, internal_energy, rho, m
        REAL(KIND=ikind), DIMENSION(num_temps, num_rho) :: internal_energies 
        
        DO rho_iter = 1, num_rho


          IF ( ABS(rho - rhos(rho_iter)*1.0e-3)/rho < tol) THEN

            T_iter=1
            DO T = lower_T, upper_T, increment_T
              
              IF (internal_energies(T_iter,rho_iter) > internal_energy) THEN
                m = (increment_T*1.0)/ &
                    (internal_energies(T_iter,rho_iter)-internal_energies(T_iter-1,rho_iter))
               
                !PRINT *,internal_energies(T_iter,rho_iter)

                GetTemperatureInternalEnergy= &
                    m*(internal_energy-internal_energies(T_iter,rho_iter))+T

                RETURN 
              
              END IF 

              T_iter = T_iter+1
            END DO
          END IF

        END DO

      END FUNCTION GetTemperatureInternalEnergy


END Module SahaHydrogenDeclarations



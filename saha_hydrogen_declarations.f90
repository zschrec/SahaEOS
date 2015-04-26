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
    INTEGER :: rho_index
    !*************************************
    !*************************************
    !Physical Constants - Begin
    !*************************************
    !kg - mass of electron
    REAL(KIND=ikind), PARAMETER :: m_e_mks = 9.10938291e-31
    !g - mass of electron
    REAL(KIND=ikind), PARAMETER :: m_e = 9.10938291e-28
    ! m^2 kg s^-1 - planck's constant
    REAL (KIND=ikind), PARAMETER :: h_mks = 6.62606957e-34
    ! cm^2 g s^-1 - planck's constant
    REAL (KIND=ikind), PARAMETER :: h = 6.62606957e-27
    !g/cm^3 - density
    REAL(KIND=ikind), PARAMETER :: rho = 1.0e-3
!    !kg/m^3 - density
!    REAL(KIND=ikind), PARAMETER :: rho_mks = rho*1.0e3
    REAL(KIND=ikind) :: rho_enter
    REAL(KIND=ikind) :: internal_energy_enter
    INTEGER :: rho_iter 
    !m^2 kg s^-2 K^-1 - boltzmann's constant
    REAL(KIND=ikind), PARAMETER :: k_mks = 1.3806488e-23
    !cm^2 g s^-2 K^-1 - boltzmann's constant
    REAL(KIND=ikind), PARAMETER :: k = 1.3806488e-16
    !eV K^-1 - boltzmann's constant
    REAL(KIND=ikind), PARAMETER :: k_eV = 8.6173324e-5 
    !Avogadro's number
    REAL(KIND=ikind), PARAMETER :: N_o = 6.0221413e23
    !pi
    REAL(KIND=ikind), PARAMETER :: pi = 4.0*atan(1.0)
    !constant in ratio calculation
    REAL (KIND=ikind), PARAMETER :: const_mks = ((2*pi*m_e_mks*k_mks)/(h_mks**2))**(3.0/2.0)
    REAL (KIND=ikind), PARAMETER :: const = ((2*pi*m_e*k)/(h**2))**(3.0/2.0)
    !*************************************
    !Physical Constants - End 
    !*************************************
 
    !arrays of constants for elements/ionization levels
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: x = (/ 1 /)
    ! g/mol
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: A = (/ 1.0080 /)
    ! kg/mol
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: A_mks = (/ 1.0080e-3 /)
    ! g
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: m_h = (/ 1.673e-24 /)
    ! kg
    REAL (KIND=ikind),DIMENSION(elements), PARAMETER :: m_h_mks = (/ 1.673e-27 /)
    REAL (KIND=ikind),DIMENSION(elements) :: v_e
    !hydrogen ionization constants
    REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h_eV = (/ 13.5984 /)
    ! g
    REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h = (/ 2.1787e-11 /)
    ! kg
    REAL (KIND=ikind), DIMENSION(1), PARAMETER :: X_h_kg = (/ 2.1787e-18 /)
    REAL (KIND=ikind), DIMENSION(2), PARAMETER :: B_h = (/ 10.0**0.3, 10.0**0.0 /)
    !iteration variables
    REAL (KIND=ikind) :: N_e = 1e20
    REAL (KIND=ikind) :: N_e1 = 0
    REAL (KIND=ikind), PARAMETER :: tol  = 1.0e-6
    INTEGER (KIND=ikind), PARAMETER :: lower_T  = 100
    INTEGER (KIND=ikind), PARAMETER :: upper_T  = 20000
    INTEGER (KIND=ikind), PARAMETER :: increment_T  = 100
    !!number of densities
    !INTEGER, PARAMETER :: num_rho = 5
    !!array of densities
    !REAL(KIND=ikind),DIMENSION(num_rho),PARAMETER::rhos=(/ 1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-4 /)
    ! orders of magnitude
    INTEGER(KIND=ikind),PARAMETER :: lower_rho = -12 ! 1.0e-8 
    INTEGER(KIND=ikind),PARAMETER :: upper_rho = -6 ! 1.0e-1
    INTEGER(KIND=ikind),PARAMETER :: increment_rho = 1 ! rho = rho*10^1
    INTEGER, PARAMETER :: num_rho = (upper_rho - lower_rho)/increment_rho + 1
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
        Pressure=(1+ion_frac)*(k/m_h(1))*rho*T
      END FUNCTION Pressure

      FUNCTION InternalEnergy(ion_frac,T)
        REAL(KIND=ikind) :: InternalEnergy, ion_frac
        INTEGER  :: T
        InternalEnergy=1.5*Pressure(ion_frac,T)+ion_frac*(X_h(1)/m_h(1))*rho
      END FUNCTION InternalEnergy

      FUNCTION Pressure_RHO(ion_frac,T,current_rho)
        REAL(KIND=ikind) :: Pressure_RHO, ion_frac,current_rho, T
        !INTEGER  :: T
        Pressure_RHO=(1+ion_frac)*(k/m_h(1))*current_rho*T
      END FUNCTION Pressure_RHO

      FUNCTION InternalEnergy_RHO(ion_frac,T,current_rho)
        REAL(KIND=ikind) :: InternalEnergy_RHO, ion_frac, current_rho, T
        !INTEGER  :: T
        InternalEnergy_RHO=1.5*Pressure_RHO(ion_frac,T,current_rho)+ &
              ion_frac*(X_h(1)/m_h(1))*current_rho
      END FUNCTION InternalEnergy_RHO

      FUNCTION SpecificHeatConstantVolume(ion_frac,T)
        REAL(KIND=ikind) :: SpecificHeatConstantVolume, ion_frac, T
        SpecificHeatConstantVolume=(k/m_h(1))*(1.5*(1.0+ion_frac) + &
              ((1.5+ (X_h(1)/(k*T)))**2.0)*((ion_frac*(1.0-ion_frac))/(2.0-ion_frac)))
      END FUNCTION SpecificHeatConstantVolume

      FUNCTION GetTemperatureInternalEnergy(internal_energies, rho, internal_energy)
        REAL(KIND=ikind) :: GetTemperatureInternalEnergy, internal_energy, rho, m
        REAL(KIND=ikind), DIMENSION(num_temps, num_rho) :: internal_energies 
        REAL(KIND=ikind) :: approx_rho, approx_T, approx_Ug, offset
        
        rho_index = 1

        DO rho_iter = lower_rho, upper_rho, increment_rho


          IF ( 10**rho_iter > rho) THEN

            T_iter=1
            DO T = lower_T, upper_T, increment_T
              
              IF (internal_energies(T_iter,rho_index) > internal_energy) THEN

                PRINT *,'RHO', rho, 10.0**rho_iter, rho_index
                PRINT *,'T', T

                approx_rho = (-1.0*increment_T)*&
                    (internal_energies(T_iter,rho_index)-internal_energies(T_iter,rho_index-1))

                !approx_T = (rhos(rho_iter) - rhos(rho_iter-1))*&
                !    (internal_energies(T_iter,rho_iter-1)-internal_energies(T_iter,rho_iter))

                approx_T = (1.0*increment_rho)*&
                    (internal_energies(T_iter,rho_index-1)-internal_energies(T_iter,rho_index))

                !approx_Ug = (increment_T*1.0)*(rhos(rho_iter) - rhos(rho_iter-1))

                approx_Ug = (increment_T*1.0)*(1.0*increment_rho)

                !offset = -1.0* &
                !    (approx_rho*rhos(rho_iter)+approx_T*T+approx_Ug*internal_energies(T_iter,rho_iter))

                offset = -1.0* &
                  (approx_rho*rho_iter+approx_T*T+approx_Ug*internal_energies(T_iter,rho_index))

                !GetTemperatureInternalEnergy= -1.0*&
                !    (offset+approx_rho*rho+approx_Ug*internal_energy)/approx_T

                GetTemperatureInternalEnergy= -1.0*&
                    (offset+approx_rho*LOG(rho)+approx_Ug*internal_energy)/approx_T

                RETURN 
              
              END IF 

              T_iter = T_iter+1
            END DO
          END IF
          
          rho_index = rho_index+1

        END DO

      END FUNCTION GetTemperatureInternalEnergy


END Module SahaHydrogenDeclarations



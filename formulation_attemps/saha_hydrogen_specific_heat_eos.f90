

PROGRAM saha_hydrogen_specific_heat_eos
  USE SahaHydrogenSpecificHeatEOSDeclarations
  IMPLICIT NONE
  
  T_initial = T
  T_last = T

  OPEN(10, FILE='outfiles/saha_hydrogen_internal_energy_eos.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure', ' ', 'Internal_Energy_Initial', ' ', &
    'Internal_Energy_Calc', ' ', 'Specific_Heat'

  ! CGS
  DO U = lower_U, upper_U, increment_U
    
    T = T_initial
    T_last = T_initial
    T_1 = 0.0

    PRINT *, U

    DO

      N_e = (1.0*U-1.5*(rho/m(1))*k*T) / (1.5*k*T + X_h(1)) 

      ! Calculate the Saha Equation for hydrogen 
      ratio_h(1) = (const * (2.0 * B_h(2) / B_h(1)) * (T)**(3.0/2.0) * & 
        EXP(-X_h_eV(1)/(k_eV*T)) ) / ( N_e )
      
      ! Calculate the degree of ionization (aka ionization fraction)
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

      ! Calculate average electron contribution 
      v_e(1) = 1.0*y_h(1)
            
      ! Compute electron density
      N_e =  (rho*N_o) * (x(1)/A(1)) * v_e(1) 
     
      ! Recompute temperature
      T_1 = (1.0*U-N_e*X_h(1)) / (1.5 * k * ( (rho/m(1) ) + N_e ) )
      
      IF (ABS((T - T_1)/T) < tol) EXIT
     
      PRINT *, T, ' ', U

!      IF(ABS((T_last - T_1)/T_last) > 1e3) THEN
!        PRINT *, 'A'
!        T_last = T
!        T = T_1
!      ELSE IF(ABS((T_last - T_1)/T_last) > 1e-2) THEN
!        PRINT *, 'B', ' ', ABS((T_last - T_1)/T_last)
!        T_last = T
!        T = (T + T_1) / 2.0
!      ELSE
!        PRINT *, 'C'
!        T_temp = T
!        T = (T + T_1 + T_last) / 3.0
!        T_last = T_temp 
!      END IF

      T_error = T_1 - T
      T = T + gain*T_error

!      IF(ABS((T_last - T_1)/T_last) > 1e3) THEN
!        PRINT *, 'A'
!        T_last = T
!        T = T_1
!      ELSE 
!        PRINT *, 'B'
!        T_last = T
!        T = (T + T_1) / 2.0
!      END IF

    END DO  

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1), ' ', Pressure(y_h(1), T), ' ', U, ' ', &
      InternalEnergy(y_h(1), T), ' ', SpecificHeatConstantVolume(y_h(1), T)
          
  END DO

END PROGRAM saha_hydrogen_specific_heat_eos


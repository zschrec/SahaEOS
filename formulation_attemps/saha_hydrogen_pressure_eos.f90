

PROGRAM saha_hydrogen_pressure_eos
  USE SahaHydrogenPressureEOSDeclarations
  IMPLICIT NONE
  
  T_initial = T
  T_last = T

  OPEN(10, FILE='outfiles/saha_hydrogen_pressure_eos.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure_Initial', ' ', 'Pressure_Calc', ' ', &
    'Internal_Energy', ' ', 'Specific_Heat'

  ! CGS
  DO P = lower_P, upper_P, increment_P
    
    T = T_initial
    T_last = T_initial
    T_1 = 0.0

    DO

      ! Calculate the Saha Equation for hydrogen 
      ratio_h(1) = (const * (2.0 * B_h(2) / B_h(1)) * (T)**(3.0/2.0) * & 
        EXP(-X_h_eV(1)/(k_eV*T)) ) / ( ((1.0*P)/(k*T)) - (rho/m(1)) )
      
      ! Calculate the degree of ionization (aka ionization fraction)
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

      ! Calculate average electron contribution 
      v_e(1) = 1.0*y_h(1)
            
      ! Compute electron density
      N_e =  (rho*N_o) * (x(1)/A(1)) * v_e(1) 
     
      ! Recompute temperature
      T_1 = (1.0*P) / (k * ( (rho/m(1) ) + N_e ) )
      
      IF (ABS((T - T_1)/T) < tol) EXIT
      
      IF(ABS((T_last - T_1)/T_last) > 1e3) THEN
        T_last = T
        T = T_1
      ELSE IF(ABS((T_last - T_1)/T_last) > 1e0) THEN
        T_last = T
        T = (T + T_1) / 2.0
      ELSE
        T = (T + T_1 + T_last) / 3.0
        T_last = T
      END IF
 
    END DO  

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1), ' ', P, ' ', Pressure(y_h(1), T), ' ', &
      InternalEnergy(y_h(1), T), ' ', SpecificHeatConstantVolume(y_h(1), T)
          
  END DO

END PROGRAM saha_hydrogen_pressure_eos




PROGRAM saha_hydrogen_eos
  USE SahaHydrogenEOSDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e
 
  OPEN(10, FILE='outfiles/saha_hydrogen_eos.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure', ' ', 'Internal_Energy', ' ', 'Specific_Heat'

  !K
  DO T = lower_T, upper_T, increment_T
    
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0
          
    DO

      ! Calculate the Saha Equation for hydrogen 
      ratio_h(1) = const * (T)**(3.0/2.0) * & 
        (2.0 * B_h(2) * EXP(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e)
      
      ! Calculate the degree of ionization (aka ionization fraction)
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

      ! Calculate average electron contribution 
      v_e(1) = 1.0*y_h(1)
            
      ! Recompute electron density
      N_e1 =  (rho*N_o) * (x(1)/A(1)) * v_e(1) 
     
      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
        
      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN 
        N_e_last = N_e
        N_e = N_e1
      ELSE 
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF
        
    END DO  

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1), ' ', Pressure(y_h(1), T), ' ', &
      InternalEnergy(y_h(1), T), ' ', SpecificHeatConstantVolume(y_h(1), 1.0*T), ' ', N_e
    
  END DO

END PROGRAM saha_hydrogen_eos




PROGRAM saha_hydrogen_pressure_eos
  USE SahaHydrogenPressureEOSDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_hydrogen_pressure_eos.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure_Initial', ' ', 'Pressure_Calc'

  ! CGS
  DO P = lower_P, upper_P, increment_P
  
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0.0

    DO

      T = (1.0*P) / (k * (n(1) + N_e))

      PRINT *, 'P ', P
      PRINT *, 'N_e ', N_e
      PRINT *, 'T ', T
            
      ! Calculate the Saha Equation for hydrogen 
      ratio_h(1) = const * (2.0 * B_h(2) / B_h(1)) * &
          (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)) * (1.0/N_e)

      ! Calculate the degree of ionization (aka ionization fraction)
      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

      ! Calculate average electron contribution 
      v_e(1) = 1.0*y_h(1)
            
      ! Recompute electron density
      N_e1 = n(1) * v_e(1) 
     
      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
      
      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN
        PRINT *, 'HERE ', (N_e_last - N_e1)/N_e_last, ' ', (N_e - N_e1)/N_e
        N_e_last = N_e
        N_e = N_e1
      ELSE
        PRINT *, 'THERE ', (N_e_last - N_e1)/N_e_last, ' ', (N_e - N_e1)/N_e
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF
 
    END DO  

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1), ' ', P, ' ', Pressure(y_h(1), T)
          
  END DO

END PROGRAM saha_hydrogen_pressure_eos


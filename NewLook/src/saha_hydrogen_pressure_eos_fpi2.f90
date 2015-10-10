

PROGRAM saha_hydrogen_pressure_eos_fpi2
  USE SahaHydrogenPressureEOSDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_hydrogen_pressure_eos_fpi2.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure_Initial', ' ', 'Pressure_Calc'

  ! CGS
  DO P = lower_P, upper_P, increment_P
  
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0.0

    here_count = 1

    DO

      T = (1.0*P) / (k * (n(1) + N_e))

      g_N_e = const * ( 2.0* B_h(2) / B_h(1) ) * T**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T))

      N_e1 = -1.0 * g_N_e/2.0 + SQRT( g_N_e*n(1) + (g_N_e / 2.0)**2.0)

      PRINT *, const
      PRINT *, (2.0 * B_h(2) / B_h(1) )
      PRINT *, (n(1) - N_e)
      PRINT *, T**(3.0/2.0)
      PRINT *, EXP(-X_h_eV(1)/(k_eV*T))
      
      PRINT *, ' '

      PRINT *, 'P ', P
      PRINT *, 'N_e ', N_e
      PRINT *, 'N_e1 ', N_e1
      PRINT *, 'T ', T
        
      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
      
      IF(ABS((N_e_last - N_e1)/N_e_last) > 1.0e3) THEN
        PRINT *, 'HERE ', ABS((N_e - N_e1)/N_e)
        N_e_last = N_e
        N_e = N_e1
      ELSE
        PRINT *, 'THERE'
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF

    END DO  

    T = (1.0*P) / (k * (n(1) + N_e1))
    
    ! Calculate the Saha Equation for hydrogen 
    ratio_h(1) = const * (T)**(3.0/2.0) * &
        (2.0 * B_h(2) * EXP(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e1)
    ! Calculate the degree of ionization (aka ionization fraction)
    y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1), ' ', P, ' ', Pressure(y_h(1), T)
          
  END DO

END PROGRAM saha_hydrogen_pressure_eos_fpi2


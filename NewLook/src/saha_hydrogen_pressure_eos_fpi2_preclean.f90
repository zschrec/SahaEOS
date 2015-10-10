

PROGRAM saha_hydrogen_pressure_eos_fpi2
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

    here_count = 1

    DO

      T = (1.0*P) / (k * (n(1) + N_e))

      ! Calculate the Saha Equation for hydrogen 
!      IF(n(1) - N_e < 0) THEN
!        N_e1 = 1.0 * ( N_e * const * (2.0 * B_h(2) / B_h(1)) * ABS( n(1) - N_e ) * & 
!          (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)))**(1.0/3.0)
!      ELSE
!        N_e1 = ( N_e * const * (2.0 * B_h(2) / B_h(1)) * ( n(1) - N_e ) * & 
!            (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)))**(1.0/3.0)
!      END IF

!     N_e1 = ( N_e * const * (2.0 * B_h(2) / B_h(1)) * ( n(1) - N_e ) * & 
!          (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)))**(1.0/3.0)

      
      g_N_e = const * ( 2.0* B_h(2) / B_h(1) ) * T**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T))

      !N_e1 = SQRT( const * (2.0 * B_h(2) / B_h(1)) * ( n(1) - N_e ) * & 
      !    (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)))

!      N_e1 = -1.0 * const * ( B_h(2) / B_h(1) ) * T**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)) + &
!        SQRT( n(1) + (const * ( B_h(2) / B_h(1) ) * T**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)))**2.0)

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
        
!      IF ( N_e1 > n(1) ) THEN
!        PRINT *,'POOP'
!        N_e1 = N_e1 / (10.0**((LOG10(N_e1) - LOG10(n(1))) * 2.0))
!      END IF
      
      !N_e1 = n(1) / & 
      !  (1.0/(const*(T)**(3.0/2.0) * (2.0*B_h(2)*EXP(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e)) + &
      !    1.0)

      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
      
      IF(ABS((N_e_last - N_e1)/N_e_last) > 1.0e3) THEN
!      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3 .AND. here_count < 100) THEN
!      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN
        PRINT *, 'HERE ', ABS((N_e - N_e1)/N_e)
        N_e_last = N_e
        N_e = N_e1
        !here_count = here_count + 1
      ELSE
        PRINT *, 'THERE'
        !here_count = 1
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
      END IF




!      here_count = here_count + 1
!      IF (here_count == 100) THEN
!        N_e = N_e * 10.0
!        here_count = 1
!      END IF

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


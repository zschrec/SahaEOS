

PROGRAM saha_hydrogen_pressure_eos_NR
  USE SahaHydrogenPressureEOSDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_hydrogen_pressure_eos_NR.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure_Initial', ' ', 'Pressure_Calc'

  ! CGS
  DO P = lower_P, upper_P, increment_P
  
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0.0

    here_count = 0
    there_count = 0 
    DO

      T = (1.0*P) / (k * (n(1) + N_e))

      IF(T < 100) THEN
        T = 100
      END IF
      ! Calculate the Saha Equation for hydrogen - this is g(Ne)
      ratio_h(1) = const * (T)**(3.0/2.0) * &
          (2.0 * B_h(2) * EXP(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e)
     

      PRINT *, const
      PRINT *, T**(3.0/2.0)
      PRINT *, (2.0 * B_h(2) / B_h(1) )
      PRINT *, EXP(-X_h_eV(1)/(k_eV*T))
      PRINT *, 1.0 / N_e
      
      PRINT *, ' '

      PRINT *, 'P ', P
      PRINT *, 'N_e ', N_e
      PRINT *, 'T ', T
      PRINT *, 'Ratio ', ratio_h(1)      
      
!      IF(ratio_h(1) == 0) THEN
!        ratio_h(1) = 0.5
!      END IF

      G_N_e = n(1) / ((1.0/ratio_h(1)) + 1.0) - N_e
  
!     g_small_prime_N_e = const * (2.0*B_h(2)/B_h(1)) * (T)**(3.0/2.0)*EXP(-X_h_eV(1)/(k_eV*T)) * &
!        ((-3.0/(2.0*N_e*(N_e+n(1)))) - ((-X_h_eV(1) * k)/(P*k_eV*N_e)) - (1.0/N_e**2))

!      G_prime_N_e = n(1) * ((ratio_h(1)**4)/((1+ratio_h(1))**2) * g_small_prime_N_e) - 1
      G_prime_N_e = n(1) * ((ratio_h(1)**5)/((1+ratio_h(1))**2) * &
          (-3.0 / (2.0*(N_e+n(1))) - (-X_h_eV(1)*k)/(P*k_eV) - 1.0/N_e )) - 1
      
      ! Calculate the Saha Equation for hydrogen 
      N_e1 = N_e - (G_N_e / G_prime_N_e) 

      PRINT *, 'G(N_e) ', G_N_e
      PRINT *, 'G_prime(N_e) ', G_prime_N_e
!      PRINT *, 'g_prime(N_e) ', g_small_prime_N_e
      PRINT *, 'N_e1 ', N_e1
 
      IF(N_e1 == 0) THEN
!      IF(N_e1 < 1.0 .AND. (N_e > 1.0e5)) THEN
        PRINT *, 'FISH ', LOG10(N_e), ' '
        !N_e = N_e / (10.0**(LOG10(N_e) - 1.0))
        N_e = 1.0e-10
        PRINT *, 'N_e ', N_e
        CYCLE
      END IF

!      READ *, T

      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
      
!      IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN
      IF(ABS((N_e_last - N_e1)/N_e_last) > 10.0**(-3+LOG10((here_count-there_count+1)/10.0))) THEN
        PRINT *, 'N_e1_Last ', N_e_last, ' ', (N_e_last - N_e1)/N_e_last
        N_e_last = N_e
        N_e = N_e1
        PRINT *, 'HERE'
        here_count = here_count + 1
      ELSE
        PRINT *, 'N_e1_Last ', N_e_last, ' ', (N_e_last - N_e1)/N_e_last
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
        PRINT *, 'THERE'
        there_count = there_count + 1
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

END PROGRAM saha_hydrogen_pressure_eos_NR


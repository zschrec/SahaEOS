

PROGRAM saha_hydrogen_internal_energy_eos_fpi2
  USE SahaHydrogenInternalEnergyEOSDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_hydrogen_internal_energy_eos_fpi2.dat')
  
  !WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'IE_Initial', ' ', 'IE_Calc'

  ! CGS
  DO U = lower_U, upper_U, increment_U
  
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0.0
    
    T_last = 1.0
    T_last_last = 1.0

    there_count = 1

!      IF(U==1100) THEN
!        N_e = 1.043e13
!      END IF

    DO

      T = (2.0*(U - N_e*X_h(1))) / (3.0 * k * (n(1) + N_e))

      IF ( T < 100) THEN
        T = 100
      END IF

      g_N_e = const * ( 2.0* B_h(2) / B_h(1) ) * T**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T))
      N_e1 = -1.0 * g_N_e/2.0 + SQRT( g_N_e*n(1) + (g_N_e / 2.0)**2.0)

      !N_e1 = ( N_e * const * (2.0 * B_h(2) / B_h(1)) * ( n(1) - N_e ) * & 
      !    (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)))**(1.0/3.0)


      !N_e1 = SQRT( const * (2.0 * B_h(2) / B_h(1)) * ( n(1) - N_e ) * & 
      !    (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)))


      PRINT *, const
      PRINT *, (2.0 * B_h(2) / B_h(1) )
      PRINT *, (n(1) - N_e)
      PRINT *, T**(3.0/2.0)
      PRINT *, EXP(-X_h_eV(1)/(k_eV*T))
      
      PRINT *, ' '

      PRINT *, 'U ', U
      PRINT *, 'N_e ', N_e
      PRINT *, 'N_e1 ', N_e1
      PRINT *, 'T ', T
        
      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
     
      IF (U == 11000) THEN
        !READ *,T
        WRITE(10,*) ' ', T, ' ', N_e, ' ', N_e1
      END IF

      IF(ABS((N_e_last - N_e1)/N_e_last) > 1.0e3) THEN
        PRINT *, 'HERE ', ABS((N_e - N_e1)/N_e)
        banana_count = 0
        ham_count = 0
        N_e_last = N_e
        N_e = N_e1
        N_e1_last = N_e1
      ELSE
        PRINT *, 'THERE ', ABS((N_e - N_e1)/N_e), ' ', ABS(T_last_last - T)/T_last_last 
        
        N_e_last = N_e
        N_e = (N_e + N_e1) / 2.0
        
    !    IF(ABS(T_last_last - T)/T_last_last < 1.0e-3) THEN
    !      PRINT *, 'BANANA ', banana_count
    !      banana_count = banana_count + 1
    !      ham_count = 0
!   !!       PRINT *, 'N_e ', N_e
!   !!       PRINT *, 'N_e1 ', N_e1
!   !!       PRINT *, 'N_e_last ', N_e_last
!   !!       PRINT *, 'N_e1_last ', N_e1_last
    !      N_e_last_temp = N_e
    !      N_e1_last_temp = N_e1
    !!      !!!!N_e = (N_e + N_e1) / 2.0

    !      N_e = ((N_e_last-2.0*N_e) * (N_e1-N_e1_last) - &
    !        (N_e-N_e_last) * ((N_e1_last-2.0*N_e1) )) / & 
    !        ((1.0-2.0) * (N_e1-N_e1_last) - &
    !        (N_e-N_e_last) * (1.0-2.0))
          
          !N_e = ((T_last*N_e_last-T*N_e) * (N_e1-N_e1_last) - &
          !  (N_e-N_e_last) * ((T_last*N_e1_last-T*N_e1) )) / & 
          !  ((T_last-T) * (N_e1-N_e1_last) - &
          !  (N_e-N_e_last) * (T_last-T))
          
    !      N_e_last = N_e_last_temp
    !      N_e1_last = N_e1_last_temp

    !    ELSE
          
    !      PRINT *, 'HAM ', ham_count
    !      ham_count = ham_count + 1
    !      banana_count = 0
    !      N_e_last = N_e
    !      N_e1_last = N_e1
    !      N_e = (N_e + N_e1) / 2.0
    !!    END IF
        
        !!IF(U==1100) THEN
        !!  N_e_last = N_e_last_temp 
        !!  N_e = (N_e_last_temp + N_e1_last_temp) / 2.0
        !!END IF 
      END IF
      T_last_last = T_last
      T_last = T
    
    END DO  

    T = (2.0*(U - N_e*X_h(1))) / (3.0* k * (n(1) + N_e1))
    
    ! Calculate the Saha Equation for hydrogen 
    ratio_h(1) = const * (T)**(3.0/2.0) * &
        (2.0 * B_h(2) * EXP(-X_h_eV(1)/(k_eV*T))) / (B_h(1) * N_e1)
    ! Calculate the degree of ionization (aka ionization fraction)
    y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    !WRITE(10,*) T,' ', y_h(1), ' ', U, ' ', InternalEnergy(y_h(1), T)
          
  END DO

END PROGRAM saha_hydrogen_internal_energy_eos_fpi2


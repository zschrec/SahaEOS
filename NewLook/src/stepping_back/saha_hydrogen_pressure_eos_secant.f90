

PROGRAM saha_hydrogen_pressure_eos_secant
  USE SahaHydrogenPressureEOSDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_hydrogen_pressure_eos_secant.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure_Initial', ' ', 'Pressure_Calc'

  ! CGS
  DO P = lower_P, upper_P, increment_P
  
    N_e = N_e_initial
    !N_e_last = 0.0
    N_e_last = 1.0e5
    N_e1 = 0.0

      
    PRINT *, 'P ', P
    
    DO

      T_N_e = (1.0*P) / (k * (n(1) + N_e))
      T_N_e_last = (1.0*P) / (k * (n(1) + N_e_last))

!      ! Calculate the Saha Equation for hydrogen 
!      ratio_h(1) = const * (2.0 * B_h(2) / B_h(1)) * &
!          (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)) * (1.0/N_e)

!      g_N_e = const * (2.0 * B_h(2) / B_h(1)) * &
!          (T_N_e)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T_N_e)) * (1.0/N_e)
      
!      g_N_e_last = const * (2.0 * B_h(2) / B_h(1)) * &
!          (T_N_e_last)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T_N_e_last)) * (1.0/N_e_last)

      g_N_e = Little_g(N_e, T_N_e)
      g_N_e_last = Little_g(N_e_last, T_N_e_last)


!      ! Calculate the degree of ionization (aka ionization fraction)
!      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)
      
!      GG_N_e = N_e - (n(1) / ((1.0/g_N_e) + 1.0))
!      GG_N_e_last = N_e_last - (n(1) / ((1.0/g_N_e_last) + 1.0))

      GG_N_e = Big_G(N_e, g_N_e)
      GG_N_e_last = Big_G(N_e_last, g_N_e_last)

  !    PRINT *, 'P ', P
  !    PRINT *, 'N_e ', N_e
  !    PRINT *, 'N_e_last ', N_e_last
  !    PRINT *, 'T_N_e ', T_N_e
  !    PRINT *, 'T_N_e_last ', T_N_e_last
  !    PRINT *, 'g_N_e ', g_N_e
  !    PRINT *, 'g_N_e_last ', g_N_e_last
  !    PRINT *, 'GG_N_e ', GG_N_e
  !    PRINT *, 'GG_N_e_last ', GG_N_e_last
            
!      ! Calculate average electron contribution 
!      v_e(1) = 1.0*y_h(1)
            
!      ! Recompute electron density
!      N_e1 = n(1) * v_e(1) 

      N_e1 = N_e - ( (G_N_e * (N_e - N_e_last)) / (G_N_e - G_N_e_last) )
     
  !    PRINT *, 'SECANT ' , ( (G_N_e * (N_e - N_e_last)) / (G_N_e - G_N_e_last) )
  !    PRINT *, 'G_SUB ', G_N_e - G_N_e_last
  !    PRINT *, 'N_e1 ', N_e1
  !    PRINT *, 'THIS ', (N_e_last - N_e1)/N_e_last, ' ', (N_e - N_e1)/N_e

      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT
      
  !    IF(ABS((N_e_last - N_e1)/N_e_last) > 1e-3) THEN
  !      PRINT *, 'HERE ', (N_e_last - N_e1)/N_e_last, ' ', (N_e - N_e1)/N_e
        N_e_last = N_e
        N_e = N_e1
  !    ELSE
  !      PRINT *, 'THERE ', (N_e_last - N_e1)/N_e_last, ' ', (N_e - N_e1)/N_e
  !      N_e_last = N_e
  !      N_e = (N_e + N_e1) / 2.0
  !    END IF
 
    END DO  

    N_e = N_e1
    T_N_e = (1.0*P) / (k * (n(1) + N_e))
    ratio_h(1) = const * (2.0 * B_h(2) / B_h(1)) * &
        (T_N_e)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T_N_e)) * (1.0/N_e)
    y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    PRINT *, N_e
    WRITE(10,*) T_N_e,' ', y_h(1), ' ', P, ' ', Pressure(y_h(1), T_N_e), ' ', N_e
          
  END DO

END PROGRAM saha_hydrogen_pressure_eos_secant


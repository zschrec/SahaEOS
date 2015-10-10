

PROGRAM saha_hydrogen_pressure_eos_fpi
  USE SahaHydrogenPressureEOSDeclarations
  IMPLICIT NONE
  
  N_e_initial = N_e
  N_e_last = N_e

  OPEN(10, FILE='outfiles/saha_hydrogen_pressure_eos_fpi.dat')
  WRITE(10,*) 'T', ' ', 'n_h1/n_h', ' ', 'Pressure_Initial', ' ', 'Pressure_Calc'

  ! CGS
  DO P = lower_P, upper_P, increment_P
  
    N_e = N_e_initial
    N_e_last = N_e_initial
    N_e1 = 0.0

    DO

      T = (1.0*P) / (k * (n(1) + N_e))

      IF (T < 0) THEN
        READ *, g_Ne
      END IF
!      ! Calculate the Saha Equation for hydrogen 
!      ratio_h(1) = const * (2.0 * B_h(2) / B_h(1)) * &
!          (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)) * (1.0/N_e)

      g_Ne = const * (2.0 * B_h(2) / B_h(1)) * &
          (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T))

!      ! Calculate the degree of ionization (aka ionization fraction)
!      y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

!      ! Calculate average electron contribution 
!      v_e(1) = 1.0*y_h(1)
            
!      ! Recompute electron density
!      N_e1 = n(1) * v_e(1) 

      PRINT *, '***************************'
      PRINT *, 'n(1) ', n(1)
      PRINT *, 'P ', P
      PRINT *, 'N_e ', N_e
      PRINT *, 'g_Ne ', g_Ne
      PRINT *, 'T ', T
      PRINT *, '***************************'
      
      ! ATTEMPT NUMBER 1
      ! N_e1 = n(1) - ( (N_e**2.0) / g_Ne ) 

      ! ATTEMPT NUMBER 2
      ! N_e1 = SQRT( g_Ne * (n(1) - N_e) ) 

      ! ATTEMPT NUMBER 3x
      ! x_ex = 5.0
      ! IF (N_e < 0) THEN
      !   PRINT *, 'A'
      !   PRINT *, N_e
      !   N_e1 = -1.0 * ( (-1.0*N_e)**x_ex * g_Ne * (n(1) - (-1.0*N_e)) )**(1.0/(x_ex+2.0))
      ! ELSE IF (g_Ne < 0) THEN 
      !   PRINT *, 'B'
      !   PRINT *, g_Ne
      !   N_e1 = -1.0 * ( N_e**x_ex * (-1.0*g_Ne) * ABS(n(1) - N_e) )**(1.0/(x_ex+2.0))
      ! ELSE IF ((n(1) - N_e)  < 0)  THEN
      !   PRINT *, 'C'
      !   PRINT *, (n(1) - N_e)
      !   N_e1 = ( N_e**x_ex * g_Ne * (-1.0*(n(1) - N_e)) )**(1.0/(x_ex+2.0))
      ! ELSE
      !   N_e1 = ( N_e**x_ex * g_Ne * (n(1) - N_e) )**(1.0/(x_ex+2.0))
      ! END IF

      ! ATTEMPT NUMBER 4
      N_e1 = g_Ne/(-2.0) + SQRT( g_Ne*n(1) + (g_Ne/2.0)**2.0)

      IF (ABS((N_e - N_e1)/N_e) < tol) EXIT


!      N_e = N_e1

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

    T = (1.0*P) / (k * (n(1) + N_e1))

    ! Calculate the Saha Equation for hydrogen 
    ratio_h(1) = const * (2.0 * B_h(2) / B_h(1)) * &
        (T)**(3.0/2.0) * EXP(-X_h_eV(1)/(k_eV*T)) * (1.0/N_e1)

    ! Calculate the degree of ionization (aka ionization fraction)
    y_h(1) = 1.0 / ((1.0/ratio_h(1)) + 1.0)

    ! Log temperature, ionization fraction, pressure, internal energy, and specific heat
    WRITE(10,*) T,' ', y_h(1), ' ', P, ' ', Pressure(y_h(1), T)
          
  END DO

END PROGRAM saha_hydrogen_pressure_eos_fpi

